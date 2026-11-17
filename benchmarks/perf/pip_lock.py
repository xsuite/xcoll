"""Build a portable pip-only lock for paired external-engine benchmarks.

The Conda explicit lock recreates all Conda-managed distributions. Standard
``pip freeze`` includes direct references into temporary Conda feedstock build
directories, which cannot be replayed in a second environment. Exclude packages
already provided by the *same version* in the original Conda environment,
while retaining pip additions and pip upgrades at their exact pinned versions.

Run this only using the trusted main workflow's code, in the fully provisioned
baseline environment, with a Conda inventory recorded before pip installs.
"""
from __future__ import annotations

import argparse
from importlib import metadata
import json
from pathlib import Path
import re
import subprocess
import sys


_REQUIREMENT = re.compile(r"^([A-Za-z0-9][A-Za-z0-9_.-]*)\s*(?:===?|@\s+)")
_LOCAL_URL = re.compile(r"(?:^|\s)@\s*file://", re.IGNORECASE)


def normalise(name):
    return re.sub(r"[-_.]+", "-", name).lower()


def portable_lock(conda_inventory, frozen, *, installed_version=None):
    """Return requirements not already provided by the identical Conda lock.

    Fail closed on local/editable pip sources, including sources that would only
    exist on the baseline machine. Do not silently substitute PyPI packages for
    Conda builds or erase pip-installed upgrades to Conda packages.
    """
    if installed_version is None:
        installed_version = metadata.version
    conda = {}
    for entry in conda_inventory:
        if not isinstance(entry, dict) or not isinstance(entry.get('name'), str):
            raise ValueError('Malformed Conda package inventory')
        if not isinstance(entry.get('version'), str):
            raise ValueError('Malformed Conda package version')
        if str(entry.get('channel', '')).lower() != 'pypi':
            conda[normalise(entry['name'])] = entry['version']
    if not conda:
        raise ValueError('Empty Conda package inventory')

    locked = []
    seen = set()
    for raw in frozen.splitlines():
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        match = _REQUIREMENT.match(line)
        if not match:
            raise ValueError(f'Unsupported pip freeze entry: {line}')
        name = normalise(match.group(1))
        if name == 'xcoll':
            continue  # Candidate Xcoll must always come from its own checkout.
        if name in seen:
            raise ValueError(f'Duplicate pip freeze entry: {name}')
        seen.add(name)
        if name in conda and installed_version(match.group(1)) == conda[name]:
            continue  # This exact distribution already exists via Conda.
        if _LOCAL_URL.search(line):
            raise ValueError(f'Non-portable local pip dependency: {line}')
        locked.append(line)
    if not locked:
        raise ValueError('Empty pip-only lock; baseline pip installation missing?')
    return ''.join(line + '\n' for line in locked)


def portable_constraints(packages):
    """Pin all installed distributions (including Conda's) without reinstalling them."""
    versions = {}
    for pkg in packages:
        if not isinstance(pkg, dict):
            raise ValueError('Malformed pip inventory entry')
        name, version = pkg.get('name'), pkg.get('version')
        if (not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', name)
                or not isinstance(version, str)
                or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9.+!_-]*', version)):
            raise ValueError('Malformed pip inventory package or version')
        key = normalise(name)
        if key == 'xcoll':
            continue
        if key in versions:
            raise ValueError(f'Duplicate pip inventory package: {key}')
        versions[key] = version
    if not versions:
        raise ValueError('Empty pip inventory')
    return ''.join(f'{name}=={version}\n' for name, version in sorted(versions.items()))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--conda-inventory', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--constraints', required=True, type=Path)
    args = p.parse_args(argv)
    inventory = json.loads(args.conda_inventory.read_text())
    if not isinstance(inventory, list):
        raise ValueError('Expected a JSON list of Conda packages')
    frozen = subprocess.check_output(
        [sys.executable, '-m', 'pip', 'freeze', '--exclude', 'xcoll'], text=True)
    lock = portable_lock(inventory, frozen)
    # Version-only constraints prevent pip from upgrading an exact Conda
    # package build while installing candidate's additional requirements.
    all_packages = json.loads(subprocess.check_output(
        [sys.executable, '-m', 'pip', 'list', '--format=json'], text=True))
    if not isinstance(all_packages, list):
        raise ValueError('Expected a JSON list of installed distributions')
    constraints = portable_constraints(all_packages)
    args.output.write_text(lock)
    args.constraints.write_text(constraints)
    print(f'Locked {len(lock.splitlines())} pip distributions not supplied by Conda '
          f'and constrained {len(constraints.splitlines())} installed packages')


if __name__ == '__main__':
    main()
