"""Stable, human-readable local machine identity for benchmark metadata.

Prefer the operating system's *computer name*, not the Bonjour/DNS hostname.
Override explicitly with BENCHMARK_MACHINE_NAME for headless or renamed systems.
"""
from __future__ import annotations

import os
import platform
import subprocess
from pathlib import Path


def _command(*argv: str) -> str | None:
    try:
        result = subprocess.run(argv, text=True, capture_output=True, check=False,
                                timeout=2)
    except (OSError, subprocess.TimeoutExpired):
        return None
    value = result.stdout.strip() if result.returncode == 0 else ''
    return value or None


def machine_name() -> str:
    """Return a stable name without using network-specific host discovery."""
    override = os.environ.get('BENCHMARK_MACHINE_NAME', '').strip()
    if override:
        return override

    system = platform.system()
    if system == 'Darwin':
        # scutil ComputerName is the friendly name (e.g. BlackMetal3).
        # HostName and LocalHostName are deliberately not used.
        value = _command('scutil', '--get', 'ComputerName')
        if value:
            return value
    elif system == 'Windows':
        value = os.environ.get('COMPUTERNAME', '').strip()
        if value:
            return value
    elif system == 'Linux':
        value = _command('hostnamectl', '--static')
        if value:
            return value
        try:
            value = Path('/etc/hostname').read_text(encoding='utf-8').strip()
        except OSError:
            value = ''
        if value:
            return value

    # An honest fallback: do not pretend a Bonjour/network hostname is stable.
    return f'{system} | {platform.machine()}'
