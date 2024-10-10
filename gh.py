import sys
import json
from pathlib import Path
from subprocess import run, PIPE

_GIT_REPO = False
_GH_INSTALLED = False
_POETRY_INSTALLED = False


def assert_git_repo():
    global _GIT_REPO
    if not _GIT_REPO:
        cmd = run(["git", "status"], capture_output=True)
        if cmd.returncode != 0:
            raise GitError(f"{Path.cwd()} is not a git repository.")
        else:
            _GIT_REPO = True

def assert_gh_installed():
    global _GH_INSTALLED
    if not _GH_INSTALLED:
        cmd = run(["gh", "--version"], capture_output=True)
        if cmd.returncode != 0:
            raise GhError("gh is not installed.")
        else:
            _GH_INSTALLED = True

def assert_poetry_installed():
    global _POETRY_INSTALLED
    if not _POETRY_INSTALLED:
        cmd = run(["poetry", "--version"], capture_output=True)
        if cmd.returncode != 0:
            raise PoetryError("poetry is not installed.")
        else:
            _POETRY_INSTALLED = True


def git_assert_working_tree_clean():
    assert_git_repo()
    cmd = run(["git", "diff", "--quiet"], capture_output=True)
    if cmd.returncode != 0:
        raise GitError("Some changes are not staged. Stage and commit first.")
    cmd = run(["git", "diff", "--staged", "--quiet"], capture_output=True)
    if cmd.returncode != 0:
        raise GitError("Some changes are staged. Commit first.")

def _run_git(cmds):
    assert_git_repo()
    cmd = run(["git", *cmds], capture_output=True)
    if cmd.returncode == 0:
        return cmd.stdout.decode('UTF-8').strip()
    else:
        raise GitError(' '.join(cmds) + ":\n" + cmd.stderr.decode('UTF-8').strip())

def git_current_branch():
    return _run_git(["symbolic-ref", "--short", "HEAD"])

def git_rename_current_branch(new_branch, set_upstream=False):
    old_branch = git_current_branch()
    if old_branch == 'main':
        raise GitError("Cannot rename the main branch.")
    out = _run_git(["branch", "-m", new_branch])
    if out: print(out)
    if set_upstream:
        out = _run_git(["push", "origin", "--delete", old_branch])
        if out: print(out)
        git_push(set_upstream=True)

def git_switch(branch, create=False):
    cc = ['-c'] if create else []
    out = _run_git(["switch", *cc, branch])
    if out: print(out)

def git_add(files):
    out = _run_git(["add", *files])
    if out: print(out)

def git_commit(message, no_verify=False):
    vv = ['--no-verify'] if no_verify else []
    out = _run_git(["commit", *vv, '-m', message])
    if out: print(out)

def git_pull():
    out = _run_git(["pull"])
    if out: print(out)

def git_push(set_upstream=False):
    uu = ["--set-upstream", "origin", git_current_branch()] if set_upstream else []
    out = _run_git(["push", *uu])
    if out: print(out)

def git_make_tag(tag):
    out = _run_git(["tag", tag])
    if out: print(out)
    out = _run_git(["push", "origin", tag])
    if out: print(out)


def _run_gh(cmds):
    assert_gh_installed()
    cmd = run(["gh", *cmds], capture_output=True)
    if cmd.returncode == 0:
        return cmd.stdout.decode('UTF-8').strip()
    else:
        raise GhError(' '.join(cmds) + ":\n" + cmd.stderr.decode('UTF-8').strip())

def gh_pr_create(base_branch, title):
    out = _run_gh(["pr", "create", "--base", base_branch, "--title", title, '--fill'])
    if out: print(out)

def gh_pr_list(base=None, head=None):
    bb = [] if base is None else ["-B", base]
    hh = [] if head is None else ["-H", head]
    out = _run_gh(["pr", "list", *bb, *hh, "--json",
                   "number,author,headRepositoryOwner,headRepository,headRefName"])
    data = json.loads(out)
    return {int(pr['number']):
            f"{pr['headRepositoryOwner']['login']}:{pr['headRepository']['name']}/{pr['headRefName']}"
            for pr in data}

def gh_pr_merge(pr_id, admin=False, delete_branch=False):
    adm = ["--admin"] if admin else []
    db = ["--delete-branch"] if delete_branch else []
    out = _run_gh(["pr", "merge", str(pr_id), "--merge", *adm, *db])
    if out: print(out)

def gh_release_create(version, title, draft=False):
    dr = ["--draft"] if draft else []
    out = _run_gh(["release", "create", version, *dr, "--generate-notes",
                   "--title", title, "--verify-tag"])
    if out: print(out)


def _run_poetry(cmds):
    assert_poetry_installed()
    cmd = run(["poetry", *cmds], capture_output=True)
    if cmd.returncode == 0:
        return cmd.stdout.decode('UTF-8').strip()
    else:
        raise PoetryError(' '.join(cmds) + ":\n" + cmd.stderr.decode('UTF-8').strip())

def poetry_bump_version(bump):
    out = _run_poetry(["version", bump])
    if out: print(out)

def poetry_get_version():
    return _run_poetry(["version"]).split()[-1]

def poetry_get_expected_version(bump):
    return _run_poetry(["version", bump, "--dry-run"]).split()[-1]

def poetry_publish(build=False):
    pp = ["--build"] if build else []
    out = _run_poetry(["publish", *pp])
    if out: print(out)


class GitError(OSError):
    pass

class GhError(OSError):
    pass

class PoetryError(OSError):
    pass
