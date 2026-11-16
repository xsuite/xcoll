# Supplementary pull-request validation

This extends the existing CPU PR workflow without changing its trigger or
concurrency group. It uses the **existing** CERN runner; no second VM is needed.

## 1. Manually rerun only Everest physics

On GitHub Actions, select **Xcoll physics (PR)** → **Run workflow**, choose
`main`, and enter the PR number. Or run:

```sh
gh workflow run xcoll-physics-pr.yml --repo xsuite/xcoll --ref main -f pr_number=123
```

Manual physics validates the PR's existing GitHub test merge, without asking
GitHub to refresh it or modifying any branch. The resolver verifies that the
merge's second parent is the current PR head and records its actual first
parent as the tested base. It fails immediately if the merge ref is missing,
conflicted, or belongs to an older PR head. Automatic PR runs retain their
original event SHAs. Independently, the resolver pins a reviewed main checkout
(the dispatch SHA for manual runs) for the accepted reference and validation
harness. Its committed manifest is mandatory; CI never regenerates a reference.
The reference source SHA, tooling SHA, and tested base SHA are distinct.
The candidate's editable Xcoll package and lattice data are used for tracking;
the trusted harness supplies the collector, comparison, and snapshot decoder.
A candidate with an incompatible physics API or changed lattice data can fail
rather than silently substituting main's physics. The external-engine workflow
continues to require a merge with current main and its reviewed head.
The three physics matrix jobs use the pinned pair; the CPU benchmark workflow
is not triggered or cancelled. A new physics dispatch may cancel an older
physics run for the same PR. Use a fresh dispatch (or rerun all jobs), not a
failed-jobs-only retry, so provenance and report artifacts belong to one attempt.

The trusted hosted publisher archives successful, failed, and incomplete
physics runs. It reads bounded JUnit data, not candidate code. Missing matrix
reports are NOT_RUN, never silently counted as passing. The dashboard includes
all 40 physical cases, their statuses and durations, the reference-source SHA,
and the screening value/limit/status for each available captured comparison.
The helper tests still contribute to the overall workflow conclusion.

A manual dispatch is a new check suite on main. When the tested merge is still
current and its base is current main, the publisher also attempts to post a distinct **Xcoll physics
(manual)** commit status on that merge. This does not replace old failed
matrix checks or satisfy a different required-check name. If a required
automatic check is still blocking merging, use the normal branch-update/PR
trigger after any long benchmark has finished, or follow your branch-protection
policy. Never post a success for a different base or head revision.

## 2. One-time administrator setup for external PR benchmarks

**Do not enable the workflow until these controls are in place.**

1. In the organisation runner group `xcoll-benchmarks`, retain the existing
   selected-workflow restrictions and add exactly:

   ```text
   xsuite/xcoll/.github/workflows/xcoll-performance-external-pr.yml@refs/heads/main
   ```

   Do not allow arbitrary workflows or PR refs to use the group. If the plan
   does not support selected-workflow restrictions, do not enable this public
   repository's persistent-runner PR execution.
2. In Xcoll repository **Settings → Environments**, create `xcoll-external-pr`.
   Restrict deployment branches to `main`, configure required reviewers from
   the trusted administrators, and disable bypass of protection rules if your
   settings support it. Do not add any secrets to this environment. If there is
   only one administrator available for approval, leave self-review allowed;
   otherwise require a second reviewer. GitHub's standard environment reviewer
   configuration requires only one of the configured reviewers to approve.
3. In repository **Settings → Secrets and variables → Actions → Variables**,
   create `XCOLL_EXTERNAL_PR_ENABLED` with value `true`. This is the explicit
   opt-in; unset/false means no external PR job can be scheduled.
4. Keep the current FLUKA paths, engine/version variables, and Miniforge setup
   from `RUNNERS.md`. The new workflow uses the same runner group, labels,
   and concurrency queue as the daily jobs.

The hosted resolver checks repository `admin` permission via GitHub's API for
both the original and triggering actor; write/maintain permission is not
sufficient. Any API error fails closed. GitHub may still let a non-admin with
write permission submit a workflow_dispatch request, but the resolver rejects
it before the self-hosted job. The external job also requires authorisation
for the current attempt, so a failed-jobs-only rerun cannot reuse an old gate.

### Security boundary

Approval is authorisation, **not sandboxing**. The exact reviewed contributor
code (including Python packaging hooks, benchmark setup, compiled code, and
its dependency declarations) executes as the existing runner user. Read-only
GitHub tokens, disposable Conda environments, and removal of state directories
are good hygiene but cannot stop hostile code from reading accessible CERN
files, reaching the network, or persisting on the VM. Only approve a revision
you trust to execute on the machine. Review its entire diff and relevant
upstream/dependency changes, not only `.github/workflows`. Do not store personal
credentials, sudo access, Kerberos tokens, or other sensitive credentials on
that account. Follow CERN policies for access to licensed FLUKA software and
network resources. If suspicious code must be tested, do not use this workflow;
manual approval alone is not adequate isolation for untrusted execution.

## 3. Dispatch paired FLUKA or Geant4 measurements

Choose the **Xcoll external performance (PR)** workflow from the Actions tab.
Enter the open PR number, the full 40-character **reviewed head SHA**, and
`fluka` or `geant4`. Select main as the workflow branch. Alternatively:

```sh
# Inspect the PR and review the exact head before copying the full SHA.
gh pr view 123 --repo xsuite/xcoll --json headRefOid,baseRefName,url

gh workflow run xcoll-performance-external-pr.yml --repo xsuite/xcoll --ref main \
  -f pr_number=123 -f engine=fluka \
  -f approved_head_sha=REPLACE_WITH_REVIEWED_40_CHARACTER_SHA
```

Repeat with `engine=geant4` if both engines are needed. Approve the environment
job once the resolver has reported the exact base/head/merge pair. If the PR
head changes before resolution, the run fails; if it changes afterwards, only
the already pinned revision can be checked out. A new revision needs a new
review and dispatch. Use a fresh dispatch for retries, not rerun-failed-jobs.

Both revisions run sequentially on the same CERN VM, with separate temporary
Conda environments and XCOLL_PATH directories. The candidate reuses the base
Conda explicit specification and pip dependency lock (except Xcoll); its
interface and kernels are rebuilt. The engine-specific benchmark mark selects
FLUKA or Geant4. The job shares `xcoll-cern-benchmark-cpu` with CPU, FLUKA, and
Geant4 scheduled benchmarks, with queueing rather than cancellation. The
existing GitHub-hosted CPU PR workflow is not affected.

These are **paired external-engine comparisons**, not comparisons against a
scheduled CPU baseline. The original `pr_compare` tool computes the results;
only a successful nonempty comparison is archived. Failed external runs are
available in Actions with diagnostic artifacts, not treated as successful
measurements. No write token or history publication executes on CERN. The
separate hosted `workflow_run` publisher validates the data and writes it.

## 4. Dashboard and history

The new publisher appends independent immutable JSON records to the existing
`benchmark-results:pr-runs/<PR>/<run-id>-<attempt>.json` directory, with the
same history-writer queue as the CPU publisher. It reuses the existing Pages
deployment workflow and its queue. Existing CPU PR records remain readable
without conversion. Choose **Pull Requests**, the PR number, and the labelled
CPU / FLUKA / Geant4 / Everest physics workflow run. Related runs for identical
head, base, and merge SHAs link to one another. If main moved between the CPU
benchmark and the manually dispatched physics run, the two records remain
separate; the dashboard can show the physics result alongside the benchmark with an explicit different-revision warning, never as the same tested pair.

The new publication schema for physics is separate from the existing benchmark
rows. It includes the PR, branch, run ID, attempt, tested revisions, reference
SHA, trusted tooling SHA, overall conclusion, and per-case statuses/checks. It excludes raw particle
snapshots, full captured logs, and full XML reports. The original artifacts
remain available in Actions for investigation. External benchmark records keep
the existing schema with a profile label, so all timing/physics detail rows
remain usable in the current dashboard.
