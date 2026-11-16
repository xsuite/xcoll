# Xcoll performance infrastructure: setup and operations

This runbook describes the workflows in `xsuite/xcoll` on `main` (reviewed at `cbf03dc7`). It covers the CERN reference machine, GitHub-hosted PR and publication jobs, external engines, history, Pages, and local development. Repository and organisation settings are not stored in Git: verify the settings below when rebuilding the service. The VM flavour, hostname, project, and FLUKA installation paths are site-specific; the examples below are recommendations, not a record of the current VM's exact settings.

For independent manual physics validation, PR dashboard results, and administrator-approved FLUKA/Geant4 PR runs on the same CERN VM, see [PR_WORKFLOWS.md](PR_WORKFLOWS.md). Complete its opt-in and runner-group setup before enabling external PR execution.

## 1. Workflow architecture

| Workflow (`.github/workflows/`) | Trigger | Runner |
| --- | --- | --- |
| `xcoll-performance.yml` | Daily 03:17 UTC; manual quick/full | CERN for main; `ubuntu-24.04` for other target refs |
| `xcoll-performance-fluka.yml` | Daily 05:37 UTC; manual | CERN |
| `xcoll-performance-geant4.yml` | Daily 08:57 UTC; manual | CERN |
| `xcoll-performance-pr.yml` | Non-draft PR to main, excluding `release/` | `ubuntu-24.04`, paired base and merge revision |
| `xcoll-performance-publish.yml` | Reusable workflow called by FLUKA/Geant4 | `ubuntu-24.04` |
| `xcoll-performance-pr-publish.yml` | Successful PR benchmark `workflow_run` | `ubuntu-24.04` |
| `xcoll-performance-dashboard.yml` | Manual redeployment | `ubuntu-24.04` |

The reference workflow also has its own publisher/deployer. All CERN jobs share `xcoll-cern-benchmark-cpu` (queue, no cancellation); history writers share `xcoll-history-writer`, and Pages deployments share `xcoll-performance-pages`. Keep one persistent reference runner; do not consume the general Xsuite runner pool. The CERN template in `benchmarks/templates/` is for reuse, not a second active schedule.

## 2. Create the CERN OpenStack VM

1. Sign in to [CERN OpenStack](https://openstack.cern.ch/) with an authorised CERN account. Select the project with the appropriate quota and ownership; record its name, VM name, flavour, image, and disk size in the private operations inventory.
2. Import an SSH public key (keep the private key on the administrator's machine). In Compute / Instances, launch an `ALMA9 - x86_64` instance; see [CERN Linux images](https://linux.web.cern.ch/openstack/). Select a stable x86-64 flavour and adequate persistent disk. As a starting recommendation, 4–8 vCPUs, 16 GiB RAM, and 80–100 GiB disk are reasonable; adjust for FLUKA/Geant4 and compiler memory. These are not the verified existing flavour. A previous run ran out of memory, so check peak memory and avoid an undersized VM.
3. Attach the appropriate CERN network and security group. Permit inbound SSH only from approved administration networks; the runner needs outbound HTTPS, DNS, and package-repository connectivity, not inbound connections from GitHub. Follow CERN networking and registration policies.
4. SSH using the image's configured administrator account and the key. Apply OS updates, install prerequisites, and create a dedicated service user. Do not grant the runner sudo or copy personal credentials into its home.

```bash
sudo dnf upgrade -y
sudo dnf install -y git curl tar gzip unzip make gcc gcc-c++ gcc-gfortran \
  libicu openssl-libs zlib policycoreutils-python-utils
sudo useradd --create-home --shell /bin/bash xcollrunner
sudo -iu xcollrunner
```

Check `uname -m`, `lscpu`, `free -h`, `df -h`, and `curl -I https://github.com`. The workflow creates compiler environments itself; system `gfortran` is also required by FLUKA's `rfluka`. Use the GitHub runner's `bin/installdependencies.sh` as administrator if the downloaded runner reports missing libraries. Keep SELinux enabled.

## 3. Install Miniforge

As `xcollrunner`:

```bash
cd ~
curl -fsSL -o miniforge.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash miniforge.sh -b -p "$HOME/miniforge3"
rm miniforge.sh
~/miniforge3/bin/conda --version
~/miniforge3/bin/mamba --version
```

The workflows expect `/home/xcollrunner/miniforge3`; change `CONDA_EXE`/`MINIFORGE_ROOT` in the workflows if necessary. No `conda init` or system Python is needed. CPU and FLUKA use Python 3.12; Geant4 uses Python 3.14. The CPU workflow pins Conda GCC/G++ 16.2; external engines have their own toolchains. Each run installs released Xsuite and the selected Xcoll checkout into a temporary environment, prebuilds serial kernels, and cleans up afterwards. Keep adequate space in `$RUNNER_TEMP` and the Conda package cache.

## 4. Configure the GitHub organisation runner group

As an `xsuite` organisation owner, go to organisation **Settings → Actions → Runner groups**. Create `xcoll-benchmarks`, allow only the `xsuite/xcoll` repository, and explicitly permit public-repository access if required. Restrict the group to selected workflows pinned to main:

```text
xsuite/xcoll/.github/workflows/xcoll-performance.yml@refs/heads/main
xsuite/xcoll/.github/workflows/xcoll-performance-fluka.yml@refs/heads/main
xsuite/xcoll/.github/workflows/xcoll-performance-geant4.yml@refs/heads/main
```

See [GitHub runner group access](https://docs.github.com/en/actions/how-tos/manage-runners/self-hosted-runners/manage-access). Verify that the organisation plan offers workflow restrictions. Repository-only restrictions and labels are not a security boundary for a public repository. If the selected-workflow restriction is unavailable, do not treat the current setup as safely isolated; use a private orchestration repository and change the workflow design before exposing a persistent CERN runner. Keep reviewed main protected, particularly `.github/workflows/` and benchmark tooling.

In organisation **Settings → Actions → Runners → New self-hosted runner**, choose Linux x64. As `xcollrunner`, create `~/actions-runner`, download/extract the version GitHub displays, and verify its published SHA256. Register using the short-lived registration token supplied by GitHub (not a PAT or `GITHUB_TOKEN`):

```bash
cd ~/actions-runner
./config.sh --url https://github.com/xsuite \
  --token '<SHORT-LIVED-REGISTRATION-TOKEN>' \
  --name xcoll-cern-benchmark \
  --runnergroup xcoll-benchmarks \
  --labels xcoll-bench-cpu \
  --unattended
```

The runner also receives default `self-hosted`, `linux`, and `x64` labels. Do not add `--no-default-labels`. The workflow selector is:

```yaml
runs-on:
  group: xcoll-benchmarks
  labels: [self-hosted, linux, x64, xcoll-bench-cpu]
```

If an existing runner was registered at repository scope, remove and re-register it at organisation scope, or use the appropriate organisation runner migration process; a repository runner cannot simply be assigned to an organisation runner group. For an organisation runner already in Default, move it into the new group.

As administrator:

```bash
cd /home/xcollrunner/actions-runner
sudo ./svc.sh install xcollrunner
sudo ./svc.sh start
sudo ./svc.sh status
```

Verify the runner is Idle/Online in the organisation and visible to Xcoll. Check `journalctl -u <runner-service-name>` for errors. On SELinux `203/EXEC`, inspect `sudo ausearch -m AVC -ts recent` and `ls -Z runsvc.sh`; use a narrowly scoped persistent `semanage fcontext` rule and `restorecon` for the denied executable if the AVC confirms a labelling issue. For example, a policy-approved executable type such as `bin_t` on the exact `runsvc.sh` path may be appropriate; confirm with the administrator rather than relabelling all of `/home` or disabling SELinux. Restart and verify after reboot. Never store the registration token in the repository or runner logs.

## 5. GitHub-hosted runners and Actions settings

The `ubuntu-24.04` jobs are provisioned automatically; no runner registration or VM setup is needed. Under repository **Settings → Actions → General**, enable Actions and permit the actions used in the workflow files (`actions/checkout`, `setup-python`, `upload-artifact`, `download-artifact`, `configure-pages`, `upload-pages-artifact`, `deploy-pages`). Set default workflow permissions to read-only. Allow the explicitly declared publication/Pages permissions in the organisation policy. Do not enable write tokens or secrets for fork PRs; PR measurement jobs are read-only and must stay on hosted machines. The trusted `workflow_run` publisher downloads PR results as data, not executable PR code.

The current workflows use the built-in, job-scoped `GITHUB_TOKEN` (`github.token`); no separate PAT, deploy key, or repository secret is needed. Benchmark jobs declare `contents: read` and disable checkout credential persistence. History publishers declare `contents: write` and `actions: read`; Pages jobs declare `contents: read`, `actions: read`, `pages: write`, and `id-token: write`. The workflow's `git push` uses the publisher checkout credentials. `GH_TOKEN` for artifact checks is set from `github.token`. The Pages environment is `github-pages`. Do not confuse the runner registration token with these per-job tokens.

## 6. External engines

FLUKA needs an authorised, working local installation accessible to the runner account. Provision it according to the FLUKA licence and CERN policy (including any server/FlukaIO and LineBuilder requirements); do not assume EOS is mounted or accessible, and do not put personal CERN credentials on the runner. Prefer an approved read-only local installation. Use repository **Settings → Secrets and variables → Actions → Variables** for these paths, not secrets:

| Variable | Value |
| --- | --- |
| `XCOLL_FLUKA_ROOT` | Root containing `bin/rfluka` |
| `XCOLL_FLUKA_SERVER` | Full path to the FLUKA server executable |
| `XCOLL_FLUKAIO_LIB` | Path to the compiled FlukaIO library |
| `XCOLL_LINEBUILDER_ROOT` | LineBuilder installation root |
| `XCOLL_FLUKA_FEDB` | Optional FEDB path |
| `BENCHMARK_FLUKA_VERSION` | Optional accurate version override |
| `BENCHMARK_LINEBUILDER_VERSION` | Optional accurate version override |
| `BENCHMARK_FLUKAIO_VERSION` | Optional accurate version override |

As `xcollrunner`, verify each path is readable, executables are executable, shared libraries resolve, and `gfortran --version` works. `benchmarks/perf/setup_external.py` configures and compiles the Python interface inside the disposable job environment and imports FEDB if provided. Avoid fabricated version overrides; record actual installed versions. Geant4 needs no external machine installation or variables: the job installs `bdsim-g4`, its pinned HDF5/expat dependencies, and the C++ toolchain in its own Conda environment, then compiles the interface. See the current workflow for exact package pins.

## 7. History branch and rulesets

The first successful history publication creates the orphan `benchmark-results` branch if absent. Existing history must not be reinitialised. It contains immutable `runs/` and `pr-runs/` records and rolling `baselines/` snapshots. Set an active branch ruleset in **Settings → Rules → Rulesets** named `Benchmark history integrity` targeting exactly `benchmark-results` (API: `refs/heads/benchmark-results`, no quotes). Enable restrict deletions, block force pushes, and require linear history. Do not require PRs or required checks for this branch: the publisher pushes directly. Do not enable restrict updates/creations unless an appropriate narrowly scoped publisher identity is explicitly allowed to bypass; the current job token must be able to create and append. The existing ruleset (ID 27218022, checked during review) has a quoted ref pattern `refs/heads/"benchmark-results"`; correct it and verify the rules are effective on the actual branch.

Protect main separately with reviewed PRs and appropriate required checks; use CODEOWNERS/review policy for workflow and publisher changes, and avoid broad bypasses. The history rules prevent deletion, non-fast-forward pushes, and merges; they do not guarantee append-only file contents against arbitrary authorised fast-forward commits. The current publishers refuse overwriting existing run records, but stronger writer isolation requires a dedicated GitHub App/identity and workflow changes. Do not introduce a broad PAT just to bypass the ruleset.

## 8. Pages and first activation

Set repository **Settings → Pages → Build and deployment → Source: GitHub Actions**. Check the `github-pages` environment allows deployments from main (the trusted workflow ref), and does not require manual approval for scheduled jobs unless intended. Pages is built from `benchmark-results` by the trusted hosted jobs; never deploy from PR code. If the history branch does not yet exist, dispatch the reference workflow on main, target `main`, profile `quick`, with `publish_history=true` to bootstrap it. Inspect the resulting record and website, then run a full reference publication to test the full path. A manual full run does not establish/replace the scheduled rolling baseline; the first successful scheduled full run does that. The dashboard-only workflow requires an existing history branch.

Smoke test in order:

1. Main quick reference, publication disabled: runner online, job completes, measurement artifact exists, no permanent run is written.
2. Main quick with publication enabled: history branch/record and Pages deployment succeed.
3. Main full, then FLUKA and Geant4 manual runs with publication disabled: real measurements and dependency metadata, no external-engine skip masquerading as success.
4. Manual feature ref (select workflow ref main, set `target_ref` to branch/tag/SHA): uses `ubuntu-24.04`, not CERN.
5. Reviewed test PR: paired base/candidate benchmarks, read-only measurement token, successful trusted publisher, PR record appears on dashboard.
6. Manual dashboard redeployment; verify all historic runs and PR comparisons are still present. Confirm scheduled jobs create compact records and advance only their scheduled profile baseline. Manual publishing defaults to off.

The published site is https://xsuite.github.io/xcoll/ (verify the configured Pages URL in repository settings). Raw PR artefacts are retained 90 days; permanent compact history is stored on the history branch. Check scheduled workflows are enabled and main is the default branch. Do not activate the template as an additional schedule.

## 9. Local development and maintenance

```bash
# In a working Xcoll environment, from repository root:
python -m pip install -e '.[benchmarks]'
cd benchmarks
python -m pytest -m quick
python -m perf.results list
python -m perf.report
# Full CPU:
python -m pytest -m 'not fluka and not geant4 and not gpu'
```

No registration token is required for local measurements. Use a separate worktree or keep the checkout unchanged during runs; run serially, not with xdist. For reproducible comparisons use pinned dependencies and the same machine. On the VM monitor `free -h`, `df -h`, runner logs, and `journalctl -k` for OOM kills. The CPU job has a 120-minute timeout, external engines 180 minutes, and the hosted PR job 150 minutes. Avoid running competing tasks or changing the flavour mid-series; increment the testbed/physical-host identifiers if a new machine or configuration should start a distinct series. Keep system updates and runner maintenance scheduled outside benchmarks; restart the service after maintenance and confirm connectivity. If replacing the VM, remove its old runner registration, provision a new one, and explicitly decide whether it is a new measurement series.
