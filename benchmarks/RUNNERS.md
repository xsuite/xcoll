# CERN OpenStack benchmark runner

The reference workflow uses a dedicated GitHub self-hosted runner with the label `xcoll-bench-cpu`. Its Linux service must run under the unprivileged `xcollrunner` account. The runner communicates with GitHub over outbound HTTPS.

## Install Miniforge once

As `xcollrunner` on an x86-64 Linux machine:

```bash
cd ~
curl -fsSL -o miniforge.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash miniforge.sh -b -p "$HOME/miniforge3"
rm miniforge.sh
~/miniforge3/bin/conda --version
```

The workflow expects `/home/xcollrunner/miniforge3/bin/conda`. Adapt `CONDA_EXE` in the workflow if the installation is elsewhere. `conda init` is unnecessary. A compiler, Git, and the system libraries required by Xsuite must also be available.

Every benchmark run creates a fresh `python=3.12` environment in `$RUNNER_TEMP`, installs the latest released Xsuite dependencies and the checked-out Xcoll, regenerates serial kernels, and saves dependency metadata. The environment is removed after the results have been uploaded. The runner itself does not require a system Python.

## GitHub runner service

Register the VM under **Settings → Actions → Runners** in the Xcoll repository, following GitHub's Linux self-hosted runner instructions. Supply the label `xcoll-bench-cpu` and use a dedicated unprivileged account. Start the registered runner as a systemd service:

```bash
cd /home/xcollrunner/actions-runner
sudo ./svc.sh install xcollrunner
sudo ./svc.sh start
sudo ./svc.sh status
```

On AlmaLinux with SELinux enabled, systemd may reject execution of `runsvc.sh` with `status=203/EXEC` if its SELinux context is inappropriate. Diagnose the denial with `ausearch`; apply a narrow, persistent file-context rule rather than disabling SELinux.

## Security

**Do not run public PR code on the persistent CERN VM.** The reference workflow explicitly checks out reviewed `main` on the CERN runner. Manual runs of other revisions and the PR workflow use GitHub-hosted machines, and only the separate publisher receives a write token.

This is necessary but not a complete security boundary. Public-repository PRs can modify workflow definitions, so an organisation administrator should restrict which workflows can target the persistent runner. For stronger isolation, place the runner behind a private orchestration repository and have it fetch only reviewed revisions. Avoid CERN authentication credentials, SSH keys, private mounts, and privileged network access in the runner account.

GitHub's central `xsuite/xsuite` tests already use self-hosted runner labels (for example, `alma-cpu-1`). The Xcoll benchmark runner uses its own label and should not consume the existing runner pool.

## Measurement consistency

Benchmark runs share one VM, but its OpenStack hypervisor and host allocation can still vary. Use a fixed VM flavour, control competing workloads, record the CPU model and OS, and avoid comparing results from different testbeds. Keep a separate performance series when the CPU flavour, compiler, environment, or workload protocol changes.

The workflow installs the latest *released* Xsuite packages at each run. This is an integration history, not an environment-locked experiment. If a performance change must be attributed to Xcoll alone, repeat with pinned dependency versions.

## Activation

The active reference workflow is `.github/workflows/xcoll-performance.yml`. The identical CERN template is also included under `benchmarks/templates/xcoll-performance-cern.yml` for reuse. **Do not enable both as separate scheduled workflows.** The PR workflow remains `.github/workflows/xcoll-performance-pr.yml`.

To publish the dashboard, set the GitHub Pages source to **GitHub Actions**. The workflow commits immutable JSON to `benchmark-results`; the Pages site is built from all records on that branch. The publishing job may need organisation permission for repository write access.
