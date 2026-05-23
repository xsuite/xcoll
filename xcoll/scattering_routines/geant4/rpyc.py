# copyright ############################### #
# This file is part of the Xcoll Package.   #
# Copyright (c) CERN, 2025                  #
# ######################################### #

import re
import os
import time
import threading
import subprocess
from shutil import which


PORT_PATTERNS = [
    re.compile(r"server started on \[[^\]]+\]:(\d+)"),        # classic banner
    re.compile(r"^\s*\d{1,3}(?:\.\d{1,3}){3}\s+(\d+)\s*$"),   # "127.0.0.1<tab>PORT"
]
def extract_port(line):
    for pat in PORT_PATTERNS:
        m = pat.search(line)
        if m:
            return int(m.group(1))


# Thread-safe launch of rpyc server on a free port
def launch_rpyc_with_port(log_path="rpyc.log", timeout=10.0):
    exe = which("rpyc_classic")
    if not exe:
        raise RuntimeError("rpyc_classic not on PATH; install rpyc in this env")

    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"

    p = subprocess.Popen(
        [exe, '-m', 'oneshot', '-p', '0'],  # -u = unbuffered
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,                # line-buffered wrapper on our side
        close_fds=True,
        env=env,
    )

    logfile = open(log_path, "a", buffering=1)
    port = None
    deadline = time.time() + timeout

    # Read until we see the "server started ..." line (or timeout/exit)
    while True:
        line = p.stdout.readline()
        if line == "":
            # EOF or just no data yet; only bail if the proc exited
            if p.poll() is not None:
                break
            if time.time() > deadline:
                break
            time.sleep(0.02)
            continue

        logfile.write(line)
        port = extract_port(line)
        if port:
            break

        if time.time() > deadline:
            break

    if port is None:
        # drain remaining output before raising (helps debugging)
        for rest in p.stdout:
            logfile.write(rest)
        logfile.close()
        raise RuntimeError("Failed to detect rpyc port (timed out or process exited)")

    # Keep teeing the rest of the output to the same log
    def _pump():
        for line in p.stdout:
            logfile.write(line)
        logfile.close()

    threading.Thread(target=_pump, daemon=True).start()
    return p, port
