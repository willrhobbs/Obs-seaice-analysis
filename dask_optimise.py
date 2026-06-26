# Detect Available Resources from an HPC instance

import os
import multiprocessing
import psutil
import dask

def get_available_resources():
    """Detect CPU, memory, and scratch space available to this job."""
    
    # --- CPUs ---
    # PBS/SLURM set these; fall back to system count
    ncpus = int(
        os.environ.get("PBS_NCPUS") or
        os.environ.get("SLURM_CPUS_PER_TASK") or
        os.environ.get("SLURM_NPROCS") or
        multiprocessing.cpu_count()
    )

    # --- Memory ---
    # Try cgroup memory limit first (most reliable on HPC)
    mem_bytes = None
    cgroup_path = "/sys/fs/cgroup/memory/memory.limit_in_bytes"
    if os.path.exists(cgroup_path):
        with open(cgroup_path) as f:
            val = int(f.read().strip())
            # Ignore sentinel "unlimited" value
            if val < 2**62:
                mem_bytes = val

    if mem_bytes is None:
        # Fall back to PBS/SLURM env vars (convert "16GB" style strings)
        raw = os.environ.get("PBS_VMEM") or os.environ.get("SLURM_MEM_PER_NODE")
        if raw:
            mem_bytes = _parse_mem_string(raw)

    if mem_bytes is None:
        # Last resort: total system RAM
        mem_bytes = psutil.virtual_memory().total

    # --- Scratch / jobFS ---
    scratch = (
        os.environ.get("PBS_JOBFS") or
        os.environ.get("TMPDIR") or
        "/tmp"
    )
    scratch_bytes = psutil.disk_usage(scratch).free

    return ncpus, mem_bytes, scratch, scratch_bytes


def _parse_mem_string(s):
    """Parse strings like '16GB', '32000MB', '1TB' to bytes."""
    s = s.strip().upper()
    units = {"KB": 2**10, "MB": 2**20, "GB": 2**30, "TB": 2**40}
    for unit, factor in units.items():
        if s.endswith(unit):
            return int(float(s[:-len(unit)]) * factor)
    return int(s)  # assume bytes if no unit


##################################################
#compute optimal Dask settings from availavble resources
def get_dask_config(memory_fraction=0.75, threads_per_worker=2, allow_spill=False):
    ncpus, mem_bytes, scratch, _ = get_available_resources()

    n_workers = max(1, ncpus // threads_per_worker)
    memory_per_worker = int(mem_bytes * memory_fraction) // n_workers

    dask.config.set({
        "distributed.worker.memory.target": 0.60,
        "distributed.worker.memory.spill": 0.95 if not allow_spill else 0.70,
        "distributed.worker.memory.pause": 0.80,
        "distributed.worker.memory.terminate": 0.90,
    })

    config = {
        "n_workers": n_workers,
        "threads_per_worker": threads_per_worker,
        "memory_limit": memory_per_worker,
    }

    if allow_spill:
        config["local_directory"] = scratch

    return config

# 