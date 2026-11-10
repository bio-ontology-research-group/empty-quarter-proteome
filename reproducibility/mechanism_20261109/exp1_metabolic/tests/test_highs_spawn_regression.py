#!/usr/bin/env python3
"""Regression test proving parent HiGHS solve followed by child HiGHS solve
under multiprocessing 'spawn' context does NOT hang/spin in highs::parallel::sync.
"""
from __future__ import annotations

import multiprocessing as mp
import time
import numpy as np
from scipy.optimize import milp, LinearConstraint


def run_milp_solve(threads: int = 1, time_limit: float = 5.0) -> tuple[bool, float, list[float]]:
    """Solve a small integer linear program with HiGHS."""
    c = np.array([-1.0, -1.0])
    A = np.array([[1.0, 2.0], [2.0, 1.0]])
    b_u = np.array([4.0, 4.0])
    b_l = np.array([0.0, 0.0])
    constraints = LinearConstraint(A, b_l, b_u)
    res = milp(
        c=c,
        constraints=constraints,
        integrality=np.ones(2),
        options={"threads": threads, "time_limit": time_limit},
    )
    if not res.success or res.fun is None:
        # Older local SciPy versions reject the 'threads' option in milp() and require OMP_NUM_THREADS=1
        res = milp(
            c=c,
            constraints=constraints,
            integrality=np.ones(2),
            options={"time_limit": time_limit},
        )
    return bool(res.success), float(res.fun), [float(x) for x in res.x]


def child_milp_worker(q: mp.Queue) -> None:
    """Child process entrypoint for spawn regression."""
    try:
        ok, val, sol = run_milp_solve(threads=1, time_limit=5.0)
        q.put({"success": ok, "fun": val, "sol": sol, "error": None})
    except Exception as e:
        q.put({"success": False, "fun": None, "sol": None, "error": str(e)})


def main() -> None:
    print("[TEST] Step 1: Solving real MILP in parent process...")
    t0 = time.time()
    p_ok, p_val, p_sol = run_milp_solve(threads=1, time_limit=5.0)
    dt_p = time.time() - t0
    print(f"       -> Parent solve: success={p_ok}, val={p_val}, sol={p_sol} in {dt_p:.3f}s")
    assert p_ok, "Parent solve failed"

    print("[TEST] Step 2: Spawning child process using 'spawn' context...")
    ctx = mp.get_context("spawn")
    q = ctx.Queue()
    child = ctx.Process(target=child_milp_worker, args=(q,))
    t_start = time.time()
    child.start()
    child.join(timeout=10.0)

    if child.is_alive():
        child.kill()
        child.join()
        raise TimeoutError("FATAL REGRESSION: Child process hung in highs::parallel::sync!")

    dt_c = time.time() - t_start
    assert not q.empty(), "Child process exited without putting result into queue"
    res = q.get()
    print(f"       -> Child solve: success={res['success']}, val={res['fun']}, sol={res['sol']} in {dt_c:.3f}s")
    assert res["success"], f"Child solve failed: {res['error']}"

    print(f"[TEST PASS] Parent HiGHS solve then child spawn HiGHS solve succeeded cleanly in {dt_c:.3f}s!")


def test_highs_spawn_regression() -> None:
    main()


if __name__ == "__main__":
    main()
