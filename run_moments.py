"""
Production runner: exact moments of the Toeplitz/Hankel LSDs via orbit-reduced exact
polytope volumes, parallelized, checkpointed, with per-orbit Monte Carlo cross-checks.

Usage:  python run_moments.py T 10     (kind 2k)
        python run_moments.py H 12 --mc 2000000
Results: results/volumes_{kind}{2k}.json  (orbit rep -> {size, volume "p/q", mc})
"""
from __future__ import annotations
import json
import os
import sys
import time
from fractions import Fraction
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tk_core  # noqa: E402

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


def worker(args):
    pm, kind, mc_samples = args
    import exact_volume
    import tk_core as tc
    t0 = time.time()
    vol = exact_volume.pairing_volume_exact(pm, kind)
    dt = time.time() - t0
    mc = None
    if mc_samples:
        rows = tc.dedupe_rows((tc.toeplitz_forms if kind == "T" else tc.hankel_forms)(pm))
        mc = tk_core.mc_volume(rows, kind, n_samples=mc_samples, seed=hash(pm) % (2**31))
    return pm, str(vol), mc, dt


def main():
    kind = sys.argv[1]
    n2k = int(sys.argv[2])
    mc_samples = 0
    if "--mc" in sys.argv:
        mc_samples = int(sys.argv[sys.argv.index("--mc") + 1])
    workers = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else max(1, os.cpu_count() - 2)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    outfile = os.path.join(RESULTS_DIR, f"volumes_{kind}{n2k}.json")

    pair_iter = tk_core.all_pairings(n2k) if kind == "T" else tk_core.parity_pairings(n2k)
    orbits = tk_core.orbit_reps(pair_iter)
    reps = sorted(orbits.keys())
    print(f"[{kind}{n2k}] {sum(orbits.values())} pairings in {len(reps)} dihedral orbits; "
          f"workers={workers} mc={mc_samples}", flush=True)

    done = {}
    if os.path.exists(outfile):
        with open(outfile) as f:
            saved = json.load(f)
        for key, rec in saved.get("orbits", {}).items():
            done[tuple(json.loads(key))] = rec
        print(f"resuming: {len(done)} orbits already computed", flush=True)

    todo = [r for r in reps if r not in done]

    def flush():
        payload = {
            "kind": kind, "n2k": n2k,
            "n_pairings": sum(orbits.values()), "n_orbits": len(reps),
            "orbits": {json.dumps(list(r)): done[r] for r in reps if r in done},
        }
        if len(done) == len(reps):
            total = sum(Fraction(done[r]["volume"]) * orbits[r] for r in reps)
            payload["moment_exact"] = str(total)
            payload["moment_float"] = float(total)
        tmp = outfile + ".tmp"
        with open(tmp, "w") as f:
            json.dump(payload, f, indent=1)
        os.replace(tmp, outfile)

    t0 = time.time()
    n_done_now = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(worker, (r, kind, mc_samples)): r for r in todo}
        for fut in as_completed(futs):
            pm, vol, mc, dt = fut.result()
            rec = {"size": orbits[pm], "volume": vol, "time_s": round(dt, 2)}
            if mc is not None:
                rec["mc"] = round(mc, 5)
                rec["mc_abs_err"] = round(abs(mc - float(Fraction(vol))), 5)
            done[pm] = rec
            n_done_now += 1
            if n_done_now % 5 == 0 or n_done_now == len(todo):
                flush()
                el = time.time() - t0
                print(f"  {len(done)}/{len(reps)} orbits  ({el:.0f}s elapsed, last {dt:.1f}s)", flush=True)
    flush()

    total = sum(Fraction(done[r]["volume"]) * orbits[r] for r in reps)
    print(f"[{kind}{n2k}] EXACT m = {total}  = {float(total):.6f}", flush=True)
    if mc_samples:
        worst = max((done[r].get("mc_abs_err", 0) for r in reps), default=0)
        print(f"[{kind}{n2k}] worst per-orbit |MC - exact| = {worst}", flush=True)


if __name__ == "__main__":
    main()
