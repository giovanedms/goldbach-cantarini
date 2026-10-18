"""Extract the numerical tables and per-N data of the paper from the certificate outputs.

Writes into this directory:
  table2_zeros.csv          T(N) versus the truncated double series T_K(N) (Table 2)
  table3_fits.csv           fits of F(N)/N in implementations A and B (Table 3)
  data_F_implementation_A.csv / data_F_implementation_B.csv   per-N values used in the fits
Run from the repository root:  python results/extract_tables.py
"""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CERT = ROOT / "certificate"
OUT = ROOT / "results"
NUM = r"[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?"


def write(name, header, rows):
    with open(OUT / name, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f"{name}: {len(rows)} rows")


def table2():
    rows, cur = [], None
    for line in (CERT / "implementation_B" / "zeros_check_N25_58_74.out").read_text().splitlines():
        m = re.match(rf"N=(\d+) T_stieltjes=({NUM})\s+P\(N\)\^2=({NUM})\s+T_stieltjes-2P\^2=({NUM}) "
                     rf"T_stieltjes\+2P\^2=({NUM})", line)
        if m:
            if cur:
                rows.append(cur)
            cur = list(m.groups())
            continue
        m = re.match(rf"\s*K=(\d+) gamma_max=\d+ T_K=({NUM})", line)
        if m and cur is not None:
            cur.append(m.group(2))
    if cur:
        rows.append(cur)
    write("table2_zeros.csv", ["N", "T", "P2", "T_minus_2P2", "T_plus_2P2",
                               "T_K1000", "T_K4000", "T_K16000", "T_K40000"], rows)


def table3():
    a = (CERT / "implementation_A" / "fit_F_5e6.out").read_text()
    b = (CERT / "implementation_B" / "fit_F_2e5_6e6.out").read_text()
    g = lambda pat, t: re.search(pat, t).group(1)
    rows = [
        ["n", g(r"n=(\d+)", a), g(r"n=(\d+)", b)],
        ["c", g(rf"fit F/N = ({NUM})", a), g(rf"fit F/N = ({NUM})", b)],
        ["beta", g(rf"\+ ({NUM}) \* P\^2/N", a), g(rf"\+ ({NUM})\*P\^2/N", b)],
        ["sd_residual", g(rf"resid sd=({NUM})", a), g(rf"resid sd=({NUM})", b)],
        ["mean_(F-P2)/N", g(rf"\(F-P\^2\)/N: mean=({NUM})", a), g(rf"\(F-P\^2\)/N: mean=({NUM})", b)],
        ["sd_(F-P2)/N", g(rf"\(F-P\^2\)/N: mean={NUM} sd=({NUM})", a), g(rf"\(F-P\^2\)/N: mean={NUM} sd=({NUM})", b)],
        ["max_P2/N", g(rf"P\^2/N max=({NUM})", a), g(rf"P\^2/N max=({NUM})", b)],
        ["minus_4_log_2pi", g(rf"-4 log\(2 ?pi\) = ({NUM})", a), g(rf"-4 log\(2 ?pi\) = ({NUM})", b)],
    ]
    write("table3_fits.csv", ["quantity", "implementation_A", "implementation_B"], rows)


def data():
    rows = []
    for line in (CERT / "implementation_A" / "exact_F_5e6.out").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        rows.append(line.split())
    write("data_F_implementation_A.csv",
          ["N", "F_over_N", "P2_over_N", "F_minus_P2_over_N", "T_over_N", "rG2N_over_2N", "special"],
          [r + [""] * (7 - len(r)) for r in rows])
    with open(CERT / "implementation_B" / "exact_F_2e5_6e6.tsv") as f:
        rd = csv.reader(f, delimiter="\t")
        header = next(rd)
        write("data_F_implementation_B.csv", header, list(rd))


if __name__ == "__main__":
    table2(); table3(); data()
