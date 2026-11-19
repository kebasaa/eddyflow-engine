"""Build the many-gas fixture: 80 gas records, past the old ceiling of 64.

MaxNumGases went from 64 to 128 so that a PTR-MS dataset - innFLUX's use case
(Striednig et al. 2020), with a tracer per mass - fits in one run. Every other
fixture has at most eight gases, so none of them can show that slots 65 and up
are read, processed and written. This one does.

From the same three hours of CH-LAE the other fixtures use, it writes a copy of
the six raw files with 72 columns appended: the MIRO's COS (column 10) delayed
by 1..72 rows, so every copy has a different, known extra lag. The project is
base_n_gas's eight gases plus 72 COS records on those columns - 80 in all, the
same species repeated, which the engine allows (base_dup).

Gate: the run completes, every FLUXNET row matches its header, the header has
a column group for COS_80, and each copy's lag exceeds the original's by its
delay (within the search window).

Usage:  python gen_many.py [--data-out DIR]

Writes base_many.eddyflow and base_many.metadata beside this script and the
raw files to data_many/ beside it, which is gitignored. The source is the
Lagern repository's own data folder, read here and never written into.
"""

import argparse
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE_PROJECT = "base_n_gas.eddyflow"
SOURCE_METADATA = "base_site.metadata"
HEADER_ROWS = 1
COS_COL = 10
N_EXTRA = 72
FILL = "-9999"
FILES = [
    f"CH-LAE_ec_preproc_10hz_20250601-{h:02d}{m:02d}.csv"
    for h in range(3)
    for m in (0, 30)
]


def ini_value(text, key):
    m = re.search(rf"^{re.escape(key)}=(.*)$", text, re.M)
    if not m:
        raise SystemExit(f"{key} not found")
    return m.group(1).strip()


def set_key(text, key, value):
    pat = re.compile(rf"^{re.escape(key)}=.*$", re.M)
    if not pat.search(text):
        raise SystemExit(f"{key} not found")
    return pat.sub(f"{key}={value}", text)


def extend(src, dst):
    rows = src.read_text(encoding="utf-8", errors="replace").splitlines()
    head, data = rows[:HEADER_ROWS], rows[HEADER_ROWS:]
    fields = [r.split(",") for r in data]
    cos = [f[COS_COL - 1] if len(f) >= COS_COL else FILL for f in fields]
    out = []
    for h in head:
        out.append(h + "".join(f",COS_D{k}" for k in range(1, N_EXTRA + 1)))
    for i, f in enumerate(fields):
        extra = [cos[i - k] if i - k >= 0 else FILL for k in range(1, N_EXTRA + 1)]
        out.append(",".join(f + extra))
    dst.write_text("\n".join(out) + "\n", encoding="utf-8", newline="")
    return len(fields[0]) if fields else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-out", default=None)
    args = ap.parse_args()

    project = (HERE / SOURCE_PROJECT).read_text(encoding="utf-8")
    metadata = (HERE / SOURCE_METADATA).read_text(encoding="utf-8")
    src_dir = Path(ini_value(project, "data_path"))
    out_dir = Path(args.data_out) if args.data_out else HERE / "data_many"
    out_dir.mkdir(parents=True, exist_ok=True)

    ncol = None
    for name in FILES:
        src = src_dir / name
        if not src.exists():
            raise SystemExit(f"missing source file {src}")
        ncol = extend(src, out_dir / name)
    print(f"{len(FILES)} files, {ncol} source columns + {N_EXTRA}")

    # Metadata: the source's columns, then the 72 copies, described like COS.
    described = max(int(n) for n in re.findall(r"^col_(\d+)_variable=", metadata, re.M))
    if described > ncol:
        raise SystemExit(f"metadata describes {described} columns, files have {ncol}")
    block = re.search(rf"^col_{COS_COL}_variable=.*?(?=^col_{COS_COL + 1}_variable=)",
                      metadata, re.M | re.S).group(0)
    #> The rows carry more fields than the metadata names (a trailing one);
    #> the engine reads by position, so those are declared ignore and the
    #> copies go after them.
    last_block = re.search(rf"^col_{described}_variable=.*?(?=^\[|\Z)",
                           metadata, re.M | re.S).group(0)
    pad = []
    for c in range(described + 1, ncol + 1):
        blk = re.sub(rf"^col_{described}_", f"col_{c}_", last_block, flags=re.M)
        pad.append(re.sub(rf"^col_{c}_variable=.*$", f"col_{c}_variable=ignore", blk, flags=re.M))
    last = ncol
    extra = list(pad)
    for k in range(1, N_EXTRA + 1):
        extra.append(re.sub(rf"^col_{COS_COL}_", f"col_{last + k}_", block, flags=re.M))
    anchor = re.search(rf"^col_{described}_variable=.*?(?=^\[|\Z)", metadata, re.M | re.S)
    meta = metadata[:anchor.end()].rstrip("\n") + "\n" + "".join(extra) + metadata[anchor.end():]
    (HERE / "base_many.metadata").write_text(meta, encoding="utf-8")

    # Project: base_n_gas's records, then one COS record per copy.
    n0 = int(ini_value(project, "gas_num"))
    cos_rec = next(i for i in range(1, n0 + 1)
                   if ini_value(project, f"gas_{i}_var") == "cos")
    recs = []
    for k in range(1, N_EXTRA + 1):
        g = n0 + k
        recs += [f"gas_{g}_var=cos",
                 f"gas_{g}_instr={ini_value(project, f'gas_{cos_rec}_instr')}",
                 f"gas_{g}_col={last + k}",
                 f"gas_{g}_moist=0", f"gas_{g}_cell=0",
                 f"gas_{g}_mw={ini_value(project, f'gas_{cos_rec}_mw')}",
                 f"gas_{g}_diff={ini_value(project, f'gas_{cos_rec}_diff')}"]
    p = set_key(project, "gas_num", str(n0 + N_EXTRA))
    anchor = re.search(rf"^gas_{n0}_diff=.*$", p, re.M)
    p = p[:anchor.end()] + "\n" + "\n".join(recs) + p[anchor.end():]
    p = set_key(p, "data_path", out_dir.as_posix())
    p = set_key(p, "proj_file", (HERE / "base_many.metadata").as_posix())
    (HERE / "base_many.eddyflow").write_text(p, encoding="utf-8")
    print(f"wrote base_many.eddyflow with {n0 + N_EXTRA} gas records")


if __name__ == "__main__":
    main()
