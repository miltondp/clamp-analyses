#!/usr/bin/env python

import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True)
    ap.add_argument("--mapping", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    s = pd.read_csv(a.summary, sep="\t", dtype={"phenotype": str})
    m = pd.read_csv(a.mapping).set_index("phenoplier_lv").metaLV
    lv = "lv" if "lv" in s.columns else "metaLV"
    assert set(s[lv]) == set(m.index), "summary LVs do not match the mapping"
    s.insert(0, "metaLV", s.pop(lv).map(m))
    s.to_csv(a.out, sep="\t", index=False)
    print(f"{s.metaLV.nunique()} metaLVs x {s.phenotype.nunique()} traits -> {a.out}")


if __name__ == "__main__":
    main()
