#!/usr/bin/env python
# Label the top ARCHS4 samples of each MGM and test which sample types are enriched.
import argparse

import numpy as np
import pandas as pd

CATEGORIES = [
    ("Monocyte/macrophage", r"macrophag|monocyt|microglia|thp-?1|u937|kupffer|bronchoalveolar|\bbalf?\b"),
    ("Endoderm/pancreas", r"endoderm|pancrea|islet|beta[ -]?cell|endoc-?beta|ngn3|pdx1|diabet"),
    ("Heart", r"heart|cardi|myocard|atrial|atrium|ventric"),
    ("Liver", r"liver|hepat|hepg2|huh-?7|hep3b|hepa"),
    ("Intestine", r"intestin|colon|colorect|rectum|duoden|ileum|jejun|enteroid|caco-?2|ht-?29|hct-?116"),
    ("Brain", r"brain|neur|astro|oligoden|glia|cortex|cerebell|hippocamp"),
    ("Adipose", r"adipo|\bfat\b|sgbs"),
    ("Lung", r"lung|bronch|alveol|a549|airway"),
    ("Blood/immune", r"pbmc|blood|lymph|t[ -]cell|b[ -]cell|\bcd4|\bcd8|neutrophil|leukocyte|\bnk\b"),
    ("Stem cell", r"embryonic stem|\besc\b|hesc|ipsc|pluripoten|\bhues|\bh1\b|\bh9\b"),
    ("Tumour/cancer", r"tumor|tumour|cancer|carcinoma"),
]


def label_samples(b):
    text = (b.source.fillna("") + " | " + b.title.fillna("") + " | " + b.characteristics.fillna("")).str.lower()
    label = pd.Series(pd.NA, index=b.index, dtype="object")
    for name, pattern in CATEGORIES[::-1]:
        label[text.str.contains(pattern, regex=True)] = name
    return label.fillna("Other")


def study_shares(b, mask):
    return b.loc[mask, ["study", "label"]].groupby("study").label.value_counts(normalize=True).unstack(fill_value=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--b", required=True)
    ap.add_argument("--selected", required=True)
    ap.add_argument("--top-fraction", type=float, default=0.01)
    ap.add_argument("--min-pct", type=float, default=5)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261114)
    ap.add_argument("--enrichment-out", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    mgms = list(pd.read_csv(a.selected).metaLV)
    b = pd.read_csv(a.b, sep="\t", usecols=["sample", "series_id", "source", "title", "characteristics"] + mgms)
    b["study"] = b.series_id.astype(str).str.split(",").str[0]
    b["label"] = label_samples(b)
    rng = np.random.default_rng(a.seed)
    rows = []
    for m in mgms:
        top = b.index.isin(b.nlargest(round(len(b) * a.top_fraction), m).index)
        t_tab, r_tab = study_shares(b, top), study_shares(b, ~top)
        it = rng.integers(0, len(t_tab), (a.n_boot, len(t_tab)))
        ir = rng.integers(0, len(r_tab), (a.n_boot, len(r_tab)))
        for c in t_tab.columns:
            if c == "Other":
                continue
            t = t_tab[c].values
            r = r_tab[c].values if c in r_tab else np.zeros(len(r_tab))
            boot = t[it].mean(1) / np.maximum(r[ir].mean(1), 1e-12)
            lo, hi = np.percentile(boot, [2.5, 97.5])
            rows.append(dict(metaLV=m, category=c, top_pct=100 * t.mean(), rest_pct=100 * r.mean(), FC=t.mean() / r.mean(),
                             ci_low=lo, ci_high=hi, n_studies=int((t_tab[c] > 0).sum())))
    res = pd.DataFrame(rows).sort_values(["metaLV", "top_pct"], ascending=[True, False])
    res["shown"] = res.top_pct.ge(a.min_pct) & res.ci_low.gt(1)
    res.to_csv(a.enrichment_out, index=False)

    label = f"top {100 * a.top_fraction:g}%"
    out = []
    for m in mgms:
        keep = res[(res.metaLV == m) & res.shown].sort_values("top_pct", ascending=False)
        out += [dict(metaLV=m, set=label, category=r.category, pct=r.top_pct, FC=r.FC) for r in keep.itertuples()]
        out.append(dict(metaLV=m, set=label, category="other", pct=100 - keep.top_pct.sum(), FC=np.nan))
    out = pd.DataFrame(out)
    out.to_csv(a.out, index=False)
    print("labelled share of ARCHS4:", round((b.label != "Other").mean(), 3))
    print(out.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
