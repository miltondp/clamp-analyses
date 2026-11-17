import math
import sys
from pathlib import Path
from schematic_lib import *

HERE = Path(__file__).parent
IC = HERE / "icons"
W, H = 180, 184
c = Canvas(W, H, "CLAMP overview: prior-informed gene modules, consensus reference and reuse")
BLUE_BOX, TEAL_BOX, PURPLE_BOX = ("#f7fafc", BLUE_LINE), ("#f3faf9", TEAL), ("#f8f6fb", "#b7aad3")
RED = "#c0392b"
BSEQ = ["#f1f6fb", "#dde9f4", "#c7dbee", "#aecae4", "#94b8d9", "#7aa5ce"]
GSEQ = ["#f4f5f6", "#e2e5e8", "#cbd0d5", "#b0b7be", "#949ea7", "#7a858f"]
GLINE = "#7a858f"
TSEQ = ["#f2f9f8", "#d6ece9", "#b5dcd7", "#8fc9c2", "#5fb0a6", "#33978c"]
DIV = ["#2f6fa8", "#7aa8d2", "#c7dbee", "#f3c7bb", "#e08670", "#c0392b"]


def box(x, y, w, h, style, title, sub=None):
    c.rect(x, y, w, h, fill=style[0], stroke=style[1], sw=0.35, rx=1)
    c.text(x + w / 2, y + 4.4, title, size=TITLE, weight="bold")
    if sub:
        c.text(x + w / 2, y + 7.4, sub, size=SMALL)


def wbox(x, y, w, h, rows, size=SMALL):
    c.rect(x, y, w, h, fill="white", stroke=INK, sw=0.25, rx=0.6)
    lh = size * PT * 1.18
    y0 = y + h / 2 - (len(rows) - 1) * lh / 2 + size * PT * 0.35
    for i, r in enumerate(rows):
        c.text(x + w / 2, y0 + i * lh, r, size=size)


def bell(x, y, w, hgt, wide):
    n = 30
    sd = w / (2.3 if wide else 7.5)
    pts = [(x + w * k / n, y - hgt * math.exp(-((w * k / n) ** 2) / (2 * sd ** 2))) for k in range(n + 1)]
    p = " ".join(f"{a:.3f},{b:.3f}" for a, b in pts)
    if wide:
        c.add(f'<polygon points="{x:.3f},{y:.3f} {p} {x + w:.3f},{y:.3f}" fill="{TEAL_FILL}" fill-opacity="0.8" stroke="none"/>')
        c.add(f'<polyline points="{p}" fill="none" stroke="{TEAL}" stroke-width="0.3"/>')
    else:
        c.add(f'<polyline points="{p}" fill="none" stroke="#9ba7ae" stroke-width="0.3"/>')


def star(cx, cy, leaves, r=7.5, hub=3.4, leaf=1.5):
    pos = []
    for k, (fill, line) in enumerate(leaves):
        a = math.radians(-90 + k * 360 / len(leaves))
        x, y = cx + r * math.cos(a), cy + r * math.sin(a)
        c.line([(cx + hub * math.cos(a), cy + hub * math.sin(a)), (x, y)], sw=0.4, stroke=line)
        c.circle(x, y, leaf, fill, line, sw=0.3)
        pos.append((x, y))
    c.circle(cx, cy, hub, BLUE, BLUE_LINE, sw=0.35)
    return pos


def network(cx, cy, r, rings=(), diamond=None, fill="#e3f1ef", big=()):
    c.circle(cx, cy, r, fill)
    nodes = [(cx + r * 0.62 * math.cos(a), cy + r * 0.62 * math.sin(a)) for a in [k * 2 * math.pi / 7 + 0.4 for k in range(7)]] + [(cx, cy)]
    for i, j in [(0, 1), (1, 2), (2, 7), (7, 4), (3, 4), (4, 5), (5, 6), (6, 0), (7, 0), (2, 3)]:
        c.line([nodes[i], nodes[j]], sw=0.2, stroke="#7d8a94", dash="0.6,0.5")
    for k, (x, y) in enumerate(nodes):
        if k == diamond:
            d = 1.25
            c.add(f'<polygon points="{x:.3f},{y - d:.3f} {x + d:.3f},{y:.3f} {x:.3f},{y + d:.3f} {x - d:.3f},{y:.3f}" fill="{ORANGE_FILL}" stroke="{ORANGE}" stroke-width="0.3"/>')
            continue
        if k in rings:
            c.circle(x, y, 1.45, "none", YELLOW_LINE, sw=0.4)
        c.circle(x, y, 1.25 if k in big else 0.85, "#3aa197", TEAL, sw=0.2)
    return nodes


def zmat(x, y, w, h, ncols, nrows, seed=0, sym=None, sym_size=LABEL, solid=False):
    c.rect(x, y, w, h, fill="white", stroke="none")
    for r in range(nrows):
        for k in range(ncols):
            v = (r * 7 + k * 3 + seed * 5 + (r * k) % 4) % 6
            c.rect(x + w * k / ncols, y + h * r / nrows, w / ncols, h / nrows, fill=TSEQ[v], stroke="none")
    c.matrix(x, y, w, h, fill="none", stroke=TEAL, ncols=ncols, nrows=nrows)
    if sym:
        c.text(x + w / 2, y + h / 2 + sym_size * PT * 0.35, sym, size=sym_size, italic=True, halo="#e6f3f1")


def ymat(x, y, w, h, ncols, nrows, seed=0, sym=None, sym_size=TITLE, rowbias=None, **kw):
    import random
    rng = random.Random(seed)
    base = [rng.choice([0, 1, 2, 3, 4, 5]) for _ in range(nrows)]
    for r in range(nrows):
        for k in range(ncols):
            v = max(0, min(5, base[r] + rng.choice([-1, 0, 0, 1])))
            if rowbias:
                v = max(0, min(5, v + rowbias(r, k)))
            c.rect(x + w * k / ncols, y + h * r / nrows, w / ncols, h / nrows, fill=GSEQ[v], stroke="none")
    c.matrix(x, y, w, h, fill="none", stroke=GLINE, ncols=ncols, nrows=nrows, **kw)
    if sym:
        c.text(x + w / 2, y + h / 2 + sym_size * PT * 0.35, sym, size=sym_size, italic=True, halo="#eef0f2")


def steps(x, y, rows):
    for i, r in enumerate(rows):
        c.text(x, y + i * SMALL * PT * 1.18, r, size=SMALL, weight="bold" if i == 0 else "normal")


# ---------------- a ----------------
c.tag(1, 4, "a")
box(4, 1, 60.7, 41, BLUE_BOX, "CLAMPbase", "Semi-NMF of gene expression")
c.add('<g transform="translate(-1.7,-4.9) scale(1.12)">')
ymat(9.5, 18, 13, 15, 6, 7, seed=1, sym="Y", top="Samples", left="Genes")
c.text(16, 35.4, "Expression", size=SMALL)
for k in range(6):
    c.rect(13 + k * 1.0, 36.4, 1.0, 1.3, fill=GSEQ[k], stroke="none")
c.rect(13, 36.4, 6.0, 1.3, stroke=GLINE, sw=0.2)
c.text(12.4, 37.5, "Low", size=SMALL, anchor="end")
c.text(19.6, 37.5, "High", size=SMALL, anchor="start")
c.op(25.2, 25.5, "≈")
zpat = [[5, 0], [3, 1], [0, 4], [1, 2], [4, 0], [0, 3], [2, 5]]
for r, row in enumerate(zpat):
    for k, v in enumerate(row):
        c.rect(29.5 + k * 2.75, 18 + r * 15 / 7, 2.75, 15 / 7, fill=TSEQ[v], stroke="none")
c.matrix(29.5, 18, 5.5, 15, fill="none", stroke=TEAL, top="Modules", left="Genes", ncols=2, nrows=7)
c.text(32.25, 25.5 + TITLE * PT * 0.35, "Z", size=TITLE, italic=True, halo="#e6f3f1")
c.text(32.25, 35.4, "Gene loading", size=SMALL)
for k in range(6):
    c.rect(29.25 + k * 1.0, 36.4, 1.0, 1.3, fill=TSEQ[k], stroke="none")
c.rect(29.25, 36.4, 6.0, 1.3, stroke=TEAL, sw=0.2)
c.text(28.65, 37.5, "0", size=SMALL, anchor="end")
c.text(35.85, 37.5, "+", size=SMALL, anchor="start")
c.text(32.25, 40.8, "Z ≥ 0", size=TITLE, italic=True)
c.op(38, 25.5, "×")
bx, by, cw, ch = 42.5, 21, 1.15, 2.0
pat = [[4, 2, 5, 1, 3, 0, 4, 5, 2, 3], [1, 5, 3, 0, 4, 2, 5, 1, 0, 4], [5, 0, 4, 3, 1, 5, 2, 0, 3, 1], [2, 4, 1, 5, 0, 3, 1, 4, 5, 2]]
for r, row in enumerate(pat):
    for k, v in enumerate(row):
        c.rect(bx + k * cw, by + r * ch, cw, ch, fill=BSEQ[v], stroke="none")
c.rect(bx, by, cw * 10, ch * 4, stroke=BLUE_LINE, sw=0.3)
c.rect(bx, by + ch, cw * 10, ch, stroke="#4f7fae", sw=0.4)
c.text(bx + cw * 10 + 0.6, by + ch + 1.6, "j", size=SMALL, italic=True, anchor="start")
c.text(bx + 5.75, by - 1.0, "Samples", size=AXIS)
c.text(bx - 1.0, by + 4, "Modules", size=AXIS, rotate=-90)
c.text(bx + 5.75, by + 4 + TITLE * PT * 0.35, "B", size=TITLE, italic=True, halo="#e8f0f8")
c.text(bx + 5.75, 35.4, "Module activity", size=SMALL)
for k in range(6):
    c.rect(bx + 2.75 + k * 1.0, 36.4, 1.0, 1.3, fill=BSEQ[k], stroke="none")
c.rect(bx + 2.75, 36.4, 6.0, 1.3, stroke=BLUE_LINE, sw=0.2)
c.text(bx + 2.15, 37.5, "Low", size=SMALL, anchor="end")
c.text(bx + 9.35, 37.5, "High", size=SMALL, anchor="start")

c.add("</g>")
box(4, 44, 60.7, 26, PURPLE_BOX, "Cell-type signal in bulk", "Activity tracks known proportions")
c.add('<g transform="translate(4.6,0)">')
c.text(8.6, 55.4, "b_{j}", size=LABEL, italic=True)
for k, v in enumerate([1, 5, 3, 0, 4, 2, 5, 1, 0, 4]):
    c.rect(11 + k * 1.8, 53.6, 1.8, 2.0, fill=BSEQ[v], stroke="none")
c.rect(11, 53.6, 18, 2.0, stroke="#4f7fae", sw=0.4)
for k, f in enumerate([0.62, 0.35, 0.55, 0.28, 0.48, 0.2, 0.66, 0.58, 0.38, 0.45]):
    x0, top, hgt = 11.2 + k * 1.8, 57, 8.5
    c.rect(x0, top, 1.4, hgt * (1 - f), fill=GREY, stroke="none")
    c.rect(x0, top + hgt * (1 - f), 1.4, hgt * f, fill=PURPLE_FILL, stroke="none")
c.text(9.6, 61.3, "Cell types", size=SMALL, rotate=-90)
c.text(20, 68.2, "Known proportions", size=SMALL)
c.arrow((31, 61), (34.6, 61), sw=0.3)
c.line([(38, 54.5), (38, 66), (53, 66)], sw=0.25)
c.line([(39.2, 64.8), (52, 55.6)], sw=0.25, stroke=PURPLE, dash="0.8,0.6")
for x, y in [(39.8, 64.2), (41.6, 63.3), (43.6, 61.2), (44.9, 61.6), (46.7, 59.8), (48.7, 58.2), (50.4, 56.6), (51.6, 56.1)]:
    c.circle(x, y, 0.55, PURPLE_FILL, PURPLE, sw=0.2)
c.text(36.8, 60.3, "Module j", size=SMALL, rotate=-90)
c.text(45.5, 68.6, "Proportion of c", size=SMALL)


c.add("</g>")
box(66, 1, 112, 69, TEAL_BOX, "CLAMPfull", "Curated gene sets set how strongly each loading is shrunk")
c.add('<g transform="translate(-10,0)">')
gx, gy, gw, gh, nr = 92, 19, 12.5, 34, 10
rowh = gh / nr
c.text(gx + gw / 2, 13.6, "Curated priors", size=LABEL, weight="bold")
c.text(gx + gw / 2, 16.6, "Gene sets", size=AXIS)
c.rect(gx, gy, gw, gh, fill=GREY, stroke="none")
hits = {(0, 1), (1, 1), (1, 3), (2, 3), (5, 1), (8, 3)}
for (r, k) in hits:
    c.rect(gx + gw * k / 5, gy + rowh * r, gw / 5, rowh, fill=TEAL_FILL, stroke="none")
c.matrix(gx, gy, gw, gh, fill="none", stroke=GREY_LINE, ncols=5, nrows=nr, left="Genes")
for k in (1, 3):
    c.rect(gx + gw * k / 5, gy, gw / 5, gh, stroke=TEAL, sw=0.3)
c.text(gx + gw / 2, gy + gh + 3.0, "C", size=TITLE, italic=True)
cur = [8.5, 9.2, 5.5, 9.8, 2.6, 4.9, 2.1, 1.7, 3.6, 1.1]
for i, L in enumerate(cur):
    c.rect(107, gy + i * rowh + 0.6, L, rowh - 1.2, fill=TEAL, stroke="none")
c.line([(107, gy), (107, gy + gh)], sw=0.25)
c.text(107, 16.6, "z_{j} (current)", size=AXIS, italic=True, anchor="start")
c.arrow((118, 36), (120.6, 36), sw=0.3)
wbox(121, 26, 23, 20, [])
c.text(132.5, 29.6, "1  Select gene sets", size=AXIS, weight="bold")
c.lines(132.5, 33.0, ["AUROC prescreen,", "cross-validated non-", "negative elastic net"], size=SMALL)
c.text(124.4, 43.4, "u_{j}", size=AXIS, italic=True)
for k in range(5):
    c.rect(126 + k * 3.2, 41.2, 3.2, 2.8, fill=TEAL_FILL if k in (1, 3) else "white", stroke=GREY_LINE, sw=0.25)
c.arrow((145, 36), (147.6, 36), sw=0.3)
c.text(156.5, 13.6, "2  Prior on z_{gj}", size=AXIS, weight="bold")
c.text(156.5, 16.6, "centred at 0", size=SMALL)
wide_rows = {0, 1, 2, 5, 8}
for i in range(nr):
    bell(149.6, gy + (i + 1) * rowh - 0.35, 14.4, rowh * (0.34 if i in wide_rows else 0.88), i in wide_rows)
c.line([(149.6, gy), (149.6, gy + gh)], sw=0.25)
c.text(149.6, gy + gh + 2.2, "0", size=SMALL)
c.arrow((165, 36), (167.6, 36), sw=0.3)
c.text(172.6, 13.6, "3  Update", size=AXIS, weight="bold")
c.text(172.6, 16.6, "z_{j} (new)", size=AXIS, italic=True)
new = [3.9, 4.4, 2.7, 3.6, 0.5, 2.6, 0.0, 0.0, 2.2, 0.0]
for i, L in enumerate(new):
    if L:
        c.rect(169, gy + i * rowh + 0.6, L, rowh - 1.2, fill=PURPLE if i == 3 else TEAL, stroke="none")
c.line([(169, gy), (169, gy + gh)], sw=0.25)
c.text(173.2, gy + 3 * rowh + 1.9, "kept", size=SMALL, anchor="start", fill=PURPLE, weight="bold")
c.line([(171, 54), (171, 56.6), (132.5, 56.6)], sw=0.3)
c.arrow((132.5, 56.6), (132.5, 46.6), sw=0.3)
c.text(151.5, 59.4, "Re-estimated during training", size=SMALL)
lx, ly = 103, 64.4
c.add(f'<polygon points="{lx},{ly + 0.6} {lx},{ly - 0.4} {lx + 1.5},{ly - 0.3} {lx + 3},{ly + 0.1} {lx + 4.5},{ly + 0.45} {lx + 4.5},{ly + 0.6}" fill="{TEAL_FILL}" stroke="{TEAL}" stroke-width="0.25"/>')
c.text(lx + 5.6, ly + 0.6, "Gene in a selected set: wide prior, shrunk less", size=SMALL, anchor="start")
c.add(f'<polyline points="{lx + 47},{ly - 2.0} {lx + 47.4},{ly - 0.8} {lx + 48.0},{ly + 0.3} {lx + 51.5},{ly + 0.6}" fill="none" stroke="#9ba7ae" stroke-width="0.3"/>')
c.text(lx + 52.6, ly + 0.6, "Other genes: base penalty", size=SMALL, anchor="start")

c.add("</g>")

# ---------------- b ----------------
c.tag(1, 77, "b")
yb = 74
box(4, yb, 174, 36, TEAL_BOX, "ARCHS4 consensus modules", "Modules reproduced across independent CLAMPfull runs")
yc0 = yb + 18.5
ylab = yb + 30.6
ymat(9.5, yc0 - 6, 18, 12, 7, 5, seed=4, sym="Y", top="Samples", left="Genes")
steps(18.5, ylab, ["ARCHS4 expression"])
c.arrow((29.2, yc0), (38.4, yc0), sw=0.3)
for k in range(3):
    zmat(41.2 + k * 1.3, yc0 - 4.5 - k * 1.1, 10, 11, 4, 5, seed=k, sym="Z" if k == 2 else None)
c.text(48.8, yc0 - 7.9, "Modules", size=SMALL)
c.text(40.1, yc0 + 1, "Genes", size=SMALL, rotate=-90)
steps(47.5, ylab, ["Independent runs"])
c.arrow((56.2, yc0), (66.6, yc0), sw=0.3)
zmat(69, yc0 - 5, 15, 10, 6, 5, seed=3)
for xx in (72.75, 80.25):
    c.line([(xx - 0.9, yc0 - 5.9), (xx + 0.9, yc0 - 4.1)], sw=0.4, stroke=RED)
    c.line([(xx + 0.9, yc0 - 5.9), (xx - 0.9, yc0 - 4.1)], sw=0.4, stroke=RED)
steps(76.5, ylab, ["Noise screen"])
c.arrow((86, yc0), (95.4, yc0), sw=0.3)
dx, dy = 98.8, yc0 - 8
den = [(0, 13.5, 0, 10.5), (2.4, 13.5, 2.4, 10.5), (0, 10.5, 2.4, 10.5), (1.2, 10.5, 1.2, 8),
       (5, 13.5, 5, 11), (7.4, 13.5, 7.4, 11), (5, 11, 7.4, 11), (6.2, 11, 6.2, 8),
       (1.2, 8, 6.2, 8), (3.7, 8, 3.7, 4.5),
       (11, 13.5, 11, 9.5), (13.4, 13.5, 13.4, 9.5), (11, 9.5, 13.4, 9.5), (12.2, 9.5, 12.2, 4.5),
       (3.7, 4.5, 12.2, 4.5), (8, 4.5, 8, 2.5)]
for x1, y1, x2, y2 in den:
    c.line([(dx + x1, dy + y1), (dx + x2, dy + y2)], sw=0.3)
c.line([(dx - 1.5, dy + 6.6), (dx + 15, dy + 6.6)], sw=0.3, stroke=RED, dash="0.8,0.5")
for i, x in enumerate([0, 2.4, 5, 7.4, 11, 13.4]):
    c.circle(dx + x, dy + 14.4, 0.7, TEAL if i < 4 else PURPLE_FILL, TEAL if i < 4 else PURPLE, sw=0.2)
steps(105.5, ylab, ["Cluster across runs"])
c.arrow((116.2, yc0), (122.6, yc0), sw=0.3)
x0 = 124.85
for k in range(4):
    zmat(x0 + k * 2.6, yc0 - 6, 2.1, 12, 1, 6, seed=k + 1)
c.text(x0 + 5, yc0 - 7.2, "Members", size=SMALL)
c.arrow((x0 + 10.9, yc0), (x0 + 16.5, yc0), sw=0.3)
c.text(x0 + 13.7, yc0 - 1.0, "median", size=SMALL)
zmat(x0 + 17.7, yc0 - 6, 2.4, 12, 1, 6, seed=5)
c.text(x0 + 18.9, yc0 - 7.2, "Z_{MGM}", size=SMALL, italic=True)
steps(134.5, ylab, ["Consensus loadings"])
c.arrow((147.2, yc0), (154.4, yc0), sw=0.3)
mgm_nodes = network(163.5, yc0, 8.2, big=(0, 2, 7, 5))
steps(163.5, ylab, ["Meta-gene module (MGM)"])


# ---------------- c ----------------
c.tag(1, 117, "c")
yc, bw, bh = 114, 42.5, 33
xs = [4 + k * (bw + 1.33) for k in range(4)]
titles = [("Annotate modules", "Pathways, cell types, tissues"), ("Project new studies", "Module activities without retraining"),
          ("Test genetic associations", "PhenoPLIER GLS on S-MultiXcan"), ("Score drug-disease pairs", "Opposing trait and drug profiles")]
for k, (t, s) in enumerate(titles):
    box(xs[k], yc, bw, bh, ("#f6f7f8", "#aab3ba"), t, s)
x = xs[0]
c.add('<g transform="translate(0.6,3.6)">')
zmat(x + 6.5, yc + 11, 2.4, 11, 1, 6, seed=2)
c.text(x + 7.7, yc + 25, "Z_{MGM} column j", size=SMALL)
c.text(x + 12.3, yc + 14.2, "Rank", size=SMALL)
c.arrow((x + 9.6, yc + 16), (x + 15.0, yc + 16), sw=0.25)
for k, (wd, op) in enumerate([(10, 1), (8.5, .85), (7, .65), (5.5, .45)]):
    c.rect(x + 16, yc + 11.4 + k * 2, wd, 1.3, fill=TEAL, stroke="none", opacity=op)
c.text(x + 21, yc + 21.6, "Top 1% of genes", size=SMALL)
c.line([(x + 27.4, yc + 13.4), (x + 29, yc + 13.4)], sw=0.25)
c.line([(x + 29, yc + 11.4), (x + 29, yc + 16.2)], sw=0.25)
c.arrow((x + 29, yc + 11.4), (x + 30.8, yc + 11.4), sw=0.25)
c.arrow((x + 29, yc + 16.2), (x + 30.8, yc + 16.2), sw=0.25)
c.text(x + 31.4, yc + 12.2, "Pathways", size=SMALL, anchor="start")
c.text(x + 31.4, yc + 17.0, "Cell types", size=SMALL, anchor="start")
c.add("</g>")

x = xs[1]
c.rect(x + 4, yc + 10.6, 6.5, 0.9, fill=CONTROL, stroke="none")
c.rect(x + 10.5, yc + 10.6, 6.5, 0.9, fill=STUDY, stroke="none")
ymat(x + 4, yc + 12.4, 13, 10, 5, 4, seed=2, sym="Y′", sym_size=LABEL, left="Genes", rowbias=lambda r, k: (2 if k >= 3 and r in (0, 2) else 0))
c.text(x + 7.25, yc + 25, "Control", size=SMALL)
c.text(x + 13.75, yc + 25, "Treated", size=SMALL)
c.arrow((x + 18.2, yc + 17.4), (x + 21.4, yc + 17.4), sw=0.3)
c.line([(x + 24.5, yc + 10.6), (x + 24.5, yc + 24), (x + 39.5, yc + 24)], sw=0.25)
c.text(x + 23.4, yc + 17.4, "Module j", size=SMALL, rotate=-90)
for px, py in [(28, 21), (29.6, 22.2), (27.2, 22.6), (30.4, 20.4), (29, 19.8)]:
    c.rect(x + px - 0.6, yc + py - 0.6, 1.2, 1.2, fill=CONTROL, stroke="none")
for px, py in [(35.4, 12.6), (37, 13.4), (34.8, 14.4), (36.2, 15.4), (37.6, 12.0)]:
    c.rect(x + px - 0.6, yc + py - 0.6, 1.2, 1.2, fill=STUDY, stroke="none")
c.text(x + bw / 2, yc + 28.4, "B′ = (Z_{MGM}^{T}Z_{MGM} + λI)^{−1}Z_{MGM}^{T}Y′", size=LABEL, italic=True)

x = xs[2]
c.line([(x + 6, yc + 10.6), (x + 6, yc + 24), (x + 24, yc + 24)], sw=0.25)
c.text(x + 4.8, yc + 17.4, "−log10 P", size=SMALL, rotate=-90)
for px, py in [(9, 20.5), (10.5, 21.6), (11.6, 20.2), (9.8, 22.4), (12.4, 21.4)]:
    c.circle(x + px, yc + py, 0.6, GREY, GREY_LINE, sw=0.2)
c.line([(x + 8.4, yc + 21.3), (x + 13.2, yc + 21.3)], sw=0.35)
for px, py in [(18, 13.5), (19.6, 12.4), (20.8, 13.8), (18.8, 11.6), (22.2, 13.0)]:
    c.circle(x + px, yc + py, 0.6, YELLOW, YELLOW_LINE, sw=0.25)
c.line([(x + 17.4, yc + 13.2), (x + 22.8, yc + 13.2)], sw=0.35)
c.arrow((x + 15.6, yc + 21), (x + 15.6, yc + 14.4), sw=0.3)
c.text(x + 16.4, yc + 18.4, "β_{s}", size=SMALL, anchor="start", italic=True)
c.text(x + 10.8, yc + 26.2, "Other", size=SMALL)
c.text(x + 20.2, yc + 26.2, "Top 1%", size=SMALL)
rx0, ry0, rs = x + 28, yc + 11, 11
YS = ["#ffffff", "#fdf4dc", "#f8e3ad", "#f1cd78", "#e4b54c", "#c99a35"]
blocks = [0, 0, 0, 1, 1, 2]
for i in range(6):
    for j in range(6):
        v = 5 if i == j else ({1: 3, 2: 2}.get(abs(i - j), 1) if blocks[i] == blocks[j] else 0)
        c.rect(rx0 + j * rs / 6, ry0 + i * rs / 6, rs / 6, rs / 6, fill=YS[v], stroke="white", sw=0.15)
c.rect(rx0, ry0, rs, rs, stroke=YELLOW_LINE, sw=0.3)
c.text(rx0 + rs / 2, ry0 - 0.8, "Genes", size=SMALL)
c.text(rx0 - 0.8, ry0 + rs / 2, "Genes", size=SMALL, rotate=-90)
c.text(rx0 + rs / 2, yc + 25, "R (gene-gene)", size=SMALL, italic=False)
c.text(x + bw / 2, yc + 28.4, "m = β_{0} + sβ_{s} + ε,  ε ~ MVN(0, σ^{2}R)", size=LABEL, italic=True)

x = xs[3]
TRAITC, DRUGC = "#ecc75c", "#86bb73"
c.rect(x + 4, yc + 11.6, 1.8, 1.8, fill=TRAITC, stroke=YELLOW_LINE, sw=0.2)
c.text(x + 6.4, yc + 13.1, "Trait", size=SMALL, anchor="start")
c.rect(x + 4, yc + 15.6, 1.8, 1.8, fill=DRUGC, stroke=GREEN_LINE, sw=0.2)
c.text(x + 6.4, yc + 17.1, "Drug", size=SMALL, anchor="start")
base = yc + 18
vals = [(5, -4), (-3, 3.4), (4.2, -3.6), (0.9, 0.4), (-4, 3.8), (2.4, -2.2)]
for k, (t, d) in enumerate(vals):
    bx0 = x + 14 + k * 4.4
    c.rect(bx0, base - max(t, 0), 1.6, abs(t), fill=TRAITC, stroke=YELLOW_LINE, sw=0.2)
    c.rect(bx0 + 1.8, base - max(d, 0), 1.6, abs(d), fill=DRUGC, stroke=GREEN_LINE, sw=0.2)
c.line([(x + 13, base), (x + 40, base)], sw=0.3)
c.text(x + 27, yc + 25, "Consensus modules", size=SMALL)
c.text(x + bw / 2, yc + 28.4, "score = −(drug · trait)", size=LABEL, italic=True)


# ---------------- d ----------------
c.tag(1, 154, "d")
yd = 151
c.rect(4, yd, 174, 31.5, fill="#f6f7f8", stroke="#aab3ba", sw=0.35, rx=1)
c.add('<g transform="translate(0,1.8)">')
yd = 151
CR, CR_FILL = "#c25b86", "#f1cadb"
c.text(22, yd + 3.4, "Lipid CRISPR-Cas9 screen", size=TITLE, weight="bold")
c.icon(IC / "liver.svg", 17.5, yd + 5.6, 9, 5.6)
wbox(7, yd + 13.5, 30, 4, ["Lipid-increasing genes"])
wbox(7, yd + 19, 30, 4, ["Lipid-decreasing genes"])
c.text(22, yd + 26.6, "HepG2 knockouts", size=SMALL)
c.arrow((39, yd + 17), (43.6, yd + 17), sw=0.3)
c.text(58, yd + 3.4, "Enrichment in modules", size=TITLE, weight="bold")
pts = [(45, 15), (47, 11.6), (49, 12.8), (52, 9.2), (55, 10.6), (59, 8.4), (63, 10.2), (67, 11.6), (71, 13.2)]
c.line([(a, yd + b) for a, b in pts], sw=0.45, stroke=TEAL)
c.line([(45, yd + 16), (71, yd + 16)], sw=0.2, stroke=GREY_LINE)
for t in [46, 47.4, 48.4, 50.2, 53.6, 58.0]:
    c.line([(t, yd + 17), (t, yd + 20)], sw=0.35, stroke=INK)
c.text(58, yd + 23.6, "Genes ranked by module loading", size=SMALL)
c.arrow((73, yd + 17), (77.6, yd + 17), sw=0.3)
c.text(91, yd + 3.4, "Characterize hit modules", size=TITLE, weight="bold")
for k, lab in enumerate(["GLS trait association", "Lipid pathways (ORA)", "Drug perturbations", "Active sample types"]):
    wbox(79, yd + 6.4 + k * 5, 24, 4, [lab])
c.arrow((105, yd + 17), (109.6, yd + 17), sw=0.3)
c.text(146, yd + 3.4, "Recovery beyond single-gene tests", size=TITLE, weight="bold")
mx, my = 121, yd + 16.5
MOD, MOD_TINT = "#2a9686", "#d3e9e6"
c.circle(mx, my, 10.2, MOD_TINT)
import random as _r
_rng = _r.Random(7)
pts_m = {"cr": (115.2, my + 0.4), "t1": (126.2, my - 4.6), "t2": (126.8, my + 4.8)}
others = [(118.0, my - 6.4), (121.6, my - 7.6), (124.8, my - 1.0), (119.6, my - 2.0), (117.6, my + 5.6), (121.4, my + 3.4), (122.8, my + 8.0), (113.6, my - 3.8)]
allp = list(pts_m.values()) + others
for i in range(len(allp)):
    for j in range(i + 1, len(allp)):
        if _rng.random() < 0.28:
            c.line([allp[i], allp[j]], sw=0.18, stroke="#9aa4ac", dash="0.5,0.45")
for x, y in others:
    c.circle(x, y, 0.9, MOD, "#1f6f64", sw=0.2)
for k in ("t1", "t2"):
    c.circle(*pts_m[k], 1.6, "none", MOD, sw=0.45)
    c.circle(*pts_m[k], 0.9, MOD, "#1f6f64", sw=0.2)
x0, y0 = pts_m["cr"]
c.add(f'<polygon points="{x0},{y0 - 1.5} {x0 + 1.5},{y0} {x0},{y0 + 1.5} {x0 - 1.5},{y0}" fill="{MOD}" stroke="{INK}" stroke-width="0.3"/>')
c.rect(136, my - 2.4, 13, 4.8, fill="#eef7f6", stroke=MOD, sw=0.35, rx=0.8)
c.text(142.5, my + 0.75, "Disease", size=AXIS)
for k, ye in (("t1", my - 1.3), ("t2", my + 1.3)):
    c.arrow((pts_m[k][0] + 1.8, pts_m[k][1]), (135.4, ye), sw=0.32, stroke=MOD)
lx = 153
d = 1.1
c.add(f'<polygon points="{lx},{yd + 9.6 - d} {lx + d},{yd + 9.6} {lx},{yd + 9.6 + d} {lx - d},{yd + 9.6}" fill="{MOD}" stroke="{INK}" stroke-width="0.3"/>')
c.text(lx + 2, yd + 10.5, "CRISPR screen gene", size=SMALL, anchor="start")
c.circle(lx, yd + 14.6, 0.8, MOD, "#1f6f64", sw=0.2)
c.text(lx + 2, yd + 15.5, "Gene (top loadings)", size=SMALL, anchor="start")
c.circle(lx, yd + 19.6, 1.3, "none", MOD, sw=0.4)
c.circle(lx, yd + 19.6, 0.7, MOD, "#1f6f64", sw=0.2)
c.text(lx + 2, yd + 20.5, "TWAS-associated gene", size=SMALL, anchor="start")
c.line([(lx - 1.4, yd + 24.6), (lx + 1.4, yd + 24.6)], sw=0.25, stroke="#9aa4ac", dash="0.5,0.45")
c.text(lx + 2, yd + 25.3, "Shared enriched pathway", size=SMALL, anchor="start")

c.add("</g>")
c.save(Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "overview.svg")
