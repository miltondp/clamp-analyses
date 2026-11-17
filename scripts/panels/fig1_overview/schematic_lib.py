"""Shared helpers for the hand-built schematics (Figure 1, Figure 4a projection scheme).

All coordinates are millimetres at final print size; font sizes are given in points and
converted, so the figures meet Nature Methods rules directly (Helvetica/Arial, 5-7 pt text,
8 pt bold lowercase panel letters).
"""
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

PT = 25.4 / 72  # mm per point
FONT = "Helvetica, Arial, 'Liberation Sans', sans-serif"
TAG, TITLE, LABEL, AXIS, SMALL = 8, 7, 6, 5.5, 5

INK = "#20262d"
BLUE, BLUE_LINE = "#bfd8ed", "#729ebf"
GREY, GREY_LINE = "#e4e8ea", "#9ba7ae"
YELLOW, YELLOW_LINE = "#f7d984", "#b38d35"
GREEN, GREEN_LINE = "#b5d6a0", "#649157"
TEAL, TEAL_FILL = "#208b82", "#aad3d0"
PURPLE, PURPLE_FILL = "#5e4b9c", "#b7aad3"
ORANGE, ORANGE_FILL = "#d9785a", "#f4c3b0"
STUDY, CONTROL = "#347ba5", "#9ea9b2"
SCREEN, SCREEN_LINE = "#eff5ea", "#9bb28d"


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


class Canvas:
    def __init__(self, w, h, title):
        self.w, self.h, self.title, self.items = w, h, title, []
        self.n_icon = 0

    def add(self, s):
        self.items.append(s)

    def text(self, x, y, s, size=LABEL, weight="normal", anchor="middle", fill=INK,
             italic=False, rotate=0, halo=None):
        """s may contain _{sub} and ^{sup} markup; baseline at y."""
        style = f"font-size:{size * PT:.3f}px;font-weight:{weight};"
        if italic:
            style += "font-style:italic;"
        if halo:  # halo drawn as its own layer so it never covers neighbouring glyphs
            self.text(x, y, s, size, weight, anchor, halo, italic, rotate)
            self.items[-1] = self.items[-1].replace('style="', f'style="stroke:{halo};stroke-width:{size * PT * 0.5:.3f}px;stroke-linejoin:round;', 1)
        body = ""
        for part in re.split(r"(_\{[^}]*\}|\^\{[^}]*\})", s):
            if part.startswith("_{"):
                body += f'<tspan dy="{size * PT * 0.25:.3f}" font-size="{size * PT * 0.72:.3f}">{esc(part[2:-1])}</tspan><tspan dy="{-size * PT * 0.25:.3f}"></tspan>'
            elif part.startswith("^{"):
                body += f'<tspan dy="{-size * PT * 0.35:.3f}" font-size="{size * PT * 0.72:.3f}">{esc(part[2:-1])}</tspan><tspan dy="{size * PT * 0.35:.3f}"></tspan>'
            elif part:
                body += esc(part)
        tr = f' transform="rotate({rotate} {x:.3f} {y:.3f})"' if rotate else ""
        self.add(f'<text x="{x:.3f}" y="{y:.3f}" text-anchor="{anchor}" fill="{fill}" style="{style}"{tr}>{body}</text>')

    def lines(self, x, y, rows, size=LABEL, lh=1.18, **kw):
        for i, r in enumerate(rows):
            self.text(x, y + i * size * PT * lh, r, size=size, **kw)

    def rect(self, x, y, w, h, fill="none", stroke=INK, sw=0.25, rx=0, dash=None, opacity=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        o = f' fill-opacity="{opacity}"' if opacity is not None else ""
        self.add(f'<rect x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" height="{h:.3f}" rx="{rx}" fill="{fill}"{o} stroke="{stroke}" stroke-width="{sw}"{d}/>')

    def circle(self, x, y, r, fill, stroke=None, sw=0.2):
        s = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
        self.add(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{r:.3f}" fill="{fill}"{s}/>')

    def line(self, pts, sw=0.3, stroke=INK, arrow=False, dash=None):
        p = " ".join(f"{x:.3f},{y:.3f}" for x, y in pts)
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.add(f'<polyline points="{p}" fill="none" stroke="{stroke}" stroke-width="{sw}" stroke-linejoin="round"{d}/>')
        if arrow:
            self.arrowhead(pts[-2], pts[-1], sw, stroke)

    def arrowhead(self, p0, p1, sw, color=INK):
        """Filled triangle at p1 pointing from p0 to p1.

        Drawn as a plain polygon rather than an SVG marker, because WeasyPrint
        (Manubot's PDF renderer) misplaces and misorients markers. Geometry
        matches the former marker: 5.5 stroke widths long, tip 0.6 units past p1.
        """
        (x0, y0), (x1, y1) = p0, p1
        length = math.hypot(x1 - x0, y1 - y0) or 1.0
        ux, uy = (x1 - x0) / length, (y1 - y0) / length
        k = sw * 5.5 / 6
        tip = (x1 + 0.6 * k * ux, y1 + 0.6 * k * uy)
        back = (x1 - 5.4 * k * ux, y1 - 5.4 * k * uy)
        half = 2.6 * k
        pts = [tip, (back[0] - half * uy, back[1] + half * ux), (back[0] + half * uy, back[1] - half * ux)]
        q = " ".join(f"{x:.3f},{y:.3f}" for x, y in pts)
        self.add(f'<polygon points="{q}" fill="{color}"/>')

    def arrow(self, *pts, **kw):
        self.line(pts, arrow=True, **kw)

    def bracket(self, x0, x1, y, depth=1.6, tip=None):
        """Under-brace joining x0..x1 at y; optional arrow from its centre down to tip."""
        self.line([(x0, y), (x0, y + depth), (x1, y + depth), (x1, y)])
        if tip:
            self.arrow(((x0 + x1) / 2, y + depth), *tip)

    def op(self, x, y, s):
        """Operator (approx., times) vertically centred on y."""
        self.text(x, y + 9 * PT * 0.32, s, size=9)

    def tag(self, x, y, letter):
        self.text(x, y, letter, size=TAG, weight="bold", anchor="start")

    def matrix(self, x, y, w, h, fill=BLUE, stroke=BLUE_LINE, cell=2.2, sym=None, name=None,
               top=None, left=None, bottom=None, title=None, hl_col=None, hl_row=None, sym_size=TITLE,
               ncols=None, nrows=None):
        """Grid matrix with dimension labels placed at fixed offsets so they align everywhere."""
        self.rect(x, y, w, h, fill=fill, stroke="none")
        nc, nr = ncols or max(1, round(w / cell)), nrows or max(1, round(h / cell))
        g = []
        for i in range(1, nc):
            xi = x + w * i / nc
            g.append(f"M{xi:.3f},{y:.3f}V{y + h:.3f}")
        for j in range(1, nr):
            yj = y + h * j / nr
            g.append(f"M{x:.3f},{yj:.3f}H{x + w:.3f}")
        self.add(f'<path d="{"".join(g)}" stroke="{stroke}" stroke-width="0.12" fill="none"/>')
        if hl_col:
            cw = w / nc
            self.rect(x + w - cw, y, cw, h, fill=hl_col[0], stroke=hl_col[1], sw=0.3)
        if hl_row:
            rh = h / nr
            self.rect(x, y, w, rh, fill=hl_row[0], stroke=hl_row[1], sw=0.3)
        self.rect(x, y, w, h, stroke=stroke, sw=0.3)
        cy = y + h / 2
        if name and sym:
            self.text(x + w / 2, cy - 0.4, name, size=LABEL, weight="bold", halo=fill)
            self.text(x + w / 2, cy + sym_size * PT * 0.95, sym, size=sym_size, italic=True, halo=fill)
        elif name:
            self.text(x + w / 2, cy + LABEL * PT * 0.35, name, size=LABEL, weight="bold", halo=fill)
        elif sym:
            self.text(x + w / 2, cy + sym_size * PT * 0.35, sym, size=sym_size, italic=True, halo=fill)
        if top:
            self.text(x + w / 2, y - 1.0, top, size=AXIS)
        if title:
            self.text(x + w / 2, y - 1.0 - (AXIS * PT + 1.0 if top else 0), title, size=LABEL, weight="bold")
        if left:
            self.text(x - 1.0, cy, left, size=AXIS, rotate=-90)
        if bottom:
            self.text(x + w / 2, y + h + 1.0 + AXIS * PT * 0.8, bottom, size=AXIS)

    def icon(self, path, x, y, w, h, thin=0):
        """Embed an exported plain SVG icon scaled into the (x, y, w, h) box.

        thin > 0 trims filled line-art with a white outline (icon units; white backgrounds only)."""
        self.n_icon += 1
        root = ET.parse(path).getroot()
        vb = root.get("viewBox")
        ns = "{http://www.w3.org/2000/svg}"
        keep = [c for c in root if c.tag not in (ns + "defs", ns + "title", ns + "desc", ns + "metadata")
                and not c.tag.endswith("namedview")]
        inner = "".join(ET.tostring(c, encoding="unicode") for c in keep)
        inner = re.sub(r"ns\d+:|xmlns(:\w+)?=\"[^\"]*\"", "", inner)
        inner = re.sub(r'id="([^"]+)"', rf'id="i{self.n_icon}_\1"', inner)
        if thin:
            inner = re.sub(r"(fill:#208b82;[^\"]*?)stroke:none;stroke-width:[0-9.]+",
                           rf"\1stroke:#ffffff;stroke-width:{thin};stroke-linejoin:round", inner)
        self.add(f'<svg x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" height="{h:.3f}" viewBox="{vb}" '
                 f'preserveAspectRatio="xMidYMid meet" overflow="visible">{inner}</svg>')

    def save(self, path):
        head = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.w}mm" height="{self.h}mm" '
                f'viewBox="0 0 {self.w} {self.h}" font-family="{FONT}">\n<title>{esc(self.title)}</title>\n'
                f'<rect width="{self.w}" height="{self.h}" fill="white"/>\n')
        Path(path).write_text(head + "\n".join(self.items) + "\n</svg>\n")
