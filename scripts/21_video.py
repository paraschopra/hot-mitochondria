"""Plain-language explainer video (~58 s, 1080x1080, burned-in captions, original score from 20_score.py).

Scenes (seconds):
  0-8    hook: mitochondria are the cell's power plants; claimed ~50 °C, >10° hotter than the rest of the cell
  8-17   proteins record heat: 13 species, living temperature vs heat their proteins can take (real Meltome medians)
  17-30  the test: expected 5-12° tougher vs found ~0-1° (pooled lysate / intact-cell estimates and their CIs)
  30-37  two more checks: same-family "sister" proteins, AI sequence reader
  37-48  alternative: the whole cell is just as hot? needs ~47-51 °C; published cytosol readings rise 1-2 °C
  48-58  one experiment would settle it; end card

Writes video/hot_mitochondria.mp4 (+ video/captions.srt). Usage: python scripts/21_video.py [WORKERS]
"""
import os
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mpl"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors as mcolors
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch, Rectangle

OUT = ROOT / "video"
W = H = 1080
DPI = 100
FPS = 30
DURATION = 58.0
REPO_URL = "github.com/paraschopra/hot-mitochondria"

BG, FG, DIM, FAINT = "#0f1216", "#f3efe7", "#9aa1ab", "#2a3038"
HOT, WARM, COOL, GOOD = "#ff6a45", "#ffb15c", "#5aa9ff", "#5fd39a"
FONT = "Avenir Next"
plt.rcParams["font.family"] = FONT

# 13 species, lysate, 0.5% NP40 (Meltome Atlas; same table as Fig. 1A): growth temperature, median protein Tm
SPECIES = [
    (15, 41.79, "Antarctic\nsea microbe"), (20, 41.16, None), (25, 45.18, None), (28, 50.90, None), (28, 42.94, None),
    (30, 44.01, None), (30, 49.50, None), (37, 49.69, None), (37, 51.67, "Humans"), (37, 50.02, None),
    (55, 70.27, None), (60, 71.57, None), (70, 84.67, "Hot-spring\nmicrobe"),
]
FIT = np.polyfit([s[0] for s in SPECIES], [s[1] for s in SPECIES], 1)

CAPTIONS = [
    (0.3, 2.8, "Mitochondria are the tiny power plants inside our cells."),
    (2.8, 5.4, "Some scientists think they run at about 50 °C,"),
    (5.4, 8.0, "more than 10 degrees hotter than the rest of the cell."),
    (8.0, 12.5, "Living things from hotter places build proteins that hold up to more heat."),
    (12.5, 17.0, "Across 13 species, from icy seas to hot springs, the pattern holds."),
    (17.0, 21.5, "So proteins working inside hot power plants should be built much tougher."),
    (21.5, 25.5, "We checked public data: over 100,000 measurements, 6 species, 4 labs."),
    (25.5, 30.0, "They're not tougher. They hold up about the same as the rest of the cell."),
    (30.0, 32.6, "Two more checks agreed:"),
    (32.6, 37.0, "matched “sister” proteins, and an AI that reads heat-toughness from a protein's recipe."),
    (37.0, 42.0, "Another possibility: maybe the whole cell gets just as hot, so all its proteins had to toughen up."),
    (42.0, 48.0, "That seems unlikely. The few measurements so far show the rest of the cell warming only 1–2 degrees."),
    (48.0, 53.5, "One simple lab test could settle it: switch the power plants off and see if their proteins suddenly look tougher."),
    (53.5, 57.6, "Summary and code: link in the post."),
]


# ---------------------------------------------------------------- helpers
def smooth(x):
    x = np.clip(x, 0, 1)
    return x * x * (3 - 2 * x)


def ramp(t, a, b):
    return smooth((t - a) / (b - a))


def window(t, a, b, fade=0.45):
    return min(ramp(t, a, a + fade), 1 - ramp(t, b - fade, b))


def mix(c1, c2, f):
    a, b = np.array(mcolors.to_rgb(c1)), np.array(mcolors.to_rgb(c2))
    return tuple(a + (b - a) * float(np.clip(f, 0, 1)))


def wrap(text, width=38):
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width and cur:
            lines.append(cur); cur = w
        else:
            cur = (cur + " " + w).strip()
    lines.append(cur)
    return "\n".join(lines)


def txt(ax, x, y, s, size=22, color=FG, alpha=1.0, weight="regular", ha="center", va="center", **kw):
    if alpha <= 0.01:
        return
    ax.text(x, y, s, fontsize=size, color=color, alpha=alpha, fontweight=weight, ha=ha, va=va, **kw)


def mito(ax, x, y, w, h, angle, heat, alpha):
    """A bean-shaped mitochondrion with folded inner membrane; heat 0..1 shifts colour cool->hot and adds glow."""
    col = mix(COOL, HOT, heat)
    for k in range(4):  # glow
        g = 1 + 0.18 * (k + 1)
        ax.add_patch(Ellipse((x, y), w * g, h * g * 1.15, angle=angle, color=col, alpha=alpha * heat * 0.07, lw=0))
    ax.add_patch(Ellipse((x, y), w, h, angle=angle, facecolor=mix(BG, col, 0.35), edgecolor=col, lw=3, alpha=alpha))
    u = np.linspace(-0.40, 0.40, 160)
    v = 0.30 * np.sin(u * np.pi * 9) * np.sqrt(np.clip(1 - (u / 0.45) ** 2, 0, 1))
    th = np.radians(angle)
    px = x + (u * w) * np.cos(th) - (v * h) * np.sin(th)
    py = y + (u * w) * np.sin(th) + (v * h) * np.cos(th)
    ax.plot(px, py, color=col, lw=2, alpha=alpha * 0.9, solid_capstyle="round")


MITOS = [(0.33, 0.66, 0.13, 0.055, 25), (0.64, 0.72, 0.12, 0.05, -20), (0.70, 0.50, 0.13, 0.055, 70),
         (0.36, 0.45, 0.12, 0.05, -35), (0.53, 0.38, 0.11, 0.048, 10)]


def cell(ax, alpha, cyto_heat=0.0, mito_heat=1.0, pulse=0.0, cx=0.5, cy=0.56, r=0.29):
    fill = mix("#1b2430", "#5a2418", cyto_heat)
    ax.add_patch(Ellipse((cx, cy), 2 * r * 1.08, 2 * r, facecolor=fill, edgecolor=mix(DIM, HOT, cyto_heat), lw=3,
                         alpha=alpha))
    k = r / 0.29
    ax.add_patch(Circle((cx + 0.02 * k, cy + 0.03 * k), 0.065 * k, facecolor=mix("#262f3b", "#6b2a1c", cyto_heat),
                        edgecolor=DIM, lw=2, alpha=alpha * 0.9))
    for (x, y, w, h, a) in MITOS:
        mito(ax, cx + (x - 0.5) * k, cy + (y - 0.56) * k, w * k, h * k, a, mito_heat * (0.85 + 0.15 * pulse), alpha)


def thermo(ax, x, y0, y1, frac, color, alpha, label_top=None):
    """Vertical thermometer from y0 to y1, filled to frac."""
    wdt = 0.05
    ax.add_patch(FancyBboxPatch((x - wdt / 2, y0), wdt, y1 - y0, boxstyle="round,pad=0,rounding_size=0.025",
                                facecolor=FAINT, edgecolor=DIM, lw=2, alpha=alpha))
    ax.add_patch(Circle((x, y0), 0.045, facecolor=color, edgecolor=DIM, lw=2, alpha=alpha))
    fh = (y1 - y0 - 0.02) * np.clip(frac, 0, 1)
    ax.add_patch(Rectangle((x - wdt / 2 + 0.011, y0), wdt - 0.022, fh, color=color, alpha=alpha, lw=0))


# ---------------------------------------------------------------- scenes
def scene_hook(ax, t, a):
    heat = ramp(t, 2.6, 4.6)
    cell(ax, a, cyto_heat=0.0, mito_heat=heat, pulse=0.5 + 0.5 * np.sin(t * 3))
    txt(ax, 0.5, 0.92, "Are mitochondria running a fever?", 30, FG, a * ramp(t, 0.2, 1.0), "bold")
    al = a * ramp(t, 3.0, 3.8)
    txt(ax, 0.84, 0.77, "~50 °C?", 34, HOT, al, "bold")
    ax.plot([0.69, 0.78], [0.735, 0.765], color=HOT, lw=2, alpha=al)
    al2 = a * ramp(t, 5.4, 6.2)
    txt(ax, 0.17, 0.25, "rest of the cell\n~37 °C", 22, COOL, al2, "demibold")


def scene_species(ax, t, a):
    t0 = 8.0
    x0, x1, y0, y1 = 0.15, 0.90, 0.33, 0.82
    X = lambda v: x0 + (v - 10) / (75 - 10) * (x1 - x0)
    Y = lambda v: y0 + (v - 35) / (90 - 35) * (y1 - y0)
    txt(ax, 0.5, 0.92, "Proteins remember the heat they live in", 28, FG, a, "bold")
    ax.plot([x0, x1], [y0, y0], color=DIM, lw=2, alpha=a)
    ax.plot([x0, x0], [y0, y1], color=DIM, lw=2, alpha=a)
    for v in (20, 40, 60):
        txt(ax, X(v), y0 - 0.025, f"{v}°", 16, DIM, a)
    for v in (40, 60, 80):
        txt(ax, x0 - 0.03, Y(v), f"{v}°", 16, DIM, a, ha="right")
    txt(ax, (x0 + x1) / 2, y0 - 0.07, "Temperature it lives at", 19, DIM, a)
    txt(ax, x0 - 0.10, (y0 + y1) / 2, "Heat its proteins can take", 19, DIM, a, rotation=90)
    order = sorted(range(len(SPECIES)), key=lambda i: SPECIES[i][0])
    for k, i in enumerate(order):
        g, tm, lab = SPECIES[i]
        al = a * ramp(t, t0 + 0.6 + 0.28 * k, t0 + 1.0 + 0.28 * k)
        col = mix(COOL, HOT, (g - 15) / 55)
        ax.add_patch(Circle((X(g), Y(tm)), 0.014, color=col, alpha=al, lw=0))
        if lab:
            dx, dy, ha = {"Humans": (-0.02, 0.06, "center"), "Hot-spring\nmicrobe": (-0.02, -0.075, "center"),
                          "Antarctic\nsea microbe": (0.03, 0.075, "center")}[lab]
            txt(ax, X(g) + dx, Y(tm) + dy, lab, 17, FG, al, "demibold", ha=ha)
    lf = ramp(t, t0 + 4.8, t0 + 6.5)
    if lf > 0:
        gx = np.linspace(12, 12 + (73 - 12) * lf, 50)
        ax.plot([X(v) for v in gx], [Y(np.polyval(FIT, v)) for v in gx], color=WARM, lw=3, alpha=a * 0.8)


def scene_test(ax, t, a):
    t0 = 17.0
    txt(ax, 0.5, 0.92, "How much tougher are power-plant proteins?", 27, FG, a, "bold")
    lo, hi = -3, 14
    x0, x1, yl = 0.10, 0.90, 0.52
    X = lambda v: x0 + (v - lo) / (hi - lo) * (x1 - x0)
    ax.plot([x0, x1], [yl, yl], color=DIM, lw=2.5, alpha=a)
    for v in (0, 5, 10):
        ax.plot([X(v), X(v)], [yl - 0.012, yl + 0.012], color=DIM, lw=2, alpha=a)
        txt(ax, X(v), yl - 0.045, f"+{v}°" if v else "0°", 18, DIM, a)
    txt(ax, X(0), yl - 0.09, "same as the\nrest of the cell", 15, DIM, a)
    # expected band (+5.5 to +12)
    ab = a * ramp(t, t0 + 1.0, t0 + 2.2)
    ax.add_patch(FancyBboxPatch((X(5.5), yl - 0.03), X(12) - X(5.5), 0.06, boxstyle="round,pad=0,rounding_size=0.012",
                                color=HOT, alpha=ab * 0.55, lw=0))
    txt(ax, (X(5.5) + X(12)) / 2, yl + 0.10, "Expected if they ran\n10–15° hotter", 20, HOT, ab, "demibold")
    # data badges
    ad = a * ramp(t, t0 + 4.6, t0 + 5.4)
    for i, (big, small) in enumerate((("100,000+", "measurements"), ("6", "species"), ("4", "labs"))):
        cx = 0.24 + 0.26 * i
        txt(ax, cx, 0.78, big, 30, FG, ad, "bold")
        txt(ax, cx, 0.735, small, 17, DIM, ad)
    # found: lysate -0.2 [-1.0, +0.6]; intact cells +1.3 [+0.7, +1.8] -> span -1.0..+1.8
    af = a * ramp(t, t0 + 8.6, t0 + 9.4)
    drop = 1 - ramp(t, t0 + 8.6, t0 + 9.4)
    yy = yl + 0.12 * drop
    ax.add_patch(FancyBboxPatch((X(-1.0), yy - 0.03), X(1.8) - X(-1.0), 0.06,
                                boxstyle="round,pad=0,rounding_size=0.012", color=GOOD, alpha=af * 0.8, lw=0))
    txt(ax, X(0.4), yl - 0.165, "What we found:\nabout the same", 21, GOOD, af, "bold")


def check_card(ax, x, y, title, body, a):
    ax.add_patch(FancyBboxPatch((x - 0.2, y - 0.17), 0.4, 0.36, boxstyle="round,pad=0,rounding_size=0.03",
                                facecolor="#171c22", edgecolor=FAINT, lw=2, alpha=a))
    ax.add_patch(Circle((x, y + 0.10), 0.045, color=GOOD, alpha=a, lw=0))
    ax.plot([x - 0.02, x - 0.004, x + 0.024], [y + 0.10, y + 0.082, y + 0.118], color=BG, lw=5, alpha=a,
            solid_capstyle="round")
    txt(ax, x, y + 0.01, title, 22, FG, a, "bold")
    txt(ax, x, y - 0.085, body, 16, DIM, a, linespacing=1.35)


def scene_checks(ax, t, a):
    t0 = 30.0
    txt(ax, 0.5, 0.92, "Two more checks", 30, FG, a, "bold")
    check_card(ax, 0.28, 0.58, "Sister proteins",
               "47 matched pairs:\none copy inside the\npower plant, one outside.\nNo difference.",
               a * ramp(t, t0 + 1.0, t0 + 1.8))
    check_card(ax, 0.72, 0.58, "AI protein reader",
               "Reads heat-toughness\nfrom a protein's recipe.\nNo sign of extra heat,\nin 5 of 5 species.",
               a * ramp(t, t0 + 2.8, t0 + 3.6))


def scene_alt(ax, t, a):
    t0 = 37.0
    split = ramp(t, t0 + 5.0, t0 + 5.6)  # cell fades out, then thermometers fade in
    ac = a * (1 - ramp(t, t0 + 4.4, t0 + 5.0))
    if ac > 0.01:
        cell(ax, ac, cyto_heat=ramp(t, t0 + 0.8, t0 + 2.6), mito_heat=1.0, pulse=0.5 + 0.5 * np.sin(t * 3))
        txt(ax, 0.5, 0.92, "Another possibility", 30, FG, ac, "bold")
    at = a * split
    if at > 0.01:
        txt(ax, 0.5, 0.92, "Would the rest of the cell be ~50 °C?", 28, FG, at, "bold")
        y0, y1 = 0.42, 0.80
        # scale 35..55 °C
        fr = lambda v: (v - 35) / (55 - 35)
        grow = ramp(t, t0 + 5.6, t0 + 7.0)
        thermo(ax, 0.32, y0, y1, fr(37 + (49 - 37) * grow), HOT, at)
        thermo(ax, 0.68, y0, y1, fr(37 + (38.5 - 37) * grow), COOL, at)
        txt(ax, 0.32, y0 - 0.075, "Needed for this idea", 20, HOT, at, "demibold")
        txt(ax, 0.32, y0 - 0.12, "~47–51 °C", 24, FG, at, "bold")
        txt(ax, 0.68, y0 - 0.075, "Measured so far", 20, COOL, at, "demibold")
        txt(ax, 0.68, y0 - 0.12, "warms 1–2 °C", 24, FG, at, "bold")
        txt(ax, 0.5, 0.86, "Brief, local spikes haven't been ruled out.", 16, DIM, at * ramp(t, t0 + 7.0, t0 + 7.8),
            style="italic")


def scene_next(ax, t, a):
    t0 = 48.0
    off = ramp(t, t0 + 2.2, t0 + 3.6)
    cell(ax, a, cyto_heat=0.0, mito_heat=1 - off, pulse=(0.5 + 0.5 * np.sin(t * 3)) * (1 - off), cy=0.60, r=0.24)
    txt(ax, 0.5, 0.92, "One experiment could settle it", 30, FG, a, "bold")
    # switch
    sx, sy = 0.5, 0.285
    ax.add_patch(FancyBboxPatch((sx - 0.07, sy - 0.03), 0.14, 0.06, boxstyle="round,pad=0,rounding_size=0.03",
                                facecolor=mix(HOT, FAINT, off), edgecolor=DIM, lw=2, alpha=a))
    ax.add_patch(Circle((sx + 0.04 - 0.08 * off, sy), 0.024, color=FG, alpha=a, lw=0))
    txt(ax, sx - 0.13, sy, "ON", 18, HOT, a * (1 - off), "bold", ha="right")
    txt(ax, sx + 0.13, sy, "OFF", 18, COOL, a * off, "bold", ha="left")


def scene_end(ax, t, a):
    txt(ax, 0.5, 0.72, "Are mitochondria\nreally 50 °C?", 44, FG, a, "bold", linespacing=1.1)
    txt(ax, 0.5, 0.55, "Their proteins show no sign of it.", 27, WARM, a * ramp(t, 54.2, 55.0), "demibold")
    txt(ax, 0.5, 0.42, REPO_URL, 22, FG, a * ramp(t, 54.8, 55.6))
    txt(ax, 0.5, 0.36, "A fresh look at public data, with AI assistance.\nNot a new lab measurement.", 16, DIM,
        a * ramp(t, 54.8, 55.6), linespacing=1.4)


SCENES = [(0.0, 8.0, scene_hook), (8.0, 17.0, scene_species), (17.0, 30.0, scene_test), (30.0, 37.0, scene_checks),
          (37.0, 48.0, scene_alt), (48.0, 53.5, scene_next), (53.5, DURATION, scene_end)]


def draw_caption(ax, t):
    for a0, b0, s in CAPTIONS:
        al = window(t, a0, b0, 0.18)
        if al > 0.01:
            w = wrap(s, 40)
            ax.text(0.5, 0.055, w, fontsize=25, color=FG, alpha=al, ha="center", va="bottom", fontweight="demibold",
                    linespacing=1.3, bbox=dict(boxstyle="round,pad=0.45,rounding_size=0.6", facecolor="black",
                                               alpha=0.55 * al, edgecolor="none"))


def frame(t, fig):
    fig.clf()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 1, 1, color=BG, zorder=-10))
    for a0, b0, fn in SCENES:
        if a0 - 0.01 <= t <= b0 + 0.01:
            fade_in = 0.0 if a0 == 0 else 0.45
            al = min(ramp(t, a0, a0 + fade_in) if fade_in else 1.0, 1 - ramp(t, b0 - 0.45, b0))
            if b0 >= DURATION:
                al = min(ramp(t, a0, a0 + 0.45), 1 - ramp(t, DURATION - 0.8, DURATION))
            if al > 0.005:
                fn(ax, t, al)
    draw_caption(ax, t)
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba())[:, :, :3].tobytes()


def render_chunk(args):
    i, f0, f1 = args
    path = OUT / f"_chunk{i:02d}.mp4"
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    p = subprocess.Popen(["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}",
                          "-r", str(FPS), "-i", "-", "-c:v", "libx264", "-preset", "medium", "-crf", "18",
                          "-pix_fmt", "yuv420p", str(path)], stdin=subprocess.PIPE)
    for f in range(f0, f1):
        p.stdin.write(frame(f / FPS, fig))
    p.stdin.close(); p.wait()
    plt.close(fig)
    return path


def srt_time(s):
    ms = int(round(s * 1000))
    return f"{ms // 3600000:02d}:{ms // 60000 % 60:02d}:{ms // 1000 % 60:02d},{ms % 1000:03d}"


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    if len(sys.argv) > 1 and sys.argv[1] == "still":  # quick preview frames: python scripts/21_video.py still 5 20 ...
        fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
        for ts in sys.argv[2:]:
            frame(float(ts), fig)
            fig.savefig(OUT / f"_still_{float(ts):05.1f}.png", dpi=DPI)
        sys.exit()
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else max(2, os.cpu_count() - 2)
    n = int(DURATION * FPS)
    edges = np.linspace(0, n, workers + 1).astype(int)
    with Pool(workers) as pool:
        chunks = pool.map(render_chunk, [(i, edges[i], edges[i + 1]) for i in range(workers)])
    lst = OUT / "_chunks.txt"
    lst.write_text("".join(f"file '{c.name}'\n" for c in chunks))
    subprocess.run([sys.executable, str(ROOT / "scripts" / "20_score.py"), str(DURATION)], check=True)
    final = OUT / "hot_mitochondria.mp4"
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", str(lst),
                    "-i", str(OUT / "score.wav"), "-c:v", "copy", "-c:a", "aac", "-b:a", "192k", "-shortest",
                    "-movflags", "+faststart", str(final)], check=True)
    for c in chunks:
        c.unlink()
    lst.unlink()
    (OUT / "captions.srt").write_text("".join(
        f"{k + 1}\n{srt_time(a0)} --> {srt_time(b0)}\n{s}\n\n" for k, (a0, b0, s) in enumerate(CAPTIONS)))
    print("wrote", final)
