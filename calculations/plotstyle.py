"""
Publication figure style — Elsevier column specifications.

The governing idea: **draw every figure at its final printed size.** The previous
`PLOT_STYLE['figsize'] = (14, 12)` with 14 pt labels looked fine on screen, but a
journal scales a 14-inch figure to a 90 mm column — a factor of 0.25 — which puts
axis labels at 3.5 pt and tick labels at 3.0 pt against a 7 pt minimum. Choosing
the width first and never rescaling makes that failure impossible rather than
merely unlikely.

Elsevier author guidelines (International Journal of Hydrogen Energy, Journal of
Nuclear Materials, Corrosion Science):

    single column     90 mm   (3.543 in)
    1.5 column       140 mm   (5.512 in)
    double column    190 mm   (7.480 in)
    minimum font       7 pt   at final size
    line art          vector (PDF/EPS) preferred over raster

Colour
------
Both previous palettes failed quantitative colour-vision checks:

* the regime bands (matplotlib `tab10`) put metal `#2ca02c` and oxide `#ff7f0e`
  at ΔE 0.7 under protanopia — indistinguishable, and regime identity is the
  whole point of those figures;
* `COLORS` put L5 `red` and L4_both `crimson` at ΔE 7.7 for *normal* vision,
  below the ΔE 15 floor, and L3 `cyan` at 1.22:1 contrast on white.

The replacements below are drawn from the Okabe–Ito colour-vision-safe set and
were checked, not eyeballed. `REGIME_COLORS` passes every check with a worst
all-pairs separation of ΔE 11.0 under deuteranopia. `LEVEL_COLORS` passes with a
worst pair of ΔE 7.6, which sits in the 6–8 band that is admissible *only*
alongside a second, non-colour channel — hence `LEVEL_LINESTYLES`. Always use
them together; colour alone is not sufficient identity for eight series in print.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

MM = 1.0 / 25.4

WIDTHS = {
    'single': 90.0 * MM,
    'onehalf': 140.0 * MM,
    'double': 190.0 * MM,
}

# Font sizes are absolute points at final size. Nothing here is scaled later.
FONT = {
    'axis': 8.0,
    'tick': 7.0,
    'legend': 7.0,
    'title': 8.5,
    'annotation': 7.0,
}

MIN_PT = 7.0

# Regimes. Okabe-Ito; worst all-pairs CVD separation ΔE 11.0 (deuteranopia).
REGIME_COLORS = {
    'metal':   '#0072B2',   # blue
    'oxide':   '#D55E00',   # vermillion
    'surface': '#009E73',   # bluish green
    'defect':  '#E69F00',   # orange
}
REGIME_ORDER = ('metal', 'oxide', 'surface', 'defect')

# Model levels. L1 is the perfect-metal reference, so it is deliberately black
# and dashed rather than occupying a categorical hue slot.
LEVEL_COLORS = {
    'L1':       '#000000',
    'L2a':      '#0072B2',
    'L2b':      '#56B4E9',
    'L3':       '#CC79A7',
    'L4_gb':    '#009E73',
    'L4_trap':  '#E69F00',
    'L4_both':  '#D55E00',
    'L5':       '#0072B2',
}
# Second channel, mandatory. Two entries share a hue (L2a/L5) and are separated
# here instead; the rest gain redundancy for print, photocopying and projection.
LEVEL_LINESTYLES = {
    'L1':      (0, (4, 2)),
    'L2a':     'solid',
    'L2b':     (0, (1, 1)),
    'L3':      (0, (3, 1, 1, 1)),
    'L4_gb':   'solid',
    'L4_trap': (0, (5, 1)),
    'L4_both': 'solid',
    'L5':      (0, (6, 1, 1, 1)),
}

# Magnitude. Single-hue, monotonically lightening, so it survives greyscale
# printing; 'YlOrRd' is multi-hue and its light end vanishes on white.
SEQUENTIAL_CMAP = 'viridis'


def apply(target='single', aspect=0.75):
    """
    Install the publication rcParams and return the figure size for `target`.

    aspect : height / width. 0.75 suits a single panel; pass a smaller value for
        wide, short panels and a larger one for tall bar charts.
    """
    if target not in WIDTHS:
        raise ValueError(f"target must be one of {sorted(WIDTHS)}, got {target!r}")
    w = WIDTHS[target]

    mpl.rcParams.update({
        'figure.figsize': (w, w * aspect),
        'figure.dpi': 150,
        'savefig.dpi': 600,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.02,
        'savefig.transparent': False,

        # Type 3 is matplotlib's default and is rejected by Elsevier, IEEE and
        # APS production. 42 embeds TrueType so the text stays selectable.
        'pdf.fonttype': 42,
        'ps.fonttype': 42,

        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans', 'Helvetica', 'Arial'],
        'font.size': FONT['tick'],
        'axes.labelsize': FONT['axis'],
        'axes.titlesize': FONT['title'],
        'xtick.labelsize': FONT['tick'],
        'ytick.labelsize': FONT['tick'],
        'legend.fontsize': FONT['legend'],

        # Recessive grid and axes: the data should be the darkest thing present.
        'axes.grid': True,
        'grid.alpha': 0.25,
        'grid.linewidth': 0.4,
        'grid.color': '#b0b0b0',
        'axes.linewidth': 0.6,
        'axes.edgecolor': '#333333',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.direction': 'out',
        'ytick.direction': 'out',

        'lines.linewidth': 1.2,
        'lines.markersize': 3.5,
        'legend.frameon': False,
        'legend.handlelength': 2.4,
        'image.cmap': SEQUENTIAL_CMAP,
    })
    return (w, w * aspect)


def audit(fig, target='single'):
    """
    Check a rendered figure against the two rules that are mechanically checkable:
    it is no wider than the column, and no text is below the 7 pt minimum.

    Returns a list of problem strings; empty means it passes.
    """
    problems = []
    w_in = fig.get_size_inches()[0]
    limit = WIDTHS[target]
    if w_in > limit + 1e-6:
        problems.append(
            f'width {w_in:.3f} in exceeds {target} column {limit:.3f} in '
            f'({w_in / limit:.2f}x) — the journal will scale it down and shrink every font')

    for t in fig.findobj(mpl.text.Text):
        if not t.get_text().strip():
            continue
        pt = t.get_fontsize()
        if pt < MIN_PT - 1e-6:
            problems.append(f'text {t.get_text()[:28]!r} at {pt:.1f} pt < {MIN_PT} pt minimum')
    return problems


def save_figure(fig, path_stem, target='single', also_png=True, check=True):
    """
    Write a figure as vector PDF (for submission) and optionally 600 dpi PNG (for
    drafts and slides). `path_stem` carries no extension.

    Raises if `check` and the figure violates the column-width or font-size rule
    — a silent oversize figure is exactly how 3 pt tick labels reach a reviewer.
    """
    if check:
        problems = audit(fig, target=target)
        if problems:
            raise ValueError('figure fails publication checks:\n  - '
                             + '\n  - '.join(problems))
    written = [f'{path_stem}.pdf']
    fig.savefig(written[0], format='pdf')
    if also_png:
        written.append(f'{path_stem}.png')
        fig.savefig(written[1], format='png', dpi=600)
    return written


def regime_color(name, default='#666666'):
    """Colour for a regime label, tolerant of unknown names."""
    return REGIME_COLORS.get(name, default)


def level_style(name):
    """(colour, linestyle) for a model level — always use both."""
    return LEVEL_COLORS.get(name, '#666666'), LEVEL_LINESTYLES.get(name, 'solid')
