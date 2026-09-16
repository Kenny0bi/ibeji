"""Prove the figure checks fail when they should and pass when they should.

Run: python3 src/test_viz_checks.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import viz_style as vs


def bars(ax):
    ax.bar([0, 1], [3, 2], color=vs.EUR)
    ax.bar([0, 1], [1, 1], bottom=[3, 2], color=vs.EUR, alpha=0.38)


def case_gray_legend():
    """The Figure 2 bug: gray swatches for blue bars."""
    fig, ax = plt.subplots()
    bars(ax)
    ax.legend([Rectangle((0, 0), 1, 1, color=vs.INK_2), Rectangle((0, 0), 1, 1, color=vs.INK_2, alpha=0.38)],
              ["solid", "light"])
    return fig


def case_matching_legend():
    fig, ax = plt.subplots()
    bars(ax)
    ax.legend([Rectangle((0, 0), 1, 1, color=vs.EUR), Rectangle((0, 0), 1, 1, color=vs.EUR, alpha=0.38)],
              ["solid", "light"])
    return fig


def case_alpha_mismatch():
    """Right hue, wrong transparency: still a mismatch."""
    fig, ax = plt.subplots()
    bars(ax)
    ax.legend([Rectangle((0, 0), 1, 1, color=vs.EUR, alpha=0.8)], ["wrong alpha"])
    return fig


def case_label_on_neighbor():
    """A long y-label from the right panel running into the left panel."""
    fig = plt.figure(figsize=(4, 2))
    a = fig.add_axes([0.1, 0.15, 0.4, 0.75])
    b = fig.add_axes([0.52, 0.15, 0.4, 0.75])
    a.plot([0, 1], [0, 1])
    b.plot([0, 1], [1, 0])
    b.set_ylabel("a deliberately long label that runs left", fontsize=11)
    return fig


def case_clean_panels():
    fig = plt.figure(figsize=(6, 2))
    a = fig.add_axes([0.1, 0.2, 0.35, 0.7])
    b = fig.add_axes([0.6, 0.2, 0.35, 0.7])
    a.plot([0, 1], [0, 1]); b.plot([0, 1], [1, 0]); b.set_ylabel("short")
    return fig


def case_overlapping_text():
    fig, ax = plt.subplots()
    ax.text(0.5, 0.5, "first label", fontsize=12)
    ax.text(0.52, 0.5, "second label", fontsize=12)
    return fig


def case_rounded_ticks():
    """The Figure 12 bug: ticks every 0.25 formatted to one decimal."""
    fig, ax = plt.subplots()
    ax.plot([0, 1], [-1.0, 0.25])
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.25))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.1f}"))
    return fig


def case_mirrored_abs_ticks():
    """Mirrored axes that print |value| are intentional and must pass."""
    fig, ax = plt.subplots()
    ax.plot([0, 1], [-0.4, 0.4])
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v):.1f}"))
    return fig


def case_sqrt_axis_correct():
    """Positions at sqrt(value), labels state the value, transform declared: must pass."""
    import numpy as np
    fig, ax = plt.subplots()
    vals = [0, 0.01, 0.1, 0.4]
    ax.set_xticks(np.sqrt(vals)); ax.set_xticklabels([f"{v:g}" for v in vals]); ax.set_xlim(0, 0.8)
    ax.xaxis.ibeji_value_of = np.square
    return fig


def case_sqrt_axis_wrong_label():
    """Same axis, but one label states the wrong value: must be caught."""
    import numpy as np
    fig, ax = plt.subplots()
    vals = [0, 0.01, 0.1, 0.4]
    ax.set_xticks(np.sqrt(vals)); ax.set_xticklabels(["0", "0.01", "0.2", "0.4"]); ax.set_xlim(0, 0.8)
    ax.xaxis.ibeji_value_of = np.square
    return fig


def case_tiny_text():
    """Text below the 7 pt print minimum must be caught."""
    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    ax.plot([0, 1], [0, 1])
    ax.text(0.2, 0.8, "too small to print", fontsize=5.5)
    return fig


EXPECT = {
    case_tiny_text: "text below",
    case_sqrt_axis_correct: None,
    case_sqrt_axis_wrong_label: "tick label misstates value",
    case_rounded_ticks: "tick label misstates value",
    case_mirrored_abs_ticks: None,
    case_gray_legend: "legend color not in plot",
    case_matching_legend: None,
    case_alpha_mismatch: "legend color not in plot",
    case_label_on_neighbor: "text on another panel",
    case_clean_panels: None,
    case_overlapping_text: "overlap",
}

failures = 0
for make, expected in EXPECT.items():
    fig = make()
    problems = vs.check_layout(fig)
    plt.close(fig)
    hit = expected is None and not problems or expected is not None and any(p.startswith(expected) for p in problems)
    failures += not hit
    print(f"{'PASS' if hit else 'FAIL'}  {make.__name__:24s} expected={expected!s:28s} got={problems}")
print("all checks behave correctly" if failures == 0 else f"{failures} check(s) misbehaving")
raise SystemExit(failures)
