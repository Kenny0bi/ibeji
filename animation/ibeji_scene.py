"""ibeji: the Shapley cube behind one cross-ancestry disagreement.

Render:  manim -qh animation/ibeji_scene.py ShapleyCube

About 50 seconds, one idea. Every number on screen comes from animation/dnajb7_numbers.json,
written by src/export_anim_numbers.py, which checks it against the decomposition output.

  1. DNAJB7 has two models: European weights above an axis, Yoruba weights below (real weights).
  2. Genetic variance of predicted expression: V = w' D^1/2 R D^1/2 w, three ingredients.
  3. Each ingredient can come from either population, so the 8 mixes are the corners of a cube,
     from all-European (000) to all-Yoruba (111). Each corner holds its real log V.
  4. There are 6 ways to walk from 000 to 111, one ingredient at a time. Each walk splits the
     total gap differently.
  5. Averaging each ingredient's step over all 6 walks gives phi_w, phi_D, phi_R, and they add
     up exactly to the gap.
"""
import json
from pathlib import Path

import numpy as np
from manim import (
    ORIGIN,
    DOWN, LEFT, RIGHT, UP, UL, UR, Create, Dot, FadeIn, FadeOut, Line, MathTex, Rectangle, Scene,
    Text, Transform, VGroup, Write, config, LaggedStart, Indicate,
)

NUM = json.loads((Path(__file__).resolve().parent / "dnajb7_numbers.json").read_text())

INDIGO = "#3056c8"
CAMWOOD = "#e0632f"
PALM = "#159b77"
GOLD = "#b98500"      # darker than the print gold so it reads on white video
CORAL = "#d4508a"
INK = "#1b1b1b"
MUTED = "#8a8882"
config.background_color = "#ffffff"
Text.set_default(font="Helvetica Neue")

PLAYER_COLOR = {"w": PALM, "D": GOLD, "R": CORAL}
PLAYER_TEX = {"w": r"\phi_w", "D": r"\phi_D", "R": r"\phi_R"}


def signed(v):
    # LaTeX typesets "-" as a true minus and rejects U+2212, so MathTex strings keep the hyphen.
    return f"{v:+.2f}"


class ShapleyCube(Scene):
    def construct(self):
        # 1. the twin models, real weights
        title = Text(f"{NUM['gene']}: one gene, two models", font_size=34, color=INK).to_edge(UP)
        self.play(Write(title), run_time=1.2)
        axis = Line(LEFT * 5.5, RIGHT * 5.5, color=MUTED, stroke_width=2)
        wE, wY = np.array(NUM["weights_eur"]), np.array(NUM["weights_yri"])
        scale = max(np.abs(wE).max(), np.abs(wY).max())
        xs = np.linspace(-5.3, 5.3, len(wE))
        stemsE = VGroup(*[Line([x, 0, 0], [x, 1.6 * np.sqrt(abs(v) / scale), 0], color=INDIGO, stroke_width=3)
                          for x, v in zip(xs, wE) if v != 0])
        stemsY = VGroup(*[Line([x, 0, 0], [x, -1.6 * np.sqrt(abs(v) / scale), 0], color=CAMWOOD, stroke_width=3)
                          for x, v in zip(xs, wY) if v != 0])
        labE = Text(f"European model, {NUM['n_nonzero_eur']} SNPs", font_size=22, color=INK).move_to(UP * 2.1 + LEFT * 3)
        labY = Text(f"Yoruba model, {NUM['n_nonzero_yri']} SNPs", font_size=22, color=INK).move_to(DOWN * 2.1 + LEFT * 3)
        self.play(Create(axis), run_time=0.6)
        self.play(LaggedStart(*[Create(s) for s in stemsE], lag_ratio=0.02), FadeIn(labE), run_time=1.6)
        self.play(LaggedStart(*[Create(s) for s in stemsY], lag_ratio=0.02), FadeIn(labY), run_time=1.6)
        self.wait(1.2)

        # 2. the variance formula and its three ingredients
        formula = MathTex(r"V", r"=", r"w^\top", r"D^{1/2}", r"R", r"D^{1/2}", r"w", color=INK, font_size=56)
        formula[2].set_color(PALM); formula[6].set_color(PALM)
        formula[3].set_color(GOLD); formula[5].set_color(GOLD)
        formula[4].set_color(CORAL)
        self.play(FadeOut(VGroup(axis, stemsE, stemsY, labE, labY)), run_time=0.8)
        self.play(Write(formula), run_time=1.5)
        keys = VGroup(
            Text("w  weights", font_size=24, color=PALM),
            Text("D  allele frequencies, 2p(1-p)", font_size=24, color=GOLD),
            Text("R  linkage disequilibrium", font_size=24, color=CORAL),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.25).next_to(formula, DOWN, buff=0.6)
        self.play(FadeIn(keys), run_time=1.0)
        note = Text("each ingredient can come from either population", font_size=24, color=INK).next_to(keys, DOWN, buff=0.5)
        self.play(FadeIn(note), run_time=0.8)
        self.wait(1.5)
        self.play(FadeOut(VGroup(formula, keys, note)), run_time=0.8)

        # 3. the cube of 8 mixes, drawn in an oblique projection
        def pos(bits):
            b = [int(c) for c in bits]
            return np.array([-2.2 + 3.2 * b[0] + 1.4 * b[2], -2.3 + 2.9 * b[1] + 1.0 * b[2], 0])

        corners = sorted(NUM["vertex_logV"])
        edges = VGroup()
        for a in corners:
            for k in range(3):
                if a[k] == "0":
                    b = a[:k] + "1" + a[k + 1:]
                    edges.add(Line(pos(a), pos(b), color=PLAYER_COLOR["wDR"[k]], stroke_width=3, stroke_opacity=0.45))
        dots = VGroup(*[Dot(pos(c), radius=0.08, color=INK) for c in corners])
        center = np.mean([pos(c) for c in corners], axis=0)

        def corner_label(c):
            direction = pos(c) - center
            direction = direction / np.linalg.norm(direction)
            lab = MathTex(f"{NUM['vertex_logV'][c]:.2f}", font_size=26, color=INK)
            lab.move_to(pos(c) + direction * 0.5)
            lab.add_background_rectangle(color="#ffffff", opacity=1, buff=0.05)
            return lab

        values = VGroup(*[corner_label(c) for c in corners])
        values.set_z_index(5)  # walk segments pass underneath the corner labels
        # Endpoint labels hang off their corner's value label, so the two can never overlap.
        value_of = dict(zip(corners, values))
        start = Text("all European", font_size=20, color=INDIGO).next_to(value_of["000"], DOWN, buff=0.08)
        end = Text("all Yoruba", font_size=20, color=CAMWOOD).next_to(value_of["111"], UP, buff=0.08)
        cube_title = Text("log V at every mix of weights, frequencies and LD", font_size=26, color=INK).to_edge(UP).shift(DOWN * 0.6)
        self.play(Transform(title, Text(f"{NUM['gene']}: the Shapley cube", font_size=34, color=INK).to_edge(UP)), run_time=0.6)
        self.play(FadeIn(cube_title), Create(edges), FadeIn(dots), run_time=1.6)
        edge_key = VGroup(*[
            VGroup(Line(ORIGIN, RIGHT * 0.5, color=PLAYER_COLOR[p], stroke_width=6),
                   Text(label, font_size=22, color=INK)).arrange(RIGHT, buff=0.2)
            for p, label in (("w", "swap weights"), ("D", "swap allele frequencies"), ("R", "swap LD"))
        ]).arrange(DOWN, aligned_edge=LEFT, buff=0.25).to_edge(RIGHT, buff=0.5).shift(UP * 0.8)
        self.play(FadeIn(values), FadeIn(start), FadeIn(end), FadeIn(edge_key), run_time=1.2)
        gap = MathTex(r"\Delta = " + signed(NUM["delta"]), font_size=40, color=INK).to_corner(RIGHT + DOWN).shift(UP * 0.4)
        self.play(Write(gap), run_time=0.8)
        self.wait(1.0)

        # 4. walk the 6 paths, one ingredient at a time
        tally = {"w": [], "D": [], "R": []}
        for path in NUM["paths"]:
            state = "000"
            walk = VGroup()
            for step in path["steps"]:
                k = "wDR".index(step["player"])
                nxt = state[:k] + "1" + state[k + 1:]
                seg = Line(pos(state), pos(nxt), color=PLAYER_COLOR[step["player"]], stroke_width=9)
                walk.add(seg)
                tally[step["player"]].append(step["change"])
                state = nxt
            order = MathTex(r" \rightarrow ".join(path["order"]), font_size=40, color=INK).to_corner(LEFT + DOWN).shift(UP * 0.4)
            self.play(Create(walk), FadeIn(order), run_time=1.3)
            self.wait(0.3)
            self.play(FadeOut(walk), FadeOut(order), run_time=0.4)

        # 5. average each ingredient's steps: the Shapley components
        self.play(FadeOut(VGroup(edges, dots, values, start, end, cube_title, edge_key)), run_time=0.8)
        rows = VGroup()
        for p in ("w", "D", "R"):
            steps = ", ".join(signed(v) for v in tally[p])
            row = MathTex(PLAYER_TEX[p] + r" = \tfrac{1}{6}\,(" + steps + r") = " + signed(NUM["phi"][p]),
                          font_size=30, color=PLAYER_COLOR[p])
            rows.add(row)
        rows.arrange(DOWN, aligned_edge=LEFT, buff=0.35).move_to(UP * 0.6)
        self.play(LaggedStart(*[Write(r) for r in rows], lag_ratio=0.5), run_time=3.0)
        total = MathTex(r"\phi_w + \phi_D + \phi_R = " + signed(sum(NUM["phi"].values())) + r" = \Delta",
                        font_size=40, color=INK).next_to(rows, DOWN, buff=0.6)
        self.play(Write(total), Indicate(gap, color=INK), run_time=1.6)

        # bar of the three components, widths proportional to |phi|
        mags = {p: abs(NUM["phi"][p]) for p in ("w", "D", "R")}
        width_total = 8.0
        bars = VGroup()
        for p in ("w", "D", "R"):
            bars.add(Rectangle(width=max(width_total * mags[p] / sum(mags.values()), 0.05), height=0.35,
                               fill_color=PLAYER_COLOR[p], fill_opacity=1, stroke_width=0))
        bars.arrange(RIGHT, buff=0.04).next_to(total, DOWN, buff=0.5)
        caption = VGroup(
            Text("most of this gene's gap comes from allele frequency", font_size=24, color=INK),
            Text("bar widths show size; the weights part points the other way", font_size=20, color=MUTED),
        ).arrange(DOWN, buff=0.15).next_to(bars, DOWN, buff=0.3)
        self.play(Create(bars), FadeIn(caption), run_time=1.4)
        self.wait(3.0)
