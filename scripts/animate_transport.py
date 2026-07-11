"""
Manim animation of the 1D advection-diffusion result produced by the
Fortran transport solver (app/main_transport.f90).

Run the Fortran solver first so results/transport1d_t*.dat exist:
    make run-transport

Then render with Manim:
    manim -qm -pql scripts/animate_transport.py TransportWave
    (or, inside a Jupyter cell:  %%manim -qm TransportWave)

This reads the actual numerical snapshots -- it is not a symbolic /
idealized redraw of the PDE, it's the real Fortran output.
"""
import glob
import numpy as np
from manim import *


class TransportWave(Scene):
    def construct(self):
        files = sorted(glob.glob("results/transport1d_t*.dat"))
        snapshots = [np.loadtxt(f) for f in files]

        x = snapshots[0][:, 0]
        u_all = np.array([s[:, 1] for s in snapshots])
        umax = u_all.max() * 1.15

        axes = Axes(
            x_range=[x.min(), x.max(), 0.2],
            y_range=[0, umax, umax / 5],
            x_length=10,
            y_length=5,
            axis_config={"include_tip": True},
        )
        labels = axes.get_axis_labels(x_label="x", y_label="u(x,t)")

        title = Text("1D advection-diffusion: du/dt + v du/dx = D d\u00b2u/dx\u00b2",
                      font_size=28).to_edge(UP)

        def make_curve(u):
            pts = list(zip(x, u))
            return axes.plot_line_graph(
                x_values=x, y_values=u, line_color=BLUE, add_vertex_dots=False
            )["line_graph"]

        curve = make_curve(u_all[0])
        self.play(Write(title), Create(axes), Write(labels))
        self.play(Create(curve))
        self.wait(0.3)

        for u in u_all[1:]:
            new_curve = make_curve(u)
            self.play(Transform(curve, new_curve), run_time=0.25)

        self.wait(1)
