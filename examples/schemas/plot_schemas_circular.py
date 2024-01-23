"""
=======================================
Understanding circular HiC contact maps
=======================================

"""

# This is not super clean. I should define the proper transformation, and be
# able to map directly from genomic coordinates to circular coordinates.

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import circhic


lengths = np.array([50])

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(polar=True)

ax.set_ylim(0, 50)
ax.set_rorigin(-20)

resolution = 1000
theta = np.linspace(0, 2*np.pi, resolution)
ones = np.ones(theta.shape)
ax.plot(theta, ones*0, linewidth=4, color="black", zorder=11)
ax.plot(theta, ones*21, linewidth=3, color="r", zorder=11)
ax.plot(theta, ones*(-21), linewidth=3, color="r", zorder=11)
ax.set_rticks([])
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.tick_params(which='major', width=1.0, length=15)
ax.tick_params(which='major', length=10)
ax.tick_params(which='minor', width=1.0, labelsize=10)
ax.tick_params(which='minor', length=5, labelsize=10, labelcolor='0.25')


# The original rectangle is at (34, 25)
theta, r = circhic.utils.convert_xy_to_thetar(
    (np.array([33.5, 34.5, 24.5, 25.5]), np.array([24.5, 25.5, 33.5, 34.5])),
    lengths=lengths)

rect = Rectangle((theta[0], r[0]), theta[1]-theta[0], 1,
                 facecolor="black", fill=True)
ax.add_artist(rect)

rect = Rectangle((theta[2], r[2]), theta[3]-theta[2], 1,
                 facecolor="black", fill=True, angle=0)
rect.get_path()._interpolation_steps = 100

ax.add_artist(rect)


# Now plot the two dotted line to indicate loci i and j
# theta, r = circhic.utils.convert_xy_to_thetar(
#           (24.5 * np.ones(resolution), np.linspace(-0, 100, resolution)),
# lengths)
#ax.plot(theta, r, linestyle=":", color="0.25", linewidth=0.5)
#theta, r = circhic.utils.convert_xy_to_thetar(
#          (25.5 * np.ones(resolution), np.linspace(-0, 100, resolution)),
#          lengths)
#ax.plot(theta, r, linestyle=":", color="0.25", linewidth=0.5)
#
#
#theta, r = circhic.utils.convert_xy_to_thetar(
#          (np.linspace(-0, 100, resolution), 33.5 * np.ones(resolution)),
#          lengths)
#ax.plot(theta, r, linestyle=":", color="0.25", linewidth=0.5)
#theta, r = circhic.utils.convert_xy_to_thetar(
#          (np.linspace(-0, 100, resolution), 34.5 * np.ones(resolution)),
#          lengths)
#ax.plot(theta, r, linestyle=":", color="0.25", linewidth=0.5)


theta, r = circhic.utils.convert_xy_to_thetar(
          (np.linspace(25, 100, resolution), 34 * np.ones(resolution)),
          lengths)
ax.plot(theta, r, linestyle=":", color="0.25", linewidth=0.5)

theta, r = circhic.utils.convert_xy_to_thetar(
          (25*np.ones(resolution), np.linspace(-0, 35, resolution)),
          lengths)
ax.plot(theta, r, linestyle=":", color="0.25", linewidth=0.5)

genomic_ticks = np.arange(0, 50, 10)
ax.set_thetagrids(genomic_ticks / 50 * 360, [0, 10, 20, 30, 40])
# Create grid manually
for i in np.arange(0, 50, 10):
    theta, r = circhic.utils.convert_xy_to_thetar(
          (i * np.ones(resolution), np.linspace(-100, 100, resolution)),
          lengths)
    ax.plot(theta, r, linestyle="--", color="0.25", linewidth=0.5)

    theta, r = circhic.utils.convert_xy_to_thetar(
          (np.linspace(-100, 100, resolution), i * np.ones(resolution)),
          lengths)
    ax.plot(theta, r, linestyle="--", color="0.25", linewidth=0.5)


def circle(ax, x, y, radius=2):
    from matplotlib.patches import Circle
    from matplotlib.patheffects import withStroke
    circle = Circle((x, y), radius, clip_on=False, zorder=12, linewidth=1,
                    edgecolor='black', facecolor=(0, 0, 0, .0125),
                    path_effects=[withStroke(linewidth=5, foreground='w')])
    ax.add_artist(circle)


def text(ax, x, y, text, zorder=10):
    ax.text(x, y, text, backgroundcolor="white",
            ha='center', va='top', weight='bold', color='blue', zorder=zorder)

ax.grid(linestyle="--", linewidth=0, color=".25")

ax.set_title("Circular contact map",
             fontsize=20, verticalalignment='bottom')

main_ax = fig.add_subplot(facecolor="none")
circle(main_ax, 0.71, 0.7, radius=0.04)
text(main_ax, 0.71, 0.63, "k-th diagonal")

circle(main_ax, 0.64, 0.46, radius=0.04)
text(main_ax, 0.64, 0.39, "Main diagonal")

circle(main_ax, 0.395, 0.316, radius=0.04)
text(main_ax, 0.395, 0.246, "Contact count\nbetween i and j")


circle(main_ax, 0.495, 0.35, radius=0.04)
text(main_ax, 0.495, 0.28, "Loci j")

circle(main_ax, 0.37, 0.43, radius=0.04)
text(main_ax, 0.37, 0.38, "Loci i")


main_ax.spines["left"].set_visible(False)
main_ax.spines["top"].set_visible(False)
main_ax.spines["bottom"].set_visible(False)
main_ax.spines["right"].set_visible(False)
main_ax.set_xticks([])
main_ax.set_yticks([])
