import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.gridspec import GridSpec
from .tools import genCircData as _generate_circular_data


class CircHiCFigure:
    name = "circhic"

    def __init__(self, lengths, origin=1, resolution=None, figure=None):
        """
        A circular HiC figure

        Parameters
        ----------
        lengths : ndarray

        origin : integer, optional, default: 1
            position of the origin. The origin is set to the east of the plot

        figure : matplotlib.figure.Figure, optional, default: None
            A Matplotlib figure. If not provided, will create it.

        """
        # If figure is not provided, create a square figure.
        self.fig = figure if figure is not None else plt.figure(figsize=(8, 8))

        # Create a gridspec : 1000 x 1000 should be enough for a high
        # resolution placements of axes.
        self._gridspec = GridSpec(nrows=1000, ncols=1000, figure=self.fig)

        self.lengths = lengths
        self.origin = origin
        self.resolution = resolution if resolution is not None else 1
        self._polar_axes = []

    def plot_hic(self, counts, inner_gdis=None, outer_gdis=None,
                 inner_radius=0, outer_radius=1):
        ax = self._create_subplot(
            outer_radius, polar=False, zorder=-99,
            label=("hic_%d" % (len(self._polar_axes)+1)))

        if outer_gdis is None:
            outer_gdis = int(np.round(counts.shape[0] / 2 * self.resolution))
        if inner_gdis is None:
            inner_gdis = int(np.round(counts.shape[0] / 2 * self.resolution))

        # Generate circular hic map
        circular_data = _generate_circular_data(
            counts, res=self.resolution,
            pos0=self.origin, r_in=inner_radius, s_in=inner_gdis,
            s_out=outer_gdis)
        im = ax.imshow(circular_data, norm=colors.SymLogNorm(1), zorder=-99)
        ax.set_axis_off()
        return (im, ax)

    def plot_marks(self, marks, s_out=1, s_in=1, r_in=0, outer_radius=1,
                   inner_radius=0, zorder=None):
        ax_m = self._create_subplot(outer_radius)

        for name, mark in marks.items():
            if 'color' not in mark:
                color = 'white'
            else:
                color = mark['color']

            if 'marker' not in mark:
                marker = 'o'
            else:
                marker = mark['marker']

            if 'ms' not in mark:
                ms = 14
            else:
                ms = mark['ms']

            theta = [np.pi/2-mark['bin']*2*np.pi / self.lengths.sum()]
            r = [r_in + s_in*(1 - r_in)/(s_out + s_in)]
            ax_m.plot(theta, r, marker, ms=ms, color=color, zorder=zorder)
            ax_m.set_rmax(1)

        ax_m.set_axis_off()
        self._polar_axes += [ax_m]

    def plot_lines(self, data, color=None, linestyle=None,
                   inner_radius=0, outer_radius=1, zorder=None):

        ax_g = self._create_subplot(
            outer_radius=outer_radius,
            label=("lines_%d" % (len(self._polar_axes)+1)))

        # Need to include the theta shift here.
        theta = np.array(
            [np.pi/2-i*2*np.pi/len(data) for i in range(len(data))])

        lines = ax_g.plot(
            np.concatenate((theta, [theta[0]])),
            np.concatenate((data, [data[0]])),
            color=color, linestyle=linestyle,
            zorder=zorder)

        # Now compute the new origin
        rorigin = (
            (np.nanmin(data) - np.nanmax(data)) * outer_radius /
            (outer_radius - inner_radius))
        ax_g.set_rmin(rorigin)

        ax_g.set_axis_off()
        self._polar_axes += [ax_g]
        return (lines, ax_g)

    def plot_bars(self, data, color=None, inner_radius=0, outer_radius=1,
                  zorder=None):
        ax = self._create_subplot(
            outer_radius=outer_radius,
            label=("bars_%d" % (len(self._polar_axes)+1)))
        theta = np.array(
            [np.pi/2-i*2*np.pi/len(data) for i in range(len(data))])
        width = theta[1] - theta[0]
        bars = ax.bar(theta, data, color=color, width=width, zorder=zorder)

        # Now compute the new origin
        rorigin = (
            (np.nanmin(data) - np.nanmax(data)) * outer_radius /
            (outer_radius - inner_radius))
        ax.set_rmin(rorigin)

        ax.set_axis_off()
        self._polar_axes += [ax]
        return (bars, ax)

    def set_genomic_ticklabels(self, ticklabels=None, tickpositions=None):
        ax = self._create_subplot(label="thetaticks")
        ax.set_rgrids([])
        if tickpositions is not None:
            tickpositions = (
                tickpositions / (self.lengths.sum() *
                                 self.resolution) *
                2 * np.pi)
            ax.set_xticks(tickpositions)

        if ticklabels is None:
            theta_ticks = (ax.get_xticks() / (2*np.pi) * self.lengths.sum())
            ticklabels = [
                "%d" % np.round(s)
                for s
                in theta_ticks]
        ax.set_xticklabels(ticklabels, fontsize="x-small")
        ax.spines["polar"].set_linewidth(0)
        ax.spines["inner"].set_linewidth(0)
        ax.xaxis.grid(False)
        return ax

    def _create_subplot(self, outer_radius=1, polar=True, label=None,
                        zorder=None):
        if outer_radius == 1:
            ax_g = self.fig.add_subplot(
                polar=polar,
                facecolor="none", label=label,
                zorder=zorder)
        else:
            nrows = int(np.round((1 - outer_radius) / 2 * 1000))
            ax_g = self.fig.add_subplot(
                self._gridspec[nrows:-nrows, nrows:-nrows],
                facecolor="none",
                polar=polar,
                label=label, zorder=zorder)
        if polar:
            theta_offset = (
                (self.origin - 1) / (self.lengths.sum() * self.resolution))
            ax_g.set_theta_zero_location("E", offset=theta_offset)
        return ax_g
