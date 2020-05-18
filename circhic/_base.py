import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.gridspec import GridSpec
from .utils import generate_circular_map as _generate_circular_data


class CircHiCFigure:
    """
    A circular HiC figure

    Parameters
    ----------
    lengths : ndarray
        array of chromosome length.

    origin : integer, optional, default: 1
        position of the origin. The origin is set to the east of the plot

    resolution : integer, optional, default: None
        The resolution of the contact count matrix.

    figure : matplotlib.figure.Figure, optional, default: None
        A Matplotlib figure. If not provided, will create it.

    Notes
    -----
    See FIXME
    """

    name = "circhic"

    def __init__(self, lengths, origin=1, resolution=None, figure=None):
        # If figure is not provided, create a square figure.
        self.figure = (
            figure if figure is not None else plt.figure(figsize=(8, 8)))

        # Create a gridspec : 1000 x 1000 should be enough for a high
        # resolution placements of axes.
        self._gridspec = GridSpec(nrows=1100, ncols=1100, figure=self.figure)

        self.lengths = lengths
        self.origin = origin
        self.resolution = resolution if resolution is not None else 1
        self._polar_axes = []

    def plot_hic(self, counts, inner_gdis=None, outer_gdis=None,
                 inner_radius=0, outer_radius=1,
                 cmap="viridis",
                 ax=None):
        """
        Plot a heatmap of the HiC contact count matrix on a circular strip.

        Parameters
        ----------

        counts : ndarray (n, n)
            The contact count matrix of shape (n, n) where
            `n = lengths.sum() / resolution`

        inner_gdis : integer, optional, default: None
            Plot up to `inner_gdis` of the diagonal of the contact count
            matrix (in genomic distance). Corresponds to the lower-diagonal on
            a typical square HiC contact count matrix.

        outer_gdis : integer, optional, default: None
            Plot up to `outer_gdis` of the diagonal of the contact count matrix
            (in genomic distance). Corresponds to the upper-diagonal part of
            the contact count matrix on a typical square contact count map.

        inner_radius : float, optional, default: 0
            Radius of the inner strip, considering that the maximum outer
            radius is 1. Should be smaller than `outer_radius`.
            Note that `inner_radius` will be ignored if ax is provided.

        outer_radius : float, optional, default: 1
            Radius of the outer strip, considering that the maximum possible
            outer radius is 1. Should be larger than `inner_radius`.

        cmap : string, optional, default : "viridis"
            A Matplotlib colormap.

        ax : matplotlib.axes.Axes object, optional, default: None
            Matplotlib Axes object. By default, will create one. Note that
            outer_radius and inner_radius will be ignored if `ax` is provided.

        Returns
        -------
        (im, ax) type of artist and axes

        """
        if ax is None:
            ax = self._create_subplot(
                outer_radius, polar=False, zorder=-99,
                label=("hic_%d" % (len(self._polar_axes)+1)))
        else:
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines["left"].set_linewidth(0)
            ax.spines["top"].set_linewidth(0)
            ax.spines["bottom"].set_linewidth(0)
            ax.spines["right"].set_linewidth(0)

            rect = ax.get_position(original=True).bounds
            # Rect is left, bottom, width, height
            # This needs to be reduced by height*outer_radius and
            # width*outer_radius
            new_width = rect[2] * outer_radius
            new_height = rect[3] * outer_radius
            new_left = rect[0] + rect[2] * (1-outer_radius) / 2
            new_bottom = rect[1] + rect[3] * (1-outer_radius) / 2
            ax = self.figure.add_axes(
                (new_left, new_bottom, new_width, new_height),
                facecolor="none")

        if outer_gdis is None:
            outer_gdis = int(np.round(counts.shape[0] / 2 * self.resolution))
        if inner_gdis is None:
            inner_gdis = int(np.round(counts.shape[0] / 2 * self.resolution))

        # Need to convert inner_radius to what _generate_circular_data
        # expects (outer_radius = 1)
        cir_inner_radius = inner_radius / outer_radius

        # Generate circular hic map
        circular_data = _generate_circular_data(
            counts, res=self.resolution,
            origin=self.origin, inner_radius=cir_inner_radius,
            inner_gdis=inner_gdis,
            outer_gdis=outer_gdis)
        im = ax.imshow(
            circular_data, interpolation=None,
            norm=colors.SymLogNorm(1, base=10), cmap=cmap)

        # We don't want to remove entirely the axis, as it means setting
        # xlabels and ylabels don't work anymore.
        # Actually, is that an issue apart from the gallery?? I'm not so sure
        # anymore…
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["left"].set_linewidth(0)
        ax.spines["top"].set_linewidth(0)
        ax.spines["bottom"].set_linewidth(0)
        ax.spines["right"].set_linewidth(0)
        return (im, ax)

    def _plot_marks(self, marks, s_out=1, s_in=1, r_in=0, outer_radius=1,
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
        """
        Plot a line chart

        Parameters
        ----------
        data : ndarray (n, )

        color : a compatible matplotlib color, optional, default: None
            The line color. Possible values:

            - A single color format string (e.g. "#000000", "black", "0").
            - A float between 0 and 1

            Defaults to `None`.

        linestyle : a compatible Matplotlib linestyle, optional, default: None

        inner_radius : float (0, 1), optional, default: 0
            The inner radius of the plot, assuming the maximum outer radius
            possible is 1. Should be smaller than `outer_radius`.

        outer_radius : float (0, 1), optional, default: 1
            The outer radius of the plot, assuming the maximum outer radius
            possible is 1. Should be larger than `inner_radius`.

        zorder : float

        Returns
        -------
        (lines, ax)
        """
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
        """
        Plot a bar chart

        Parameters
        ----------
        data : ndarray (n, )

        color : a compatible matplotlib color, optional, default: None
            The line color. Possible values:

            - A single color format string (e.g. "#000000", "black", "0").
            - A float between 0 and 1

            Defaults to `None`.

        linestyle : a compatible Matplotlib linestyle, optional, default: None

        inner_radius : float (0, 1), optional, default: 0
            The inner radius of the plot, assuming the maximum outer radius
            possible is 1. Should be smaller than `outer_radius`.

        outer_radius : float (0, 1), optional, default: 1
            The outer radius of the plot, assuming the maximum outer radius
            possible is 1. Should be larger than `inner_radius`.

        zorder : float

        Returns
        -------
        (artists, ax)
        """
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

    def set_genomic_ticklabels(self, ticklabels=None, tickpositions=None,
                               ax=None):
        """
        Set the circular tick labels

        Parameters
        ----------

        ticklabels : array-like of strings
            the list of strings to plot. Should be the same length as the
            number of ticks.

        tickpositions : array of floas
            the positions of the ticks. Should be the same length as the tick
            labels.

        ax : matplotlib.axes.Axes object, optional, default: None
            Matplotlib Axes object. By default, will create one. Note that
            outer_radius and inner_radius will be ignored if `ax` is provided.

        """
        if ax is None:
            ax = self._create_subplot(label="thetaticks")
        else:
            rect = ax.get_position().bounds
            ax = self.figure.add_axes(rect, polar=True, facecolor="none")

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

    def set_colorbar(self, mappable, orientation="vertical"):
        """
        Set a colorbar on the plot

        Parameters
        ----------
        mappable : matplotlib.cm.ScalarMappable
            The matplotlib.cm.ScalarMappable (i.e., Image, ContourSet, etc.)
            described by this colorbar.

        orientation : {"vertical", "horizontal"}, default: "vertical"
            Whether to plot a vertical or horizontal colorbar.
        """
        if orientation == "vertical":
            ax = self.figure.add_subplot(
                self._gridspec[:1000, 1070:1100])
        else:
            ax = self.figure.add_subplot(
                self._gridspec[1070:1100, :1000])

        cab = self.figure.colorbar(mappable, cax=ax, orientation=orientation)
        return cab

    def _create_subplot(self, outer_radius=1, polar=True, label=None,
                        zorder=None):
        nrows = int(np.round((1 - outer_radius) / 2 * 1000))
        ax_g = self.figure.add_subplot(
            self._gridspec[nrows:-nrows-100, nrows:-nrows-100],
            facecolor="none",
            polar=polar,
            label=label, zorder=zorder)

        if polar:
            theta_offset = (
                (self.origin - 1) / (self.lengths.sum() * self.resolution) *
                360)
            ax_g.set_theta_zero_location("E", offset=theta_offset)
        return ax_g
