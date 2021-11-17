import numpy as np
import warnings
import matplotlib  # Need this to check the version
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.gridspec import GridSpec
from matplotlib import patches
from matplotlib.container import BarContainer
from matplotlib import rcParams
from matplotlib import ticker
from .utils import generate_circular_map
from .utils import generate_borders


class CircHiCFigure:
    """
    A circular HiC figure

    Parameters
    ----------
    lengths : ndarray
        array of chromosome length.

    origin : integer, optional, default: 1
        position of the origin. The origin is set to the North of the plot

    chromosome_type : {"circular", "linear"}, optional, default: "circular"
        whether to plot a circular or a linear chromosome.

    figure : matplotlib.figure.Figure, optional, default: None
        A Matplotlib figure. If not provided, will create it.

    Notes
    -----
    See FIXME
    """

    name = "circhic"

    def __init__(self, lengths, origin=1, chromosome_type="circular",
                 figure=None):
        # If figure is not provided, create a square figure.
        self.figure = (
            figure if figure is not None else plt.figure(figsize=(8, 8)))

        # Create a gridspec : 1000 x 1000 should be enough for a high
        # resolution placements of axes.
        self._gridspec = GridSpec(nrows=1150, ncols=1150, figure=self.figure)

        self.lengths = lengths
        self.origin = origin
        self._polar_axes = []
        if chromosome_type not in ["circular", "linear"]:
            raise ValueError(
                "Mode %s unknown. Expected 'circular' or 'linear'" %
                chromosome_type)
        self.chromosome_type = chromosome_type

    def plot_hic(self, counts, inner_gdis=None, outer_gdis=None,
                 inner_radius=0, outer_radius=1,
                 resolution=1,
                 cmap="viridis",
                 mode="reflect",
                 vmin=None,
                 vmax=None,
                 alpha=1,
                 border_thickness=0.02,
                 ax=None,
                 granularity=0.5):
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

        resolution : integer, optional, default: None
            Resolution of the HiC contact count map. By default, the function
            will estimate the resolution given the lengths of the chromosome
            and the shape of the contact count matrix.

        cmap : string, optional, default : "viridis"
            A Matplotlib colormap.

        mode : {"reflect", "distant"}, optional, default: "reflect"

            - if `"reflect"`, the contact count matrix will be plotted from
              `inner_gdis` to 0, then 0 to `outer_gdis`.
            - if `"distant"`, the contact count matrix will be plotted from
              `inner_gdis` to `outer_gdis`: this option is useful to visualize
              contact counts far away from the diagonal.

        vmin, vmax : float, optional, default: None
            `vmin` and `vmax` define the data range that the colormap covers.
            By default, the colormap covers the complete value range of the
            supplied data.


        ax : matplotlib.axes.Axes object, optional, default: None
            Matplotlib Axes object. By default, will create one. Note that
            outer_radius and inner_radius will be ignored if `ax` is provided.

        Returns
        -------
        (im, ax) type of artist and axes
        """

        if mode not in ["reflect", "distant"]:
            raise ValueError(
                "mode %s is unknown. Possible values for mode are "
                "'reflect', and 'distant'.")
        n = counts.shape[0]
        if resolution is None:
            resolution = self.lengths.sum() / n
        else:
            estimate_shape = self.lengths.sum() / resolution
            if not ((n - 1) < estimate_shape < (n + 1)):
                warnings.warn(
                    "The resolution provided does not match the shape of the "
                    "contact count matrix and the length of the chromosome")

        if ax is None:
            ax = self._create_subplot(
                outer_radius, polar=False, zorder=-99,
                resolution=resolution,
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

        resolution = resolution if resolution is not None else 1

        if outer_gdis is None:
            outer_gdis = int(np.round(counts.shape[0] / 2 * resolution))
        if inner_gdis is None:
            inner_gdis = int(np.round(counts.shape[0] / 2 * resolution))

        # Need to convert inner_radius to what _generate_circular_data
        # expects (outer_radius = 1)
        cir_inner_radius = inner_radius / outer_radius

        # Generate circular hic map
        circular_data = generate_circular_map(
            counts, resolution=resolution,
            granularity=granularity,
            origin=self.origin, inner_radius=cir_inner_radius,
            mode=mode,
            chromosome_type=self.chromosome_type,
            inner_gdis=inner_gdis,
            outer_gdis=outer_gdis)
        if vmin is None:
            if matplotlib.__version__ < "3.2.0":
                norm = colors.SymLogNorm(1, vmax=vmax)
            else:
                norm = colors.SymLogNorm(1, base=np.e, vmax=vmax)
        else:
            if matplotlib.__version__ < "3.2.0":
                norm = colors.SymLogNorm(vmin, vmax=vmax)
            else:
                norm = colors.SymLogNorm(vmin, base=np.e, vmax=vmax)

        if outer_gdis != 0:
            extent = (-outer_gdis, outer_gdis, -outer_gdis, outer_gdis)
        else:
            # I'm not sure this makes sense?
            extent = (-inner_gdis, inner_gdis, -inner_gdis, inner_gdis)

        im = ax.imshow(
            circular_data,
            interpolation="none",
            alpha=alpha,
            norm=norm, cmap=cmap,
            extent=extent,
            )
        if border_thickness != 0:
            border_im = generate_borders(
                counts, granularity=granularity,
                chromosome_type=self.chromosome_type,
                mode=mode,
                resolution=resolution,
                origin=self.origin, inner_radius=cir_inner_radius,
                inner_gdis=inner_gdis,
                outer_gdis=outer_gdis,
                thick_r=border_thickness,
                thick_extreme=border_thickness)
            ax.imshow(
                border_im, interpolation="none",
                extent=extent,
                cmap="Greys_r")
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

        # Now, keep stuff in memory for raxis
        contact_count_maps = {
            "outer_gdis": outer_gdis,
            "inner_gdis": inner_gdis,
            "outer_radius": outer_radius,
            "inner_radius": inner_radius,
            "mode": mode,
            "resolution": resolution
        }
        self._contact_count_maps = contact_count_maps
        return (im, ax)

    def plot_raxis(self):
        """Plot the r-axis, corresponding to the genomic distance
        """
        if not hasattr(self, "_contact_count_maps"):
            raise ValueError(
                "The raxis can only be plotted if a contact map is plotted")
        sspines = self._plot_raxis(
            self._contact_count_maps["outer_radius"],
            self._contact_count_maps["inner_radius"],
            self._contact_count_maps["outer_gdis"],
            self._contact_count_maps["inner_gdis"],
            self._contact_count_maps["resolution"],
            mode=self._contact_count_maps["mode"])
        return sspines

    def _plot_raxis(self, outer_radius, inner_radius, outer_gdis, inner_gdis,
                    resolution=None, mode="reflect"):
        if resolution is None:
            resolution = 1

        nrows = int(np.round((1 - outer_radius) / 2 * 1000))
        side_spine = self.figure.add_subplot(
            self._gridspec[50+nrows:550, 10:550],
            facecolor="none")
        side_spine.spines["right"].set_linewidth(0)
        side_spine.spines["top"].set_linewidth(0)
        side_spine.spines["bottom"].set_linewidth(0)
        side_spine.set_xticks([])

        if mode == "reflect":
            rorigin = (
                -inner_gdis - inner_radius *
                (outer_gdis+inner_gdis) / (outer_radius-inner_radius))
            yticks = [-inner_gdis, 0, outer_gdis]
            ticklabels = ["%d" % (inner_gdis * resolution),
                          "0",
                          "%d" % (outer_gdis * resolution)]

        else:
            rorigin = (
                inner_gdis - inner_radius *
                (outer_gdis - inner_gdis) / (outer_radius - inner_radius))
            yticks = [inner_gdis, outer_gdis]
            ticklabels = ["%d" % (inner_gdis * resolution),
                          "%d" % (outer_gdis * resolution)]

        y_bottom_lim = rorigin

        side_spine.set_ylim((y_bottom_lim, outer_gdis))
        if mode == "reflect":
            # Dude !! Use the grid option
            side_spine.axhline(0, linestyle="--",
                               linewidth=rcParams['axes.linewidth'],
                               color=rcParams['axes.edgecolor'])

            side_spine.axhline(
                -inner_gdis, linestyle="--",
                linewidth=rcParams["axes.linewidth"],
                color=rcParams['axes.edgecolor'])
            side_spine.spines["left"].set_bounds(-inner_gdis, outer_gdis)
        else:
            side_spine.axhline(
                inner_gdis, linestyle="--",
                linewidth=rcParams["axes.linewidth"],
                color=rcParams["axes.edgecolor"])
            side_spine.spines["left"].set_bounds(inner_gdis, outer_gdis)

        side_spine.axhline(
            outer_gdis, linestyle="--", linewidth=rcParams["axes.linewidth"],
            color=rcParams["axes.edgecolor"])
        side_spine.tick_params(
            colors=rcParams["axes.edgecolor"])
        side_spine.set_yticks(yticks)

        # Need to set the label coordinate such that it is in the middle, like
        # in a normal Matplotlib axis

        # I would much rather have the same default as matplotlib, ie a
        # vertical label, at the midpoint. In order to do this, we either need
        # to reduce the size of the yaxis to match exactly what we need or we
        # need to set the coordinate of the label, but taking in account the
        # size of the ticklabels.

        # This isn't used, but let's keep for now.
        # y_mid_point = 1 - (
        #     (outer_gdis - rorigin) / (outer_gdis - y_bottom_lim) * 0.5)

        side_spine.set_yticklabels(ticklabels)
        side_spine.yaxis.label.set_verticalalignment("bottom")
        side_spine.yaxis.label.set_rotation("horizontal")
        side_spine.yaxis.set_label_coords(-0.1, 1.02)

        return side_spine

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
        if self.chromosome_type == "circular":
            theta = np.array(
                [i*2*np.pi/len(data) for i in range(len(data))])
            lines = ax_g.plot(
                np.concatenate((theta, [theta[0]])),
                np.concatenate((data, [data[0]])),
                color=color, linestyle=linestyle,
                zorder=zorder)
        else:
            theta = np.array(
                [(len(data)/2-i)*(2*np.pi*0.7)/len(data)
                 for i in range(len(data))])
            lines = ax_g.plot(
                theta, data,
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
        if self.chromosome_type == "circular":
            theta = np.array(
                [i*2*np.pi/len(data) for i in range(len(data))])
        else:
            theta = np.array(
                [(len(data)/2-i)*2*np.pi*0.7/(len(data))
                 for i in range(len(data))])

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

    def plot_bands(self, begin, end, colors=None,
                   inner_radius=0, outer_radius=1):
        """
        Plot bands

        Parameters
        ----------
        begin : ndarray (l, )
            List of the all the beginnings of each band.

        end : ndarray (l, )
            List of all the ends of each band.

        colors : ndarray (l, )
            List of colors. Should be the same size as `begin` and `end`

        Returns
        -------
        (artists, ax)
        """
        ax = self._create_subplot(
            outer_radius=outer_radius,
            label=("bands_%d" % (len(self._polar_axes)+1)))
        ax.set_axis_off()

        n_bins = self.lengths.sum()
        left = begin
        right = end
        if self.chromosome_type == "circular":
            # Convert the left hand side of the rectangle to the correct
            # angular form.
            left = np.array(
                [i*2*np.pi/n_bins for i in left])
            # Do the same with the end of the band
            right = np.array(
                [i*2*np.pi/n_bins for i in end])
        else:
            # here, we tackle the case of the linear chromosome
            left = np.array(
                [(n_bins/2 - i)*2*np.pi*0.7/n_bins for i in left])
            # Do the same with the end of the band
            right = np.array(
                [(n_bins/2 - i)*2*np.pi*0.7/n_bins for i in end])

        width = right - left
        height = 1
        bottom = 0

        _patches = []

        # Now compute the new origin
        rorigin = _compute_rorigin(0, 1, inner_radius, outer_radius)
        if matplotlib.__version__ == "3.2.0":
            ax.set_rmin(rorigin)
        else:
            ax.set_rorigin(rorigin)

        for i, (l, w) in enumerate(zip(left, width)):
            if colors is not None:
                c = colors[i]
            else:
                c = "C%d" % i
            r = patches.Rectangle(
                xy=(l, bottom), width=w, height=height,
                facecolor=c,
                linewidth=0,
                label='_nolegend_',
                )
            r.get_path()._interpolation_steps = 100
            r.sticky_edges.y.append(0)
            r.sticky_edges.y.append(1)

            ax.add_patch(r)
            _patches.append(r)

        bar_container = BarContainer(_patches)
        ax.add_container(bar_container)

        # Now compute the new origin
        rorigin = _compute_rorigin(0, 1, inner_radius, outer_radius)
        ax.set_rmin(rorigin)

        if matplotlib.__version__ < "3.2.0":
            ax.autoscale_view()
        else:
            ax._request_autoscale_view()

        return bar_container, ax

    def set_genomic_ticklabels(self, outer_radius=1, ticklabels=None,
                               tickpositions=None,
                               fontdict=None,
                               ax=None):
        """
        Set the circular tick labels

        Parameters
        ----------

        ticklabels : array-like of strings
            the list of strings to plot. Should be the same length as the
            number of ticks.

        tickpositions : array of floats
            the positions of the ticks. Should be the same length as the tick
            labels.

        fontdict : dict, optional
            A dictionary controlling the appearance of the ticklabels.
            The default `fontdict` is:

                {'fontsize': rcParams['axes.titlesize'],
                 'fontweight': rcParams['axes.titleweight'],
                 'verticalalignment': 'baseline',
                 'horizontalalignment': loc}

            See the Maplotlib documentation for more information on the
            fontdict.

        ax : matplotlib.axes.Axes object, optional, default: None
            Matplotlib Axes object. By default, will create one. Note that
            outer_radius and inner_radius will be ignored if `ax` is provided.

        """
        if ax is None:
            ax = self._create_subplot(label="thetaticks",
                                      outer_radius=outer_radius)
        else:
            rect = ax.get_position().bounds
            ax = self.figure.add_axes(rect, polar=True, facecolor="none")

        ax.set_rgrids([])
        if tickpositions is not None:
            if self.chromosome_type == "circular":
                tickpositions = (
                    tickpositions / self.lengths.sum() *
                    2 * np.pi)
            else:
                tickpositions = (
                    (tickpositions - self.lengths.sum()/2) /
                    self.lengths.sum() *
                    2 * np.pi * 0.7)

            ax.set_xticks(tickpositions)

        if ticklabels is None:
            if self.chromosome_type == "circular":
                theta_ticks = (
                    ax.get_xticks() / (2*np.pi) * self.lengths.sum())
            else:
                xticks = np.array(
                    [(self.lengths.sum()/2-i) * (2*np.pi*0.7) /
                     self.lengths.sum()
                     for i in np.linspace(0, self.lengths.sum(), 6)])
                ax.set_xticks(xticks)
                ax.set_thetalim(-np.pi, np.pi)
                theta_ticks = np.linspace(0, self.lengths.sum(), 6)
            ticklabels = [
                "%d" % np.round(s)
                for s
                in theta_ticks]
        ax.set_xticklabels(ticklabels, fontdict=fontdict)
        ax.spines["polar"].set_linewidth(0)
        ax.spines["inner"].set_linewidth(0)
        ax.xaxis.grid(False)
        return ax

    def set_colorbar(self, mappable, orientation="vertical",
                     tick_formatter=None):
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
                self._gridspec[50:1050, 1130:1150])
        else:
            ax = self.figure.add_subplot(
                self._gridspec[1130:1150, 50:1050])
        ax.tick_params(axis='both', which='major', labelsize="x-small")
        if tick_formatter is None:
            tick_formatter = ticker.ScalarFormatter()
        cab = self.figure.colorbar(mappable, cax=ax, orientation=orientation,
                                   format=tick_formatter)
        return cab

    def _create_subplot(self, outer_radius=1, polar=True, label=None,
                        zorder=None, resolution=None):
        resolution = resolution if resolution is not None else 1
        nrows = int(np.round((1 - outer_radius) / 2 * 1000))
        ax_g = self.figure.add_subplot(
            self._gridspec[50+nrows:-nrows-100, 50+nrows:-nrows-100],
            facecolor="none",
            polar=polar,
            label=label, zorder=zorder)

        if polar:
            if self.chromosome_type == "circular":
                theta_offset = (
                    (self.origin - 1) / (self.lengths.sum() * resolution) *
                    360)
                ax_g.set_theta_direction(-1)
            else:
                theta_offset = 0
            ax_g.set_theta_zero_location("N", offset=theta_offset)

        return ax_g


def _compute_rorigin(min_data, max_data, inner_radius, outer_radius):
    origin_data = (
        min_data - inner_radius * (max_data - min_data) /
        (outer_radius - inner_radius))
    return origin_data
