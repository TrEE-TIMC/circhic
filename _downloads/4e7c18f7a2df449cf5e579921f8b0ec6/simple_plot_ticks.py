import circhic
import numpy as np

from iced.normalization import ICE_normalization


# Start by loading the data
data = circhic.datasets.load_ecoli()
counts = data["counts"]
nbins = data["nbins"]

# Normalize the data using ICE, and keep the biases
counts, bias = ICE_normalization(counts, output_bias=True)


# Now instantiate the circhic Figure
circhicfig = circhic.CircHiCFigure(lengths=nbins)

# Compute the extreme values
vmax = np.max([counts[i, (i+1) % counts.shape[0]]
               for i in range(counts.shape[0])])
vmin = np.min(counts[counts > 0]) * 10

# define the inner genomid distances and the outer genomic distance plotted
inner_radius = 0.1
inner_gdis, outer_gdis = 200, 60

im, ax = circhicfig.plot_hic(counts, cmap="bone_r", border_thickness=0.01,
                             vmin=vmin, vmax=vmax, inner_radius=inner_radius,
                             inner_gdis=inner_gdis, outer_gdis=outer_gdis)

# Add the colorbar as a vertical colorbar
cab = circhicfig.set_colorbar(im, orientation="horizontal")
cab.set_label("Normalized contact counts", fontweight="bold", color="0.3")

# Now, let us add the ticks and tick labels
ticklabels = ["%dkb" % (i * 500) for i in range(9)]
tickpositions = [int(i*50) for i in range(9)]
ticklabels[0] = "OriC"
ax = circhicfig.set_genomic_ticklabels(
    tickpositions=tickpositions,
    ticklabels=ticklabels,
    outer_radius=1, fontdict={'fontsize': "small"})
ax.tick_params(colors="0.3")

rax = circhicfig.plot_raxis()
rax.set_yticklabels(["200kb", "0kb", "60kb"], fontsize="small")

rax.set_ylabel("Genomic distance", color="0.3",
               fontweight="bold")
rax.tick_params(colors="0.3")
