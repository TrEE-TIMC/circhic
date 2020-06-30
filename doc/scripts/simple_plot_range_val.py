import circhic

# Start by loading the data
data = circhic.datasets.load_ecoli()
counts = data["counts"]
nbins = data["nbins"]

# Now instantiate the circhic Figure
circhicfig = circhic.CircHiCFigure(lengths=nbins)

# Compute the extreme values
vmax = np.max([counts[i, (i+1) % counts.shape[0]]
               for i in range(counts.shape[0])])
vmin = np.min(counts[counts > 0])

im, ax = circhicfig.plot_hic(counts, cmap="bone_r", border_thickness=0.01,
                             vmin=vmin, vmax=vmax)

# Add the colorbar as a vertical colorbar
cab = circhicfig.set_colorbar(im, orientation="horizontal")
cab.set_label("Normalized contact counts", fontweight="bold", color="0.3")
