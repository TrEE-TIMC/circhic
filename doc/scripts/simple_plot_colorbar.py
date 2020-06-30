import circhic

# Start by loading the data
data = circhic.datasets.load_ecoli()
counts = data["counts"]
nbins = data["nbins"]

# Now instantiate the circhic Figure
circhicfig = circhic.CircHiCFigure(lengths=nbins)
im, ax = circhicfig.plot_hic(counts, cmap="bone_r", border_thickness=0.01)

# Add the colorbar as a vertical colorbar
cab = circhicfig.set_colorbar(im, orientation="horizontal")
cab.set_label("Normalized contact counts", fontweight="bold", color="0.3")
