import circhic

# Start by loading the data
data = circhic.datasets.load_ecoli()
counts = data["counts"]
nbins = data["nbins"]

# Now instantiate the circhic Figure
circhicfig = circhic.CircHiCFigure(lengths=nbins)
circhicfig.plot_hic(counts, cmap="bone_r", border_thickness=0.01)
