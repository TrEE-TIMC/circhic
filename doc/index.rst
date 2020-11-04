.. circHiC documentation master file, created by
   sphinx-quickstart on Wed May  6 15:04:36 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

############################################################################
circHiC: circular visualization of Hi-C data and integration of genomic data
############################################################################

circHiC is a `Python <https://www.python.org>`_ library, built upon the widely
used `Matplotlib <https://matplotlib.org/>`_ library, to display Hi-C-like data in circular strips
together with the possibility to overlay genomic data (plots and heat maps).
Tools are light and fast, aiming to facilitate the exploration
and understanding of bacterial chromosome structuring data.

.. image:: auto_examples/biological_ccrescentus/images/sphx_glr_plot_ccres_mappability_001.png
	:width: 49.5%
	
.. image:: auto_examples/biological_ecoli/images/sphx_glr_plot_ecoli_macrodomain_001.png
	:width: 49.5%	
   	 
Genome wide contact frequencies obtained using Hi-C-like experiments have
raised novel challenges in terms of visualization and rationalization of
chromosome structuring phenomena. In bacteria, one difficulty consists in
displaying data in a way that is congruent with the circularity of
chromosomes. Standard representations of Hi-C data under the form of square
matrices or horizontal bands are indeed not adapted to periodic conditions as
those imposed by (most) bacterial chromosomes. circHiC fills this gap.

    
Citing circHiC
--------------

If you use circHiC in a scientific publication, we would appreciate citations
to the following paper:

  Junier, I., & Varoquaux, N. circHiC: circular visualization of Hi-C data and integration of genomic data. `bioRxiv <http://doi.org/10.1101/2020.08.13.249110>`_

Bibtex entry::

       @article {Junier2020.08.13.249110,
       		author = {Junier, Ivan and Varoquaux, Nelle},
       		title = {circHiC: circular visualization of Hi-C data and integration of genomic data},
       		elocation-id = {2020.08.13.249110},
       		year = {2020},
       		doi = {10.1101/2020.08.13.249110},
       		publisher = {Cold Spring Harbor Laboratory},
       		URL = {https://www.biorxiv.org/content/early/2020/08/14/2020.08.13.249110},
       		eprint = {https://www.biorxiv.org/content/early/2020/08/14/2020.08.13.249110.full.pdf},
       		journal = {bioRxiv}}


Licence Information
--------------------

Conditions on the use and redistribution of this package.

.. literalinclude:: ../LICENSE.txt

