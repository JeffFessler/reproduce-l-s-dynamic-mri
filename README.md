# reproduce-l-s-dynamic-mri

https://github.com/JeffFessler/reproduce-l-s-dynamic-mri

Matlab code for reproducing the results in the paper:
"Efficient Dynamic Parallel MRI Reconstruction
for the Low-Rank Plus Sparse Model,"
IEEE Trans. on Computational Imaging,
5(1):17-26, 2019,
by Claire Lin (mailto::yililin@umich.edu) and Jeffrey A. Fessler,
EECS Department, University of Michigan

See
[http://doi.org/10.1109/TCI.2018.2882089]

To download:
`git clone https://github.com/JeffFessler/reproduce-l-s-dynamic-mri.git`

This code requires the Matlab version of
the Michigan Image Reconstruction Toolbox (MIRT)
from
[http://web.eecs.umich.edu/~fessler/code/index.html]

Please set up MIRT before running the examples.

The authors also thank Ricardo Otazo
for sharing the cardiac MRI data.

Before running the code you must populate the data/ directory
with five .mat files from
[http://web.eecs.umich.edu/~fessler/irt/reproduce/19/lin-19-edp/]

The data is stored separately to keep the git repo light.

The following scripts are in the example folder:
* `example_PINCAT_phantom.m` Figs. 4, 5
* `example_abdomen_dce_ga.m` Supplement Fig. 6 (non-Cartesian)
* `example_cardiac_cine.m` Figs. 1, 3
* `example_cardiac_perf.m` Figs. 1, 2


## Julia version

There is a
[Julia language version](https://juliaimagerecon.github.io/Examples/generated/mri/5-l-plus-s)
in the
[Examples part of JuliaImageRecon](https://github.com/JuliaImageRecon/Examples).
