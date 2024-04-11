![spac3tools](./images/spac3tools_2x.gif)

# spac3tools

The R package `spac3tools` contains a suite of landscape or space
 (environmental dynamic layers) that are used as input for the general
 engine for eco-evolutionary simulations
 [`gen3sis`](https://github.com/project-gen3sis/R-package).
 
This package can:

1.  convert `landscapes.rds` (`gen3sis`) to `space.rds` (`gen3sis_rf`)

2.  compress a `gen3sis_space` object, i.e. `spaces.rds` for optimal storage

3.  create `gen3sis_spaces` object and access utility functions for plotting and input creation

4.  decompress a spaces.rds (creates distances) prior to running a simulation
