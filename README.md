# muMatScale

Code for modeling solidification microstructure based on the Cellular Automata method for binary/pseudo-binary alloys.

This code is derived from the open-source software uMatIC developed at the
Imperial College London
https://www.imperial.ac.uk/engineering-alloys/research/software

## Authors ##
(Authors for the original version are not listed)

Lang Yuan, University of South Carolina, Email: langyuan@cec.sc.edu

## Contributors ##
(Contributors for the original version are not listed)

Jean-Luc Fattebert, Oak Ridge National Laboratory, Email: fattebertj@ornl.gov

Adrian Sabau, Oak Ridge National Laboratory, Email: sabaua@ornl.gov

## Dependencies ##

* [HDF5](https://support.hdfgroup.org/HDF5)

## References ##


The OpenMP offload implementation for GPU is described in:

A.S. Sabau, L. Yuan, J.-L. Fattebert, J.A. Turner,
"An OpenMP GPU-offload implementation of a non-equilibrium solidification
cellular automata model for additive manufacturing",
Computer Physics Communications 284 (2023), 108605
https://doi.org/10.1016/j.cpc.2022.108605

The model implemented in muMatScale and its example applications are described in the following reference:

W. Wang, P.D. Lee, M. McLean,
"A model of solidification microstructures in nickel-based
superalloys: predicting primary dendrite spacing selection",
Acta Materialia 51 (2003), 2971-2987
https://doi.org/10.1016/S1359-6454(03)00110-1

L. Yuan, P. D. Lee,
"A new mechanism for freckle initiation based on microstructural level simulation",
Acta Materialia, 60 (2012), 4917-4926
https://doi.org/10.1016/J.ACTAMAT.2012.04.043

A. Prasad, L. Yuan, P.D. Lee, M. Patel, D. Qiu, M. Easton, D. StJohn
"Towards understanding grain nucleation under Additive Manufacturing solidification conditions",
Acta Materialia, 195 (2020), 392-403
https://doi.org/10.1016/j.actamat.2020.05.012


## License ##

muMatScale is released under the BSD 3-clause license. See the [LICENSE](./LICENSE)
and the original license from uMatIC [uMatIC_LICENSE](./uMatIC_LICENSE).
All new contributions must be made under the BSD 3-clause license.
