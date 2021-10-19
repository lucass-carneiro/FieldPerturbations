# Cactus Framework Arrangement "FieldPerturbations"
## Author(s):
* Lucas Timotheo Sanches <lucas.t@ufabc.edu.br>

## Maintainer(s):
* Lucas Timotheo Sanches <lucas.t@ufabc.edu.br>

## Licence:
GNU GPLv3+ (see <https://www.gnu.org/licenses/>).

## Purpose
This arrangement contains [EinsteinToolkit](https://einsteintoolkit.org/index.html) thorns for the evolution of perturbation fields in arbitrary background geometries. Here's a brief description of the current thorns and their functionalities:

1. `KleinGordon`: Evolves the Klein-Gordon equation on top of an arbitrary background, without taking into account the geometry's back-reaction and thus having a null contribution to the energy-momentum tensor. This thorn is compatible with the [Carpet](https://bitbucket.org/eschnett/carpet/src/master/) AMR infrastructure.
2. `KleinGordonX`: [CarpetX](https://bitbucket.org/eschnett/cactusamrex/src/master/) compatible version of `KleinGordon`. This thorn sets the energy-momentum tensor and it should be possible to perform full evolutions using it. Note however that the energy-momentum tensor support was added very recently and bugs are to be expected.
3. `KerrSchildX`: This is an adaptation of the [EinsteinExact](https://github.com/barrywardell/EinsteinExact) Kranc-generated thorn `KerrSchild` for the `CarpetX` driver. I'ts used to setup single BH geometric initial data for the perturbation thorns. If the user desires, this thorn also enforces the exact solution throughout the evolution.
4. `MinkowskiX`: This thorn is similar in nature to the [EinsteinExact](https://github.com/barrywardell/EinsteinExact) Kranc-generated thorn `Minkowski` for the `CarpetX` driver. This thorn is used to setup Minkowski initial data for the perturbation thorns. If the user desires, this thorn also enforces the exact solution throughout the evolution.

All the tensor quantities in this thorn were expanded/computed with the help of Wolfram Mathematica. The equations implemented can be found on the compressed notebook file *equations.nb.gz*

The thorns in this branch use multipatch infrastructures. Here's a status of the multipatch conversion for each thorn:

1. `KleinGordon`: Not yet started.
2. `KleinGordonX`: Not yet started. CarpetX patch system not yet implemented.
3. `KerrSchildX`: Not yet started. CarpetX patch system not yet implemented.
4. `MinkowskiX`: Not yet started. CarpetX patch system not yet implemented.

**WARNING:** These thorns are under heavy development and testing.
