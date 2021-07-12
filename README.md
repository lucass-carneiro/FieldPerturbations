# Cactus Framework Arrangement "FieldPerturbations"
## Author(s):
* Lucas Timotheo Sahches <lucas.t@ufabc.edu.br>

## Maintainer(s):
* Lucas Timotheo Sahches <lucas.t@ufabc.edu.br>

## Licence:
GNU GPLv3+ (see <https://www.gnu.org/licenses/>).

## Purpose
This arrangement contains Cactus thorns for the evolution of perturbation fields in arbitrary background geometries. Here's a brief description of the current thorns and their functionalities:

1. KleinGordon: Evolves the Klein-Gordon equation on top of an arbitrary background, without taking into account the geometry's back-reaction and thus having a null contribution to the energy-momentum tensor. This thorn is compatible with the [Carpet](https://bitbucket.org/eschnett/carpet/src/master/) AMR infrastructure.
2. KleinGordonX: [CarpetX](https://bitbucket.org/eschnett/cactusamrex/src/master/) compatible version of KleinGordon.

All the tensorial quantities in this thorn were expanded/computed with the help of Wolfram Mathematica. The equations implemented can be found on the compressed notebook file *equations.nb.gz*


**WARNING:** These thorns are under heavy development and testing. Everything is still subject to change (even their names).
