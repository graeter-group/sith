# SITH
Splitting Intramolecular Tension via stretcHing

SITH is a novel method that decomposes the total electronic energy change of a
stretched molecule into contributions from individual degrees of freedom —such
as bond lengths, angles, and dihedrals— using numerical integration of the
work-energy theorem.

API is documented and examples are attached, but full documentation of the sith
package is still active, albeit slowly. Thank you for your patience and please
do not hesitate to direct any questions to the maintainer at
daniel.sucerquia@h-its.org, mfarrugi@nd.edu.

## Version 2.0
- Every energy distribution analysis method is a different module in the package.
- using sith is now divided in two steps. Reading and computing QM data by sith.readers. Reading is done to allow any DFT software to do the computation. So far, only g09 files are included, but the idea is to include other sofwares (QChem, ORCA and GPAW).
- fchk reader methods are now in sith.g09_reader. This will allow users of others softwares add their own readers.
- in sith.g09_reader.FileReader, every block defined by ..._header is read by _fill_aray method or directly in one line (this avoided long if else blocks).
- sith.Geometry.dim_indices.shape = (#dofs, 4) and dim_indices is now an numpy array. Distances, that only have two indices, the rest are zeros (remember g09 counts from 1). indeed all the arrays have that shape.
- Geometry.build_RIC was completely removed. This is done by readers now.
- Lmatrix class is not necessary anymore, completely removed and replaced by _hessian_reader in g09_reader
- UnitConverter class in .Utilities is not necessary, completely removed and replaced by ase.units
- all attributes called "._deformed" renamed as "structures"
