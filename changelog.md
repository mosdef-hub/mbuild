# Change Log #

All big and breaking changes for [mBuild](https://mosdef-hub.github.io/mbuild/) will be recorded here.
This project adheres to [Semantic Versioning](http://semver.org/).

## 0.9.1 (unreleased)
### Breaking Changes
### Features
### Misc and Bugfixes

## 0.8.3 (unreleased)
### Breaking Changes
* When writing hoomdxml files, units will now be in kJ/mol & nm instead of kcal/mol & ang, so particle positions will be different by a factor of 10.
### Features
* When saving hoomdxml files, `auto_scale=True` will scale reference units from max forcefield parameters.
### Misc and Bugfixes

## 0.8.2 (2018-1-8)
### Features
* Special Pair Support (1-4 pair information) to GSD writers (#473)
    * GSD files now include 1-4 special pairs for use in OPLS

### Misc and Bugfixes
* Dependency requirements have been updated (#457)
    * A dependency loop between `Foyer` and `mBuild` has been resolved
* Fixed a bug that prevented Appveyor builds from running (#477)
* Temporary PDB files left behind by packing functions are now properly removed (#471)
    * Packing.py uses temporary files which were previously never closed. This sometimes caused the program to reach the limit of open files for a process set by the OS
* `pytest-ignore-flaky` has been replaced in favor of `xfail` (#471)
* Additonal fixes for PACKMOL input files (#474)
    * Input files are now closed by `mBuild` in order to ensure it can be read by PACKMOL
    * Error reporting is now caught when the subprocess returns an error code
* Microsoft VSCode extraneous files are now ignored by git (#478)

## 0.8.1 (2018-11-28)
### Features
* Packing functions can optionally constrain the rotation of `Compounds` when using `fill_box` (#407)
* Additional lammps datafile support (#412)
    * Add functionality for `atomic`, `charge`, and `molecular` atom_styles
    * Fix `atomic` and `molecular` atom-styles
    * Add optional `atom-style` argument to `save` function
    * Add tests to check for correct `Atoms` format
* A `Compound` can be generated from a SMILES string if the user has [Open Babel](http://openbabel.org/wiki/Main_Page) installed (#430)
* The [website](https://mosdef-hub.github.io/mbuild/) was updated with details how to properly cite mBuild (#421)
* OpenMM can now be used for energy minimization (#416)
* A simple xyz file reader was added (#423)
* Defaults in `Compound.visualize` have been improved (#438)
* mBuild boxes can now store angles (#448)
* mBuild boxes can now be passed to various writers (#448)
* A changelog is now included in the root directory (#454)

### Misc and Bugfixes
* Switched from OpenMM to MDTraj to check for element existence (#408)
* Changed bilayer notebook to use `Compound` methods for object manipulation (#406)
* Added test to ensure that users can provide a custom forcefield XML file when applying a forcefield to a `Compound` (#431)
* An error is now generated if the miniconda MD5 sum does not match when performing tests (#409)
* LAMMPS box values are now written appropriately in Angstroms (#459)
* Coordinates in HOOMDXML and GSD files are now correctly written in the range [L /2 to L] (#452)
* A bug in the ordering of some Bravais angles in non-rectangular lattices has been fixed (#450)

## 0.8.0 (2018-01-18)
### Features
* Improved packing API
	* `fill_box` method now supports user-specified densities (#372)
	* Support for non-cubic boxes (#385)
	* Added edge buffer for pseudo-support of periodic boundaries (#385)
	* Allow users to access PACKMOL raw output (#385)
	* Improve documentation and removed repetitious code (#385)
* Proper support for triclinic lattices (#386)
* More intuitive `Port` behavior when adding/removing bonds
	* `Ports` are now added along the bond vector when bonds are removed
	* `Ports` are removed when using `force_overlap` when attribute `add_bond` is `True` (#390)
### Misc and Bugfixes
* Increased precision of PACKMOL overlap tolerance (#367)
* Increased robustness for `packing.py` argument types (#368)
* Continuous Integration (CI) fixes (#369, #375, #392, #402, #405)
* Documentation updates (#371)
* Combining rules for non-bonded interactions can now be specified when saving compounds (#377)
* Remove ambigious data types for `Box` attributes (#384)
* Fixed changing of basis for non-cubic lattices (#386)
* Fixed issue where `Ports` were not aligned properly along user-specified direction (#390)
* Add support for controlling `Foyer` warning verbosity when saving `Compounds` (#391)
* Added more information to `Port.__repr__` (#396)
* Fixed bug in unit conversion for periodicity in `Compound.from_parmed()` (#401)
* Rounding `Lattice.lattice_point` positions to 0. when below a certain threshold (#404)

## 0.7.1 (2017-06-09)

Merge pull request (#352) from summeraz/use_parmed
Parmed loaders by default, remove Mdtraj dependency

## 0.7.0 (2017-06-09)

Bump version to 0.7.0

## 0.6.1 (2017-02-14)

Tag 0.6.1

## 0.6.0 (2017-02-14)

Tag 0.6.0

## 0.5.2 (2015-08-27)

Merge pull request (#131) from ctk3b/dev-package
Fix requirements.
