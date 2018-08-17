# Change Log #

All big and breaking changes for [mBuild](https://mosdef-hub.github.io/mbuild/) will be recorded here.
This project adheres to [Semantic Versioning](http://semver.org/).

## 0.9.1 (unreleased)
### Breaking Changes
### Features
### Misc and Bugfixes

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
