"""Coarse-graining and backmapping tools built on CGsmiles.

This package bridges mBuild and the CGsmiles library
(https://github.com/gruenewald-lab/CGsmiles). It operates on
coarse-grained systems described either by an ``mbuild.path.Path``
(bead names, coordinates, and a bond graph) or by an
``mbuild.Compound`` whose leaf particles are the CG beads (e.g. a CG
structure loaded from a file or built by hand).

- ``to_cgsmiles_graph`` / ``to_cgsmiles`` — express the CG system as a
  CGsmiles meta graph or string.
- ``resolve`` — resolve each bead to molecular detail using CGsmiles
  fragment strings (SMILES + bonding descriptors) and/or tagged mBuild
  compounds.
- ``backmap`` — the full pipeline: returns an atomistic
  ``mbuild.Compound`` that retains the CG conformation. Works for any
  bond graph topology, including branch points and crosslinks that
  ``Polymer.build_from_path`` cannot handle.
- ``compound_to_fragment_graph`` — let a tagged mBuild Compound define
  a CGsmiles fragment.

If you use this feature please cite:

CGsmiles: A Versatile Line Notation for Molecular Representations
across Multiple Resolutions. Fabian Grünewald, Leif Seute, Riccardo
Alessandri, Melanie König, and Peter C. Kroon. Journal of Chemical
Information and Modeling 2025 65 (7), 3405-3419.
DOI: 10.1021/acs.jcim.5c00064
"""

# ruff: noqa: F401
from .backmap import backmap, resolve
from .cg_map import CGMapping, coarse_grain
from .convert import to_cgsmiles, to_cgsmiles_graph
from .fragments import compound_to_fragment_graph

__all__ = [
    "backmap",
    "coarse_grain",
    "CGMapping",
    "resolve",
    "to_cgsmiles",
    "to_cgsmiles_graph",
    "compound_to_fragment_graph",
]
