===============
Data Structures
===============

The primary building blocks in an mBuild hierarchy inherit from the
:py:class:`mbuild.Compound` class.  ``Compounds`` maintain an ordered set of ``children``
which are other ``Compounds``.  In addition, an independent, ordered dictionary
of ``labels`` is maintained through which users can reference any other
``Compound`` in the hierarchy via descriptive strings.  Every ``Compound``
knows its parent ``Compound``, one step up in the hierarchy, and knows which
``Compounds`` reference it in their ``labels``.  :py:class:`mbuild.Port` is a special type
of ``Compound`` which are used internally to connect different ``Compounds``
using the equivalence transformations described below.

``Compounds`` at the bottom of an mBuild hierarchy, the leaves of the tree, are
referred to as ``Particles`` and can be instantiated as ``foo =
mbuild.Particle(name='bar')``. Note however, that this merely serves to illustrate
that this ``Compound`` is at the bottom of the hierarchy; ``Particle`` is simply
an alias for ``Compound`` which can be used to clarify the intended role of an
object you are creating. The method :py:meth:`mbuild.Compound.particles` traverses the
hierarchy to the bottom and yields those ``Compounds``. :py:meth:`mbuild.Compound.root`
returns the compound at the top of the hierarchy.

Compound
--------

.. autoclass:: mbuild.Compound
	:members:


Box
---

.. autoclass:: mbuild.Box
	:members:

Lattice
-------

.. autoclass:: mbuild.Lattice
	:members:


Port
----

.. autoclass:: mbuild.Port
	:members:
