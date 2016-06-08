==============
Data Structure
==============

The primary building blocks in an mBuild hierarchy inherit from the
``Compound`` class.  ``Compounds`` maintain an ordered set of ``children``
which are other ``Compounds``.  In addition, an independent, ordered dictionary
of ``labels`` is maintained through which users can reference any other
``Compound`` in the hierarchy via descriptive strings.  Every ``Compound``
knows its parent ``Compound``, one step up in the hierarchy, and knows which
``Compounds`` reference it in their ``labels``.  ``Ports`` are a special type
of ``Compound`` which are used internally to connect different ``Compounds``
using the equivalence transformations described below.

``Compounds`` at the bottom of an mBuild hierarchy, the leafs of the tree, are
referred to as ``Particles`` and can be instantiated as ``foo =
mb.Particle(name='bar')``. Note however, that this merely serves to illustrate
that this ``Compound`` is at the bottom of the hierarchy; ``Particle`` is simply
an alias for ``Compound`` which can be used to clarify the intended role of an
object you are creating. The method ``Compound.particles()`` traverses the
hierarchy to the bottom and yields those ``Compounds``. ``Compound.root()``
returns the compound at the top of the hierarchy.

Compound
--------

.. autoclass:: mbuild.compound.Compound
    :members:

Port
----

.. autoclass:: mbuild.port.Port
    :members:

