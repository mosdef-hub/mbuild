=====
Units
=====

mBuild automatically performs unit conversions in its reader and writer functions.
When working with an :py:class:`mbuild.Compound`, mBuild uses the following units:

+----------+----------+
| Quantity |   Units  |
+==========+==========+
| distance |    nm    |
+----------+----------+
|   angle  | radians* |
+----------+----------+

\* :py:class:`mbuild.Lattice` and :py:class:`mbuild.Box` use degrees.

See also `foyer unit documentation <https://foyer.mosdef.org/en/stable/units.html>`_ and `ele documentation <https://ele-ment.readthedocs.io/en/latest/>`_.
