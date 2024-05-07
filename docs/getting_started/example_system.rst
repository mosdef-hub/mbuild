Example System
===============


Components in dashed boxes are drawn by hand using, e.g., `Avogadro <https://avogadro.cc>`_ or generated elsewhere.
`mBuild <https://mbuild.mosdef.org/en/stable/>`_ builds up complex systems from simple building blocks through simple attachment sites, called a ``Port`` (i.e., connection points). Each building block is a python class that can be customized or created through the pre-built options in the ``mBuild`` library ( ``mbuild.lib`` ). A hierarchical structure of parents and children is created through these classes, which can be easily parsed or modified.
This allows `mBuild <https://mbuild.mosdef.org/en/stable/>`_ to generate chemical structures in a piecemeal fashion by creating or importing molecular sections, adding ports, and connecting the ports to form bonds.
Together with `Signac <https://signac.io>`_, this functionality enables an automatic and dynamic method for generating chemical systems, allowing large-scale chemical and materials screening with minimal user interaction.

Ultimately, complex systems can be created with just a line or two
of code. Additionally, this approach seamlessly exposes tunable parameters within
the hierarchy so you can actually create whole families of structures
by adjusting a variable or two::

    pattern = Random2DPattern(20)  # A random arrangement of 20 pieces on a 2D surface.
    brush_layer = BrushLayer(chain_lenth=20, pattern=pattern, tile_x=3, tile_y=2)

.. figure:: ../images/pmpc.png
    :width: 100 %
    :align: center

    **Zwitterionic brushes on beta-cristobalite substrate.** Example system that can be created using mBuild.
    Components in dashed boxes are created from some external tool like Avogadro or SMILES strings.
    Components in solid boxes are created from these smaller dashed components and then constructed into larger,
    more complex systems using mBuild functionality.

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: http://opensource.org/licenses/MIT

Various sub-portions of this library may be independently distributed under
different licenses. See those files for their specific terms.
