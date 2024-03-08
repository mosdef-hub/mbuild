Box
========================


Import the required mbuild package.

.. code:: ipython

    import mbuild as mb


Orthogonal Box
------------------------

Build an empty orthogonal **mBuild** Box (i.e., the angle in degrees are 𝛼 = 90, 𝛽 = 90, 𝛾 = 90) measuring 4.0 nm in all the x, y, and z-dimensions.

.. note::
    Note: if the angles are not specified, the system will default to an orthogonal box
    (i.e., the angle in degrees are 𝛼 = 90, 𝛽 = 90, 𝛾 = 90).

.. code:: ipython

    empty_box = mb.Box(lengths=[4.0, 4.0, 4.0], angles=[90, 90, 90])


Non-Orthogonal Box
------------------------

Build an empty non-orthogonal **mBuild** Box (i.e., the angle in degrees are 𝛼 = 90, 𝛽 = 90, 𝛾 = 120) measuring 4.0 nm in the x and y-dimensions, and 5.0 nm in the z-dimension.

.. code:: ipython

    empty_box = mb.Box(lengths=[4.0, 4.0, 5.0], angles=[90, 90, 120])
