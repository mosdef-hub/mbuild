
Building a Simple Alkane
========================

The purpose of this tutorial is to demonstrate the construction of an
alkane polymer and provide familiarity with many of the underlying
functions in mBuild. Note that a robust polymer construction recipe
already exists in mBuild, which will also be demonstrated at the end of
the tutorial.

Setting up the monomer
----------------------

The first step is to construct the basic repeat unit for the alkane,
i.e., a :math:`CH_2` group, similar to the construction of the
:math:`CH_3` monomer in the prior methane tutorial. Rather than
importing the coordinates from a pdb file, as in the previous example,
we will instead explicitly define them in the class. Recall, that
distance units are nm in mBuild.

.. code:: ipython3

    import mbuild as mb
    
    class CH2(mb.Compound):
        def __init__(self):
            super(CH2, self).__init__()
    
            self.add(mb.Particle(name='C', pos=[0,0,0]), label='C[$]')
    
            self.add(mb.Particle(name='H', pos=[-0.11, 0, 0.0]), label='HC[$]')    
            self.add(mb.Particle(name='H', pos=[0.11, 0, 0.0]), label='HC[$]')
            
            port1 = mb.Port(anchor=self[0])
            self.add(port1, label='up')
            # Place the port at approximately half a C-C bond length.
            mb.translate(self['up'], [0, -0.07, 0]) 
            
            port2 = mb.Port(anchor=self[0])
            self.add(port2, label='down')
            # Place the port at approximately half a C-C bond length.
            mb.translate(self['down'], [0, 0.07, 0]) 
    
    monomer = CH2()
    monomer.visualize(show_ports=True)

the configuration of the monomer is not a particularly realistic
conformation. One could use this monomer to construct a polymer and then
apply an energy minimization scheme, or, as we will demonstrate here, we
can mBuild's rotation commands to provide a more realistic starting
point.

Below, we use the same basic script, but now apply a rotation to the
hydrogen atoms. Since the hydrogens start 180 degrees apart and we know
they should be ~108 degrees, each should be rotated approximately 35
degrees closer to each other around the 'y' axis, i.e., ~0.63 radians.

Similarly, the ports should be rotated around the 'x' axis by ~0.59
radians to produce a ~112 degree angle in the backbone.

.. code:: ipython3

    import mbuild as mb
    
    class CH2(mb.Compound):
        def __init__(self):
            super(CH2, self).__init__()
    
            self.add(mb.Particle(name='C', pos=[0,0,0]), label='C[$]')
    
            self.add(mb.Particle(name='H', pos=[-0.11, 0, 0.0]), label='HC[$]')    
            self.add(mb.Particle(name='H', pos=[0.11, 0, 0.0]), label='HC[$]')
            mb.rotate_around_y(self['HC'][0], 0.63)
            mb.rotate_around_y(self['HC'][1], -0.63)
            
            port1 = mb.Port(anchor=self[0])
            self.add(port1, label='up')
            # Place the port at approximately half a C-C bond length.
            mb.translate(self['up'], [0, -0.07, 0]) 
            mb.rotate_around_x(self['up'], 0.59)
    
    
            port2 = mb.Port(anchor=self[0])
            self.add(port2, label='down')
            # Place the port at approximately half a C-C bond length.
            mb.translate(self['down'], [0, 0.07, 0]) 
            mb.rotate_around_x(self['down'], -0.59) 
    
    
    monomer = CH2()
    monomer.visualize(show_ports=True)

Defining the polymerization class
---------------------------------

With a basic monomer construct, we can now construct a polymer by
connecting the ports together. Note, here, we first instantiate one
instance of the CH2 class (as last\_monomer), then use the clone
function to make a copy. The force\_overlap function is used to connect
the 'up' port from the current\_monomer to the 'down' port of the
last\_mononer.

.. code:: ipython3

    class AlkanePolymer(mb.Compound):
        def __init__(self):
            super(AlkanePolymer, self).__init__()
    
            last_monomer = CH2()
            self.add(last_monomer)
    
            for i in range (3):
                current_monomer = mb.clone(last_monomer)
        
                mb.force_overlap(move_this=current_monomer, 
                                 from_positions=current_monomer['up'], 
                                 to_positions=last_monomer['down'])
                self.add(current_monomer)
                last_monomer=current_monomer
    
    polymer = AlkanePolymer()
    polymer.visualize(show_ports=True)

Visualization of this structure should demonstrate a problem; the
polymer curls up on itself. This is a result of the fact that ports not
only define the location in space, but orientation. This can be
trivially fixed, by first flipping the port over (i.e., rotate around
'y' by 180 degrees, ~3.14 radians).

We can also add a variable 'chain\_length' both to the for loop and
'init', that will allow the length of the polymer to be adjusted when
the class is instantiated.

.. code:: ipython3

    import mbuild as mb
    
    class CH2(mb.Compound):
        def __init__(self):
            super(CH2, self).__init__()
    
            self.add(mb.Particle(name='C', pos=[0,0,0]), label='C[$]')
    
            self.add(mb.Particle(name='H', pos=[-0.11, 0, 0.0]), label='HC[$]')    
            self.add(mb.Particle(name='H', pos=[0.11, 0, 0.0]), label='HC[$]')
            mb.rotate_around_y(self['HC'][0], 0.63)
            mb.rotate_around_y(self['HC'][1], -0.63)
            
            port1 = mb.Port(anchor=self[0])
            self.add(port1, label='up')
            # Place the port at approximately half a C-C bond length.
            mb.translate(self['up'], [0, -0.07, 0]) 
            mb.rotate_around_x(self['up'], 0.59)
    
    
            port2 = mb.Port(anchor=self[0])
            self.add(port2, label='down')
            # Place the port at approximately half a C-C bond length.
            mb.translate(self['down'], [0, 0.07, 0]) 
            mb.rotate_around_y(self['down'], 3.14) 
            mb.rotate_around_x(self['down'], -0.59) 
    
    
    
    class AlkanePolymer(mb.Compound):
        def __init__(self, chain_length=1):
            super(AlkanePolymer, self).__init__()
    
            last_monomer = CH2()
            self.add(last_monomer)
    
            for i in range (chain_length-1):
                current_monomer = mb.clone(last_monomer)
        
                mb.force_overlap(move_this=current_monomer, 
                                 from_positions=current_monomer['up'], 
                                 to_positions=last_monomer['down'])
                self.add(current_monomer)
                last_monomer=current_monomer

.. code:: ipython3

    polymer = AlkanePolymer(chain_length = 10)
    polymer.visualize(show_ports=True)

Using mBuild's Polymer Class
----------------------------

mbuild provides a prebuilt class to perform this basic functionality.
Since it is designed to be more general, it takes as an argument not
just the chain length, but also the monomer and the port labels (i.e.,
'up' and 'down', since these labels are user defined).

.. code:: ipython3

    polymer = mb.Polymer(CH2(), 10, port_labels=('up', 'down'))
    polymer.visualize()

Building a System of Alkanes
----------------------------

A system of alkanes can be constructed by simply cloning the polymer
constructed above and translating the alkanes in space (and/or rotating
them). mBuild provides many routines that can be used to create
different patterns, to which the polymers can be shifted to.

.. code:: ipython3

    #define a compound to hold all the polymers
    system = mb.Compound()
    
    
    #we will first rotate the chain such that the long dimension is oriented along the z-axis
    mb.rotate_around_x(polymer, 3.14/2.0)
    
    pattern_disk = mb.DiskPattern(50)
    
    #patterns are generated betwee 0 and 1, 
    #and thus need to be scaled to provide appropriate spacing
    pattern_disk.scale(5) 
    
    for pos in pattern_disk:
        
        current_polymer = mb.clone(polymer)
        mb.translate(current_polymer, pos)
        system.add(current_polymer)

.. code:: ipython3

    system.visualize()

Other patterns can be used, e.g., the Grid3DPattern. We can also use the
rotation commands to randomize the orientation.

.. code:: ipython3

    #define a compound to hold all the polymers
    system = mb.Compound()
    import random
    
    #we will first rotate the chain such that the long dimension is oriented along the z-axis
    mb.rotate_around_x(polymer, 3.14/2.0)
    
    pattern_disk = mb.Grid3DPattern(5,5,5)
    
    
    for pos in pattern_disk:
        
        #scale the pattern 
        pos[0] = pos[0]*8.0
        pos[1] = pos[1]*8.0
        pos[2] = pos[2]*8.0
    
        
        current_polymer = mb.clone(polymer)
        #give the polymers random perturbation from their initial orientation
        mb.rotate_around_z(current_polymer, random.uniform(0,3.14))
        mb.rotate_around_x(current_polymer, random.uniform(0,3.14))
        mb.rotate_around_z(current_polymer, random.uniform(0,3.14))
    
    
    
        mb.translate(current_polymer, pos)
        system.add(current_polymer)

.. code:: ipython3

    system.visualize()

mBuild also provides an interface to Packmol, allowing the creation of a
randomized configuration.

.. code:: ipython3

    polymer = mb.Polymer(CH2(), 5, port_labels=('up', 'down'))
    system = mb.fill_box(polymer, n_compounds=100, overlap=1.5, box=[10,10,10]) 

.. code:: ipython3

    system.visualize()

Variations
----------

Rather than a linear chain, the Polymer class we wrote can be easily
changed such that small perturbations are given to each port. To avoid
accumulation of deviations from the equilibrium angle, we will clone an
unperturbed monomer each time (i.e., monomer\_proto) before applying a
random variation.

We also define a variable 'delta' which will control the maximum amount
of perturbation, where clearly a smaller value results in a more linear
conformation. Note, large values may result in the chain overlapping
itself, as mBuild does not currently include routines to exclude such
overlaps.

.. code:: ipython3

    import mbuild as mb
    
    import random
    
    class AlkanePolymer(mb.Compound):
        def __init__(self, chain_length=1, delta=0):
            super(AlkanePolymer, self).__init__()
    
            monomer_proto = CH2()
    
            last_monomer = CH2()
            mb.rotate_around_x(last_monomer['down'], random.uniform(-delta,delta))
            mb.rotate_around_y(last_monomer['down'], random.uniform(-delta,delta))
    
    
            self.add(last_monomer)
    
            for i in range (chain_length-1):
                current_monomer = mb.clone(monomer_proto)
                mb.rotate_around_x(current_monomer['down'], random.uniform(-delta,delta))
                mb.rotate_around_y(current_monomer['down'], random.uniform(-delta,delta))
    
                mb.force_overlap(move_this=current_monomer, 
                                 from_positions=current_monomer['up'], 
                                 to_positions=last_monomer['down'])
                self.add(current_monomer)
                last_monomer=current_monomer

.. code:: ipython3

    polymer = AlkanePolymer(chain_length = 200, delta=0.4)
    polymer.visualize()


