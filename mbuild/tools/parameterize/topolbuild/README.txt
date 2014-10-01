                              topolbuild
                              Version 1.3
                            September 1, 2009
                           Bruce D. Ray, Ph.D.
                           IUPUI Physics Dept.
                           402 N. Blackford St.
                          Indianapolis, IN 46202
                                   USA
USAGE:
     topolbuild -n modelname -dir path_to_interp. -ff forcefield [-r resid_name] \
                [-last ##] [-resnum ##] [-renumall] [-meas] [-gro] [-rename] \
                [-purge [level]] [-move] [-charge] [-h] [-v]


DESCRIPTION:
Read a Tripos .mol2 file with charges and generate something approximating
a gromacs .gro and .itp from it based on a force field for which atom type
definitions tables and suitable force field parameter files are available.
Requires that the .mol2 file have syntactically correct Tripos atom types.
The user is advised that while the format of the topology will be correct,
this does not guarantee correct parameterization for a molecule that does
not have any similarity to the molecules for which the chosen force field's
parameters were developed.  Please check parameterization and use the
topology produced by any automatic topology generator, not just topolbuild,
with discretion, advanced knowledge of molecular dynamics, understanding of
the issues involved in parameterization, and caution.
Note that the standard for a syntactically correct .mol2 file is the Tripos
Toolkit Utilities Manual.
Uses the following parameters and flags:
     -n                Required common name for the .mol2, .gro, .itp, .top,
                       and .log files (the path to and name of the input file
                       stripped of the .mol2) a revised .mol2 file will be
                       output to common_nameMOL.mol2 as well.

     -dir              Required absolute path to a directory that contains
                       appropriate force field data.  For amber type force
                       fields including glycam, this would be a path to a
                       directory that contains directories dat/antechamber
                       and dat/leap/parm (these in turn must contain the
                       appropriate antechamber and leap data files).
                       For tripos forcefield this would be the value of Tripos
                       variable $TA_ASCTABLES
                       For gromacs force fields, this would be a directory
                       that contains ATOMTYPE_GMX1.DEF, ATOMTYPE_gmx####.def,
                       ffgmx####.dat (where the #### represents the particular
                       force field's designation), ATOMTYPE_OPLSAA1.DEF,
                       ATOMTYPE_oplsaa.DEF and ffoplsaa.dat

     -ff               Required force field to use, either tripos, gaff, glycam##
                       or amber## (where ## designates the year of the glycam
                       or amber file (the number between parm or glycam and
                       the .dat in dat/leap/parm)), oplsaa, or gmx#### (where
                       the #### represents the particular, force field's
                       designation). Tables for gmx43a1, gmx43a2, gmx43b1,
                       gmx45a3, gmx53a5, and gmx53a6 only are provided.  Of
                       these, only the gmx53a6 has received any checking.
                       The glycam option is also very experimental and may not
                       give a correct result. For gromacs force fields, the
                       resultant topology is a united atoms topology.  All other
                       cases result in all atoms topologies.

     -r                Optional residue name to be used for atoms in the output.

     -last             Optional last number of a topology to which this topology
                       is to be appended.  The first number for an atom in this
                       topology will be one more than the number given.  A
                       -resnum is required along with this option.

     -resnum           Optional last residue number of a topology to which this
                       topology is to be appended.  The first number for the
                       first residue in this topology will be one more than the
                       number given.  A -last is required along with this option.

     -renumall         Optional flag to renumber topology as well as coordinates
                       to permit appending topology elements to a main topology
                       to permit addition of inter-molecular restraints. This
                       option requires that both -last and -resnum be given.
                       Usually a -r should be given as well.

     -meas             Optional flag to use measured bond lengths, angles and
                       torsions if set.  Otherwise, use the values found in the
                       selected forcefield.  Note that neither for oplsaa, nor
                       for standard Ryckaert-Bellemans dihedrals is recalculation
                       of the torsion parameters to the measured values done.

     -gro              Only write out gromacs coordinates (.gro) file.  The -r,
                       -last, -resnum, -renumall, and -rename options will act
                       upon the coordinate file appropriately.  The -meas option
                       is inappliacble to gromacs coordinate files and will be
                       ignored if specified.

     -rename           Optional flag to rename atoms for the output.

     -purge [level]    Optional flag to decrease the number of dihedral angles
                       included in the topology.  A single digit between 0 and
                       5 can be specified to indicate the extent of the decrease
                       with 0 as unpurged and 5 as most stringent purge.  Default
                       value is 1.  The recommended settings are:
                            0 for OPLS-AA, amber, gaff, and glycam force fields;
                            1 for gromacs (gmx) force fields.

     -move             Optional flag to translate the molecule so that it is
                       centered on its center of mass.  This option conflicts
                       with -renumall because it is assumed that -renumall is
                       done for a molecule aligned with another to which this
                       topology is to be added.

     -charge           Optional flag to assign charges from the force field
                       atom types to replace the .mol2 charges.  Currently
                       only works with oplsaaa.  Note that the charges assigned
                       will not be adjusted to give an integral charge to the
                       molecule.

     -h                Prints program version information and the usage message.

     -H                Same as -h.

     -v                Prints program version information.

     -V                Same as -v.


METHOD OF OPERATION:
The program first reads the coordinates and bonds information from a tripos
format file (*.mol2, where the * denotes the common name) with the specified
common name and analyzes the structure given to determine all the bonds, rings,
angles, torsions, and impropers in the structure.  It then determines a main
chain for the structure and revises the atom naming accordingly if requested.
Once this is complete it reads either the appropriate ATOMTYPE_*.DEF from the
specified directory, or Tripos force field files, TAFF*, from the specified
tripos directory.  When the program uses ATOMTYPE*.DEF files, it determines
atom types for the atoms in the structure.  With that information, it can then
read the .dat of an appropriate force field parameters file from the specified
directory to set forces for the structure's bonds, angles, torsions, and
impropers, and to set the van der Waals parameters related to pair interactions
within the structure.  When the program reads Tripos force field files, only the
structure's force field parameters need be set because Tripos bond and atom
types are already known from the input. The program then writes a gromacs coordinate
file (*.gro, where the * denotes the common name), followed by gromacs topology
files (*.top, ff*nb.itp, posre*.itp, and ff*.itp, where the * denotes the common
name). A log of the process is also written (*.log, where the * denotes the common
name) along with a tripos format file with the revised naming (*MOL.mol2, where
the * denotes the common name). Paramters for which values are not found in the
selected force field are marked with asterisks for their values in the log file
and the fields are left blank in the gromacs topology files to permit manual
editing to enter appropriate values.  In general, post topology generation
editing will be necessary.  In particular, some of the force field parameter
files in the antechamber data are missing some combinations of atom types for
dihedral angles and some van der Waals parameters for some ions.  Correct van der
Waals parameters for most ions are provided in ffusernb.itp and ions_gaff.itp
provided in the water_and_ions subdirectory.

Restraints in a tripos format file are converted to appropriate gromacs restraints
for the topology generated and can be activated with the following defines in the
gromp.mdp file:
     -DDISTREST        simple distance restraints
     -DDISTRANGE       distance range restraints
     -DDIHEDREST       dihedral restraiints
Simple Tripos distance restraints are modeled as harmonic potential type 6.  Tripos
range restraints are modeled as gromacs distance restraints with the second upper
bound set to 1 nm further than the first upper bound, index values set to consecutive
numbers to avoid simultaneous treatment of restraints, the second type column set to
1, and the multiplying factor set to 1.0.

If residue renaming is specified, all residue names are replaced by the single
residue name specified on the command line and the residue number will be set
to 1 for all atoms in the molecule.  If residue renaming is not specified,
the original residue names will be edited to remove tailing numbers.  Also,
because Sybyl does not necessarily keep all of the atoms in a residue together
so that residue numbers will be continuously increasing, the order of atoms
is be adjusted so that the residue number in the topology starts at 1 and
increases monotonically.  This process also eliminates zero mass particles that
Sybyl inserts as lone pairs, and assures that a residue number of zero is not
used.

See the file PROGRAM_TESTING_DONE for how this program has been tested to this
point.


Method Used for Gromacs Force Fields
Gromacs force fields have a limited set of atom types that are used in multiple
chemical environments. This works quite well when one can use a residue types
file to generate a topology. This presents problems when the goal is to determine
the topology of an arbitrary molecule for which residue type information does
not exist. Therefore, conversion to gromacs force fields requires double
determination of atom types. A first determination of atom types is performed
with a rich set of atom types derived from a combination of the gaff and amber
atom type definitions. This assignment is used to determine force constants,
bond lengths, angles, dihedral angles, and improper angles. A second
determination of atom types is performed to assign the true gromacs atom types
to match the atom types in the appropriate ffG####nb.itp from the gromacs
distribution.  A major problem of this double conversion is that for something
found in a residue type file, topolbuild does not always select the exact same
parameters as are found in the residue type file’s parameterization.  However,
tests so far give the same parameterization better than 95% of the time.


Method Used for OPLS-AA Force Fields
With the exception of van der Waals parameters and default atom charges, oplsaa
is designed such that force field parameters are associated with atom types
from a modified and amplified version of the Kollman atom types.  However,
OPLS-AA also uses a much richer set of atom types tto determine van der Waals
parameters and default atom charges.  Therefore, generation of an OPLS-AA
topology requires double determination of atom types.  A first determination
of atom types is performed with the expanded version of the Kollman atom types.
This assignment is used to assign force constants, bond lengths, angles,
dihedral angles, and improper angles.  A second determination of atom types
is performed to assign the final oplsaa atom types to match the atom types in
ffoplsaanb.itp from the gromacs distribution.  A major problem of this double
conversion is topolbuild does not always select the exact same atom type for a
residue listed in ffoplsaa.rtp.  The following table gives the differences
discovered:
     residue    atom    rtp value      topolbuild            usage
       HISA     ND1     opls_503       opls_557     imidazole N1
       HISA     HD1     opls_504       opls_562     imidazole H1
       HISA     HE1     opls_146       opls_563     imidazole H2
       HISB     HD2     opls_146       opls_565     imidazole H5
       HISB     HE1     opls_146       opls_563     imidazole H2
       HISB     NE2     opls_503       opls_557     imidazole N1
       HISB     HE2     opls_504       opls_562     imidazole H1
       TRP      HD1     opls_146       opls_597     indole H2
       TRP      CE3     opls_145       opls_590     indole C4
       TRP      HE3     opls_146       opls_599     indole H4
       TRP      CZ2     opls_145       opls_593     indole C7
       TRP      HZ2     opls_146       opls_602     indole H7
       TRP      CZ3     opls_145       opls_591     indole C5
       TRP      HZ3     opls_146       opls_600     indole H5
       TRP      CH2     opls_145       opls_592     indole C6
       TRP      HH2     opls_146       opls_601     indole H6

Because of the manner in which OPLS-AA constants were developed originally,
in all of these cases, the same bond, angle, dihedral, and improper constants
are assigned.  The van der Waals parameters are identical between the original
rtp atom type and the topolbuild choice as well.  Because topolbuild uses
charges from the mol2 file rather than the OPLS-AA atom type charges, differences
in atom type charges are moot.  Thus, the effects of these changes are minimal.


EXAMPLE:
     topolbuild -dir /local/amber10 -ff gaff -n myown -r RESID -meas -rename -move

Analyzes the file myown.mol2 and generates the log file myown.log, the revised
tripos format file myown MOL.mol2, a gromacs coordinates file myown.gro, and the
gromacs topology files myown.top, ffmyownnb.itp, posremyown.itp, and ffmyown.itp
with bond lengths, angles, and torsions determined from the coordinates, and atoms
renamed.  The molecule is centered on its center of mass.


USE IN DOCKING:
Ligand topologies generated by topolbuild can have atoms and residues numbered
so as to be sequential to the atoms in another topology.  Manual editing is
then necessary to include the atoms of the ligand topology in the files for
the other topology.  In doing this, remember to remove duplicate non-bonded
interaction parameters, remove the include duplications, properly edit
the molecules line, change the ligand residue numbers, and do not mix force
fields.  The last might require that both the ligand and its receptor be
processed by topolbuild.  Otherwise, processing for docking of a ligand
is similar to the processing used for ligand topologies derived from
other sources.


REVISIONS:
Sept. 1, 2009:  Version 1.3:  Corrected errors in setting some categories
                of rings
                Corrected Makefile to be compatible with more variants of
                make
                Added OPLS-AA force field tables
                Increased initial estimate of maximum number of rings for
                ring detection
                Added lines to log, topology, and include files to give
                version and command line of invocation
                Made internal rearrangements to ease addition of and handling
                of other force fields
                Set default prune of excess dihedral angles equivalent to
                -purge 1
                Added option to translate molecule coordinates to center of
                mass when renumbering is not requested
                Added elements to internal atomic masses / numbers table
                Revised method of correction of Tripos out of plane pyramid
                height improper force constants to cosine function force
                constants
                Altered Tripos force field setting of dihedral phases
                Changed error messages in mol2 file reading to give clearer
                statement of problems
                Changed error for not ending the bonds section of the mol2
                file with @<TRIPOS> from fatal error to warning
                Added atom type default entries to the amber atom type
                definition tables
                Corrected errors in gromacs topology defaults line settings.
                Changed preferences in dihedral angles purge to choose the
                maximum number of heavy atoms possible for dihedrals selected
                to be retained in the topology output.
                Added option to set charges from the force field atom type
                charges.


PROBLEMS NOTED, CORRECTIONS, AND POSSIBLE BUGS:
 1. Tripos force field has an out of plane bending parameter with a force
    constant in units of kcal/(mol Angstrom^2) instead of angle based
    impropers.  The program attempts to properly scale and set improper
    force values, but improper settings from Tripos force field may require 
    post-conversion editing.
 2. Ion topology and water topology files for amber and Tripos force fields
    must be separately generated.  Examples are provided in the water_and_ions
    directory.
 3. With large molecules such as proteins, polynucleotides, and peptides,
    round-off errors in the summation of charges result in small deviations
    in total system charge from integral charge of the order of 0.004.  Of
    course, such systems would be better processed with pdb2gmx when possible.
    Small ligands do not appear to suffer from this problem.
 4. The various published force field parameter files do not include parameters
    for dihedrals, impropers, and angles for every combination of atom types
    that might occur. The user is strongly advised to examine the topology
    produced for omissions. To assist in this, the log file generated by the
    program marks missing parameters with strings of asterisks, "****", and
    the topology file produced leaves spaces where the missing parameters
    ought to appear.
 5. Charge group number assignment will not be optimal in all cases and may
    need to be edited.
 6. Parameter assignments are based on the structure as given in the mol2
    input file and on the force field chosen.  The resulting assignments
    will need to be evaluated to be sure that it matches expectations as
    well as experimental data or quantum chemical modeling results correctly.
 7. Gromacs ignores case for atom types, but gaff is case sensitive.  To
    preclude confusion of the gaff nitrogen atom type na, with sodium, ions
    have their charge appended to their name.
 8. Some problems with atom types and bond types have been observed with mol2
    files generated by non-Tripos programs as follows:
    a. Charges generated are larger than charges generated by Sybyl for the
       same structure by factors varying between 2 and 20.  An adjustment
       factor has been added to the code to correct for some of this, but
       there will be cases where the excessive charge causes incorrect atom
       type determinations.
    b. Sybyl uses atom type "C.cat" for CZ of Arginine, but other programs
       choose atom type "C.2" instead.  Manual editing of the input is
       necessary to correct this.
    c. Sybyl uses bond type "ar" for all bonds to CZ of Arginine, but other
       programs use bond type "2".  This causes the atom types for the
       nitrogens to be incorrect.  The only way to correct this is to manually
       edit the bond types in the mol2 file.
    d. Sybyl uses bond types "1" and "2" within the Histidine ring, but other
       programs use bond type "ar" for all bonds within the Histidiine ring.
       This causes atom types within the ring to be incorrectly set.  Manual
       editing of the mol2 input file is needed to give correct input.
    e. Some programs terminate the file after the last bond in the bonds
       section instead of having a line that begins with @<TRIPOS> following
       the line with the last bond in it.
    f. Apparently, babel does not always assign complete Sybyl atom types
       on conversion of pdb to mol2.
 9. Several changes and additions to the definitions files and tables supplied
    with AmberTools were made as follows: 
    a. With the standard atom type definitions file, the CZ of arginine becomes
       amber type CM.  However, the type assignment was supposed to be CA.  To
       correct this, a new line must be added to the file, ATOMTYPE_AMBER.DEF,
       between the lines that define type CA and CD that reads:
          ATD  CA    *   6   3   *   *   *               (N3,N3,N3)    &
    b. A default entry is lacking for the two bond nitrogen atom type.  A new
       line was added to ATOMTYPE_AMBER.DEF to provide this as follows:
          ATD  N2    *   7   2   &
    c. Amber atom type definitions as supplied in AmberTools defines the iodine
       atom type as a fluorine.  The line that reads:
          ATD  F     *   53  1   &
       was changed to read:
          ATD  I     *   53  1   &
    d. Amber parameters lack parameterization for the C -CT-NT angle.  To
       correct this, a new line must be added to the files, parm**.dat, after
       the definition of the C -CT-N3 angle that reads:
          C -CT-NT    80.0      111.20    AA amino terminal residues
10. In the OPLS-AA force field, constants dervied from amber parameters were
    added for missing dihedral constants for polyphosphate chains.


ACKNOWLEDGEMENTS:
Much of this is derivative work based on study of:
   1. The routines found GROMACS version 3.3.1 written by David van der Spoel,
      Erik Lindahl, Berk Hess, and others and  copyright (c) 1991-2000,
      University of Groningen, The Netherlands.  Copyright (c) 2001-2004,
      The GROMACS development team
   2. The routines found in antechamber 1.27 particularly in:
      a. program atomtype, version 1.0, dated October, 2001 by Junmei Wang,
         Department of Pharmaceutical Chemistry, School of Pharmacy,
         University of California, San Francisco, CA  94143
      b. programs parmchk and prepgen, version 1.0, dated October, 2001 by
         Junmei Wang, Department of Pharmaceutical Chemistry, School of
         Pharmacy, University of California, San Francisco, CA  94143
      c. function file rings.c, by Junmei Wang, Department of Pharmaceutical
         Chemistry, School of Pharmacy, University of California, San
         Francisco, CA  94143
      d. program charmgen, version 2.0, dated June, 2004 by Victor E.
         Bazterra, Center for High Performance Computing, University
         of Utah, Salt Lake City, UT  84114
      e. program bondtype, version 1.0, dated October, 2001 by Junmei Wang,
         Department of Pharmaceutical Chemistry, School of Pharmacy,
         University of California, San Francisco, CA  94143
      f. function file mol2.c  by Junmei Wang, Department of Pharmaceutical
         Chemistry, School of Pharmacy, University of California, San
         Francisco, CA  94143
  3. Routines found in AMBCONV by Filip Ryjacek, 2002 as distributed at the
     GROMACS web site.
  4. The routines found in Dock 6.1 file dockmol.cpp version dated December 2006,
     Molecular Design Institute, University of California, San Francisco, CA  94143

