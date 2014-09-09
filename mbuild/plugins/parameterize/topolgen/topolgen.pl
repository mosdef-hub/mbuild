#!/usr/bin/perl

#
# TopolGen - written by Justin Lemkul
#
# Version 1.0 - 3/21/2009
#
# Development To-do's
#	1. Assign atom types (currently writes "opls_XXX" for all atoms as a placeholder)
#	2. Assign charges (currently writes 0.000 for all charges as a placeholder)
#		* real charge groups might be nice
#	3. Add dihedral check for amide groups (i.e., GLN sidechain)
#		* no dihedrals involving amides, only impropers = grompp complains 
#

use strict;

# Date and version info for printing to topology
chomp(my $date = `date`);
my $v = "1.0 (3/21/2009)";

# open the input (provided by command line)
unless (@ARGV) {
	die "Usage: perl $0 <input.pdb>\n";
}

my $input = $ARGV[0];

# Define some pre-determined function types required by OPLS
my $funct_bond = 1;	# Harmonic bond potential
my $funct_angle = 1;	# Harmonic angle
my $funct_dihedral = 3;	# Ryckaert-Bellemans dihedral
my $funct_pairs = 1;	# Lennard-Jones pair potential
my $funct_improper = 1;	# "Proper" dihedral representing the improper

# Generate the [ atoms ] section and header information
open(ATOMS_IN, $input) or die "$!\n";
my @input_file = <ATOMS_IN>;
close(ATOMS_IN);

# strip out only the HETATM records, write to new array
my @atoms;

foreach $_ (@input_file) {
        #if ($_ =~ /^HETATM/) {
        if ($_ =~ /^ATOM/) {
                push(@atoms, $_);
        }
}

#
# Fields in the PDB file:
# HETATM atom# atom_name residue res# x y z beta occ type
#
# Can throw away HETATM, x/y/z coords, beta, and occupancy
#

my @line;
my @lines_edit;

foreach $_ (@atoms) {
	chomp(@line = split(" ", $_));
	shift(@line);					# remove HETATM
	my @coords_b_occ = splice(@line, 4, 5); 	# remove x/y/z, beta, occupancy
							# store in case we need it
	my $string = join(" ", @line);
	# print "$string\n";
	push(@lines_edit, $string);
}

#
# So now, in each element of @line, there is a string that contains
# atom# atom_name residue res# type
#

chomp(my @temp_info = split(" ", $lines_edit[0]));
my $res_name = $temp_info[2];

# Open output file
open(ATOMS_OUT, ">>section_atoms");
print ATOMS_OUT ";\n";
print ATOMS_OUT "; \tOPLS-AA topology, built by TopolGen version $v\n";
print ATOMS_OUT "; \tScript written by: Justin Lemkul (jalemkul\@vt.edu)\n";
print ATOMS_OUT "; \tThis is your molecule's topology\n";
print ATOMS_OUT "; \tCheck it carefully for any errors. It is not necessarily perfect!\n";
print ATOMS_OUT ";\n";
print ATOMS_OUT "; \tTopology written on $date\n";
print ATOMS_OUT ";\n";
print ATOMS_OUT "; Include force field\n";
print ATOMS_OUT "#include \"ffoplsaa.itp\"\n";
print ATOMS_OUT "\n";
print ATOMS_OUT "[ moleculetype ]\n";
printf ATOMS_OUT "; Name%18s\n", "nrexcl";
printf ATOMS_OUT "%-6s%15d\n", $res_name, 3;
print ATOMS_OUT "\n";
print ATOMS_OUT "[ atoms ]\n";
printf ATOMS_OUT ";%5s%11s%7s%8s%6s%7s%11s%11s%7s%11s%11s\n", "nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass", "typeB", "chargeB", "massB";
close(ATOMS_OUT);

foreach $_ (@lines_edit) {
	chomp(my @info = split(" ", $_));

	my $atom_num = $info[0];
	my $atom_name = $info[1];
	my $residue_name = $info[2];
	my $res_num = $info[3];
	my $atom_type = $info[4];
	my $atom_mass;
	my $atom_charge;

	if ($atom_type eq "C") {
		$atom_mass = 12.011;
	} elsif ($atom_type eq "H") {
		$atom_mass = 1.008;
	} elsif ($atom_type eq "O") {
		$atom_mass = 15.9994;
	} elsif ($atom_type eq "N") {
		$atom_mass = 14.0067;
	} elsif ($atom_type eq "S") {
		$atom_mass = 32.06;
	} elsif ($atom_type eq "P") {
		$atom_mass = 30.97376;
	} elsif ($atom_type eq "F") {
		$atom_mass = 18.9984;
	} elsif ($atom_type eq "SI") {
		$atom_mass = 28.08;
	} else {
		print "Unknown atom found.  Please check the output!\n";
		$atom_mass = 0.000;
	}

	# Print to output
	open(ATOMS_OUT, ">>section_atoms");
	printf ATOMS_OUT "%6d%11s%7d%7s%7s%7d%11.3f%11.5f%7s%11s%11s\n", $atom_num, "opls_XXX", $res_num, $residue_name, $atom_name, 1, $atom_charge, $atom_mass;
	close(ATOMS_OUT);

	# Print to dummy output, to be used later in improper determination
	# open(ATOMS_OUT_RAW, ">>section_atoms_raw");
        # printf ATOMS_OUT_RAW "%6d%11s%7d%7s%7s%7d%11.3f%11.5f%7s%11s%11s\n", $atom_num, "opls_XXX", $res_num, $residue_name, $atom_name, 1, $atom_charge, $atom_mass;
        # close(ATOMS_OUT_RAW);
	
	open(ATOMS_OUT_RAW, ">>section_atoms_raw");
	print ATOMS_OUT_RAW "@info\n";
	close(ATOMS_OUT_RAW);	

}

# Generate the [ bonds ] section

open(IN, $input) or die "$!\n";
my @input_file = <IN>;
close(IN);

# strip out only the CONECT records, write to new array
my @conect;

foreach $_ (@input_file) {
	if ($_ =~ /^CONECT/) {
		push(@conect, $_);
	}
}

#
# Within CONECT section, "duplicate" entries exist - C-H bonds, then H-C bonds
#
# Steps:
#  1. Push each line into a new array, split by whitespace
#  2. Delete "CONECT" from each line (shift)
#  3. Loop over length of line to assign bonds: 1-2, 1-3, 1-4, 1-5
#      If $_[0] > $_[n], ignore element n 
#  4. Write to output
#

# Define new array to hold lines
my @line_edit;

# Open filehandle for bonds output
# Write the heading information for the [ bonds ] section
open(BONDS_OUT, ">>section_bonds");
printf BONDS_OUT "[ bonds ]\n";
printf BONDS_OUT ";%4s%6s%6s\n", "ai", "aj", "funct";
close(BONDS_OUT);

foreach $_ (@conect) {
	chomp(@line_edit = split(" ", $_));
	
	# Remove CONECT string
	shift(@line_edit);
	
	# Determine length of array (how long to loop)
	# At this point, the @line_edit array should hold one CONECT entry
	# Each atom number is an element in the array
	my $length = scalar(@line_edit);

	# Need to parse thru each line and satisfy these criteria
	# First time thru: $line_edit[0] = 1, $line_edit[1] = 2, etc.
	# Need to write out bonds based on first atom
	# If $line_edit[0] > $line_edit[$i], do nothing; check next element 

	my $string;

	foreach $_ (@line_edit) {
		if ($line_edit[0] >= $_) {
			# Do nothing! Includes self check - not pretty, but effective
		} else {
			open(BONDS_OUT, ">>section_bonds");
			printf BONDS_OUT "%5d%6d%6d\n", $line_edit[0], $_, $funct_bond;
			close(BONDS_OUT);
		}
	}

}

# Now, parse through bond entries and determine angles
#
# Defining angles based on input text file
#
# 	1. Open "section_bonds" as an array
# 	2. shift off headers ([ bonds ], ; ai aj...)
#	3. Split each line by space to @new_line array
#	4. pop off function type
#	5. $atom_j = $new_line[0], $atom_i = $new_line[1]
#	6. Need to find $atom_k for given angle, will be $new_line[1] of subsequent (not necessarily
#		consecutive) line
#		* Need to devise some sort of counter to run thru array elements...
#		* For example, in propane, atoms 1-5-8 are bonded (C-C-C)
#		* If bond(1-5) exists, search for 5 as $atom_j in following lines
#			- if $new_line[0] == $atom_j, then $atom_k = $new_line[1]
#	7. Write angles: $atom_i, $atom_j, $atom_k (i-j-k) with j being the vertex of the angle
#

open(IN, "<section_bonds") or die "$!\n";
my @angle_input_file = <IN>;
close(IN);

# shift off headers
shift(@angle_input_file);
shift(@angle_input_file);
chomp(@angle_input_file);

# Write headers for output
open(ANGLES_OUT, ">>section_angles");
printf ANGLES_OUT "[ angles ]\n";
printf ANGLES_OUT ";%4s%6s%6s%6s\n", "ai", "aj", "ak", "funct";
close(ANGLES_OUT);

my $length = scalar(@angle_input_file);

for (my $i=0; $i<$length; $i++) {
	chomp(my @line = split(" ", $angle_input_file[$i]));
	# remove function value
	pop(@line);
	
	# atom numbers need to be held in variables so they can be written to the final output
	# angle is defined as i-j-k, where j is the vertex of the angle
	my $atom_i = $line[0];
	my $atom_j = $line[1];
	my $atom_k;

	# Find matches on remaining lines
	# Initialize counter
	my $j = $i+1;

	# Start loop
	for ($j; $j<=$length; $j++) {

		# if $next_line[0] == $atom_j, then $atom_k = $next_line[1]
		chomp(my @next_line = split(" ", $angle_input_file[$j]));
		pop(@next_line);

		if ($next_line[0] == $atom_j) {
			$atom_k = $next_line[1];

			# Write output	
			open(ANGLES_OUT, ">>section_angles");
                        printf ANGLES_OUT "%5d%6d%6d%6d\n", $atom_i, $atom_j, $atom_k, $funct_angle; 
                        close(ANGLES_OUT);

		} elsif ($next_line[0] == $atom_i) {
			$atom_k = $next_line[1];

			# Write output
			open(ANGLES_OUT, ">>section_angles");
			printf ANGLES_OUT "%5d%6d%6d%6d\n", $atom_j, $atom_i, $atom_k, $funct_angle;
			close(ANGLES_OUT);
		}
	}

}

# Defining dihedrals based on input text file
#
# 	1. Open "section_angles" as an array
# 	2. shift off headers ([ angles ], ; ai aj...)
#	3. Split each line by space to @new_line array
#	4. pop off function type
#	5. $atom_i = $next_line[0], $atom_j = $next_line[1], $atom_k = $next_line[2]
#		* For a given angle, bonds exist between $atom_i-$atom_j, $atom_j-$atom_k
#	6. Need to find $atom_l for given angle, provided that a bond exists between $atom_k & $atom_l 
# 		* Need to parse thru "section_bonds" and test for presence of 
#			$atom_k as $line[0] => $atom_l = $line[1]
#		* For example, in propane, atoms 2-1-5 is an angle (H-C-C)
#		* If bond(5-8) exists, search for 5 as $atom_k  as first element in bonds line
#			- if $line[0] == $atom_k, then $atom_l = $line[1]
#	7. Write dihedrals: $atom_i, $atom_j, $atom_k, $atom_l (i-j-k-l)
#

# Open file containing bonds
open(IN2, "<section_bonds") or die "$!\n";
my @dihedral_input_file_bonds = <IN2>;
close(IN2);

# shift off headers
shift(@dihedral_input_file_bonds);
shift(@dihedral_input_file_bonds);
chomp(@dihedral_input_file_bonds);

# Open file containing angles
open(IN3, "<section_angles") or die "$!\n";
my @dihedral_input_file_angles = <IN3>;
close(IN3);

# shift off headers
shift(@dihedral_input_file_angles);
shift(@dihedral_input_file_angles);
chomp(@dihedral_input_file_angles);

# Write headers for dihedral output
open(DIHEDRALS_OUT, ">>section_dihedrals");
printf DIHEDRALS_OUT "[ dihedrals ]\n";
printf DIHEDRALS_OUT ";%4s%6s%6s%6s%6s\n", "ai", "aj", "ak", "al", "funct";
close(DIHEDRALS_OUT);

# Write headers for pairs output
open(PAIRS_OUT, ">>section_pairs");
printf PAIRS_OUT "[ pairs ]\n";
printf PAIRS_OUT ";%4s%6s%6s\n", "ai", "aj", "funct";
close(PAIRS_OUT);

my $length_angles = scalar(@dihedral_input_file_angles);

for (my $i=0; $i<$length_angles; $i++) {
	chomp(my @line = split(" ", $dihedral_input_file_angles[$i]));
	# remove function value
	pop(@line);
	
	# atom numbers need to be held in variables so they can be written to the final output
	# dihedral is defined as i-j-k-l, where angle i-j-k exists, as well as bond k-l
	my $atom_i;
	my $atom_j;
	my $atom_k;
	my $atom_l;

	# define the values for the atoms, taken from angle file
	$atom_i = $line[0];
	$atom_j = $line[1];
	$atom_k = $line[2];

	# Search for existing bond starting with $atom_k
	my $length_bonds = scalar(@dihedral_input_file_bonds);
	for (my $j=0; $j<$length_bonds; $j++) {
		chomp(my @line2 = split(" ", $dihedral_input_file_bonds[$j]));
		pop(@line2);

		if ($line2[0] == $atom_k) {
			$atom_l = $line2[1];

			# Write dihedral output	
			open(DIHEDRALS_OUT, ">>section_dihedrals");
                        printf DIHEDRALS_OUT "%5d%6d%6d%6d%6d\n", $atom_i, $atom_j, $atom_k, $atom_l, $funct_dihedral; 
                        close(DIHEDRALS_OUT);
		
			# Write pairs output
			open(PAIRS_OUT, ">>section_pairs");
			printf PAIRS_OUT "%5d%6d%6d\n", $atom_i, $atom_l, $funct_pairs;
			close(PAIRS_OUT);
	
		}
	}
}

# Defining impropers based on input text file
#
# 	1. Open "section_angles" and "section_bonds" as arrays
# 	2. shift off headers ([ angles ], ; ai aj...)
#	3. Split each line by space to @new_line array
#	4. pop off function type
#	5. For each atom ($atom_x), count how many times it appears in section_bonds
#		* if $line[0] == $atom_x || $line[1] == $atom_x, $counter++
#		* if $counter == 3, define an improper around that atom
#	6. If $atom_x is "i" for a given improper, then $atom_j, $atom_k will come from section_angles 
# 		* $atom_x = $line[1] => $atom_i; $atom_j = $line[0], $atom_k = $line[2]  
#		* Read thru bonds, if $line[1] != $atom_j && $line[1] != $atom_k, then $atom_l = $line[1] 
#	7. Write impropers: $atom_i, $atom_j, $atom_k, $atom_l (i-j-k-l), function type
#		* Also need to write out "improper_A_B_C_D" on each line
#			- impropers will have to be based on atom types
#			- THIS FEATURE SHOULD NOT BE CONSIDERED 100% RELIABLE 

# Open file containing bonds
open(IN4, "<section_bonds") or die "$!\n";
my @improper_input_file_bonds = <IN4>;
close(IN4);

# shift off headers
shift(@improper_input_file_bonds);
shift(@improper_input_file_bonds);
chomp(@improper_input_file_bonds);

# Open file containing angles
open(IN5, "<section_angles") or die "$!\n";
my @improper_input_file_angles = <IN5>;
close(IN5);

# shift off headers
shift(@improper_input_file_angles);
shift(@improper_input_file_angles);
chomp(@improper_input_file_angles);

# Write headers for improper output
open(IMPROPERS_OUT, ">>section_impropers");
printf IMPROPERS_OUT "[ dihedrals ]\n";
printf IMPROPERS_OUT ";%4s%6s%6s%6s%6s\n", "aj", "ak", "ai", "al", "funct";
close(IMPROPERS_OUT);

my $length_angles = scalar(@improper_input_file_angles);

# initialize counter to move through array
for (my $i=0; $i<$length_angles; $i++) {
	chomp(my @line = split(" ", $improper_input_file_angles[$i]));
	# remove function value
	pop(@line);
	
	# atom numbers need to be held in variables so they can be written to the final output
	# improper is defined as i-j-k-l, angle is i-j-k, but here I am re-defining the numbering
	# so that $atom_i is always "i"; angles will be interpreted as j-i-k, just for ease of writing
	# the improper output 
	my $atom_i;
	my $atom_j;
	my $atom_k;
	my $atom_l;

	# define the values for the atoms, taken from angle file
	$atom_j = $line[0];
	$atom_i = $line[1];
	$atom_k = $line[2];

	# Initialize a counter that determines how many bonds $atom_i participates in
	my $counter = 0;

	# Search for bonds containing $atom_i
	my $length_bonds = scalar(@improper_input_file_bonds);
	for (my $j=0; $j<$length_bonds; $j++) {
		chomp(my @line2 = split(" ", $improper_input_file_bonds[$j]));
		pop(@line2);

		if ($line2[0] == $atom_i || $line2[1] == $atom_i) {
			$counter++;
			
			if ($line2[0] == $atom_i && ($line2[1] != $atom_j && $line2[1] != $atom_k)) {
				$atom_l = $line2[1];
			} elsif ($line2[1] == $atom_i && ($line2[0] != $atom_j && $line2[0] != $atom_k)) {
				$atom_l = $line2[0];
			}
		}

	}

	if ($counter == 3 && $atom_l != 0) {
	
		# Write improper output	
		open(IMPROPERS_OUT_RAW, ">>section_impropers_raw");
		printf IMPROPERS_OUT_RAW "%5d%6d%6d%6d%6d\n", $atom_i, $atom_j, $atom_k, $atom_l, $funct_improper; 
		close(IMPROPERS_OUT_RAW);
	}
}

# Splice out repeats

open(IMPROPERS_FIX, "<section_impropers_raw");
my @impropers_fix = <IMPROPERS_FIX>;
close(IMPROPERS_FIX);

my $impropers_fix_length = scalar(@impropers_fix);

for (my $i=0; $i<$impropers_fix_length; $i++) {
	# Loop thru array, if next line (or any subsequent line) contains $atom_i as its first element,
	# splice that element out of @impropers_fix
	chomp(my @line3 = split(" ", $impropers_fix[$i]));
	pop(@line3);

	my $atom_i = $line3[0];
	my $atom_j = $line3[1];
	my $atom_k = $line3[2];
	my $atom_l = $line3[3];

	# Find additional matches in subsequent lines
	# Initialize a new counter
	my $j = $i+1;

	for ($j; $j<$impropers_fix_length; $j++) {
		chomp(my @next_line = split(" ", $impropers_fix[$j]));
		pop(@next_line);

		if ($next_line[0] == $atom_i) {
			splice(@impropers_fix, $j, 1);
			# Re-set array length and counter since array has been shortened
			$impropers_fix_length--;
			$j--;
		}
		
	}

	# Decision structure to decide on improper types, based on global atom name
	# Loop thru "section_atoms", parse out atom type from field
	my $atom_type_i;
	my $atom_type_k;
	my $atom_type_l;

	open(ATOM_SEARCH, "<section_atoms_raw");
	my @atom_search = <ATOM_SEARCH>;
	close(ATOM_SEARCH);
	
	foreach $_ (@atom_search) {
		chomp(my @atom_search_parse = split(" ", $_));

		if ($atom_i == $atom_search_parse[0]) {
			$atom_type_i = $atom_search_parse[4]; 
		}

		if ($atom_l == $atom_search_parse[0]) {
			$atom_type_l = $atom_search_parse[4];
		}

		if ($atom_k == $atom_search_parse[0]) {
			$atom_type_k = $atom_search_parse[4];
		}	

	}

	# Declare output string for improper
	my $improper_string;

	if ($atom_i != 0) {

		if ($atom_type_i eq "C" && $atom_type_l ne "O" && $atom_type_k ne "O") {
			$improper_string = "improper_Z_CA_X_Y";
		} elsif ($atom_type_i eq "N") {
			$improper_string = "improper_Z_N_X_Y";
		} elsif ($atom_type_i eq "C" && ($atom_type_l eq "O" || $atom_type_k eq "O")) {
			$improper_string = "improper_O_C_X_Y";
			# weird switching necessary to get improper in right order
			($atom_k, $atom_l) = ($atom_l, $atom_k);
		} else {
			print "Unknown improper found. Check your topology!\n";
		}
		
		# Write final output
		open(IMPROPERS_OUT, ">>section_impropers");
		# This is the original print statement
		# printf IMPROPERS_OUT "%5d%6d%6d%6d%6d\n", $atom_i, $atom_j, $atom_k, $atom_l, $funct_improper;
		# apparently, the OPLS improper format requires something like: j k i l
		printf IMPROPERS_OUT "%5d%6d%6d%6d%6d%4s%-16s\n", $atom_j, $atom_k, $atom_i, $atom_l, $funct_improper, "    ", $improper_string;         
		close(IMPROPERS_OUT);
	}

}

unlink "section_impropers_raw";

# Clean up the output to write some newlines between the sections
# This makes the concatenated output a bit prettier

open(ATOMS_OUT, ">>section_atoms");
printf ATOMS_OUT "\n";
close(ATOMS_OUT);

open(BONDS_OUT, ">>section_bonds");
printf BONDS_OUT "\n";
close(BONDS_OUT);

open(ANGLES_OUT, ">>section_angles");
printf ANGLES_OUT "\n";
close(ANGLES_OUT);

open(DIHEDRALS_OUT, ">>section_dihedrals");
printf DIHEDRALS_OUT "\n";
close(DIHEDRALS_OUT);

open(PAIRS_OUT, ">>section_pairs");
printf PAIRS_OUT "\n";
close(PAIRS_OUT);

open(IMPROPERS_OUT, ">>section_impropers");
printf IMPROPERS_OUT "\n";
close(IMPROPERS_OUT);

# Concatenate the final output into one single file

# Open output filehandle
open (TOPOLOGY_OUT, ">>ffoplsaa_TopolGen_${res_name}.top") or die "$!\n";

# Open input files to be written to output
# [ atoms ] section
open(FIRST_IN, "<section_atoms") or die "$!\n";

while (my $line_out = <FIRST_IN>) {
	print TOPOLOGY_OUT $line_out;
}
close(FIRST_IN);

# [ bonds ] section
open(SECOND_IN, "<section_bonds") or die "$!\n";

while (my $line_out = <SECOND_IN>) {
	print TOPOLOGY_OUT $line_out;
}
close(SECOND_IN);

# [ pairs ] section
open(THIRD_IN, "<section_pairs") or die "$!\n";

while (my $line_out = <THIRD_IN>) {
	print TOPOLOGY_OUT $line_out;
}
close(THIRD_OUT);

# [ angles ] section
open(FOURTH_IN, "<section_angles") or die "$!\n";

while (my $line_out = <FOURTH_IN>) {
	print TOPOLOGY_OUT $line_out;
}
close(FOURTH_OUT);

# proper [ dihedrals ] section
open(FIFTH_IN, "<section_dihedrals") or die "$!\n";

while (my $line_out = <FIFTH_IN>) {
	print TOPOLOGY_OUT $line_out;
}
close(FIFTH_IN);

# improper [ dihedrals ] section
# only write it out if it contains relevant lines
open (TEST, "<section_impropers") or die "$!\n";

chomp(my @test_array = <TEST>);

my $test_length = scalar(@test_array);

if ($test_length > 3) {
	open(SIXTH_IN, "<section_impropers") or die "$!\n";

	while (my $line_out = <SIXTH_IN>) {
		print TOPOLOGY_OUT $line_out;
	}
	close(SIXTH_IN);
}

print TOPOLOGY_OUT "[ system ]\n";
print TOPOLOGY_OUT "; Name\n";
print TOPOLOGY_OUT "$res_name, generated by TopolGen\n";
print TOPOLOGY_OUT "\n";
print TOPOLOGY_OUT "[ molecules ]\n";
printf TOPOLOGY_OUT ";%9s%13s\n", "Compound", "\#mols";
printf TOPOLOGY_OUT "%-6s%15d\n", $res_name, 1;

close(TOPOLOGY_OUT);

# Clean up intermediate files
unlink "section_atoms";
unlink "section_atoms_raw";
unlink "section_bonds";
unlink "section_angles";
unlink "section_pairs";
unlink "section_dihedrals";
unlink "section_impropers";
