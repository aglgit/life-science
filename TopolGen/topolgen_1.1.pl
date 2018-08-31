#!/usr/bin/perl

use strict;

#
# TopolGen - written by Justin Lemkul
#
# Version 1.1_dev - 10/14/2009
#
# Updates since version 1.0
#	1. Made command line options into a hash, requiring flags
#	2. Added printing of atom labels to bonded parameters sections
#	3. Added renumbering option
#	4. Added compatibility with chain identifiers and notes about atom type field
#	5. Added primitive atom typing and charge assignment based on neighboring atoms and connectivity
#	6. Added charge group assignment for each new C or N atom type found
#	7. Added dihedral check - proper dihedrals no longer erroneously generated for amides
#	8. Added ability to output either a .itp or a .top version of the topology
#

# Date and version info for printing to topology
chomp(my $date = `date`);
my $v = "1.1_dev (10/14/2009)";

# Store command line arguments in a hash
my %arguments = @ARGV;

# Parse through command line arguments
if (exists($arguments{"-h"})) {
	print "
#
# TopolGen, version $v
# Written by Justin Lemkul
#
# The program takes an all-atom .pdb file and generates a GROMACS-compatible topology.
# The only force field option is OPLS-AA at this time.  Usage:
#
# perl topolgen.pl -f input.pdb -o output.top [-type itp/top] [-renumber]
#
# Questions, comments, bug reports?  Email jalemkul[at]vt[dot]edu
#\n\n";
	exit;
}

# Define input file
my $input;

if (defined($arguments{"-f"})) {
	$input = $arguments{"-f"};
} else {
	die "Input not specified. Usage: perl $0 -f <input.pdb> -o <output.top> [-type itp/top] [-renumber]\n";
}

# Define output type
my $out_type;

if (defined($arguments{"-type"})) {
	$out_type = $arguments{"-type"};
} else {
	print "\nNo output type specified. Using default type of .top\n";
	$out_type = "top";
}

# Define output file
my $output;

if (defined($arguments{"-o"})) {
	$output = $arguments{"-o"};
} else {

	if ($out_type eq "top") {
		print "\nNo output filename specified. Using default output filename \"ffoplsaa_TopolGen.top.\"\n";
		$output = "ffoplsaa_TopolGen.top"; 
	} elsif ($out_type eq "itp") {
		print "\nNo output filename specified. Using default output filename \"ffoplsaa_TopolGen.itp\"\n";
		$output= "ffoplsaa_Topolgen.itp";
	} else {
		printf("\nUnknown output type %s\n", $arguments{"-type"});
	}
}

# Check to make sure output filename and type are consistent
my @check = split('\.', $output);

if ($check[1] ne $out_type) {
	printf("\nOutput filename and type are inconsistent: file name is of type %s, specified type is %s\n\n", $check[1], $out_type);
	exit;
}

# Define some pre-determined function types required by OPLS
my $funct_bond = 1;	# Harmonic bond potential
my $funct_angle = 1;	# Harmonic angle
my $funct_dihedral = 3;	# Ryckaert-Bellemans dihedral
my $funct_pairs = 1;	# Lennard-Jones pair potential
my $funct_improper = 1;	# "Proper" dihedral representing the improper

# Generate the [ atoms ] section and header information
open(ATOMS_IN, $input) or die "Problem while reading $input for header and atoms: $!\n";
my @input_file = <ATOMS_IN>;
close(ATOMS_IN);

# strip out only the HETATM records, write to new array
my @atoms;

foreach $_ (@input_file) {
        if ($_ =~ /^HETATM/) {
                push(@atoms, $_);
        }
}

#
# Fields in the PDB file:
# HETATM atom# atom_name residue res# x y z beta occ type
#
# Can throw away HETATM, x/y/z coords, beta, and occupancy
#
# Chain identifier may be present!
#	atom# atom_name residue chain res# type
#

my @line;
my @lines_edit;

foreach $_ (@atoms) {
	chomp(@line = split(" ", $_));
	shift(@line);					# remove HETATM

        # Check to see if chain identifier is present
        if (scalar(@line) > 10) { 	# line contains chain identifier
                my @coords_b_occ = splice(@line, 5, 5);         # remove x/y/z, beta, occupancy
                                                                # store in case we need it
        } else {			# no chain identifier
                my @coords_b_occ = splice(@line, 4, 5);
        }

	my $string = join(" ", @line);
	push(@lines_edit, $string);
}

#
# So now, in each element of @line, there is a string that contains
# atom# atom_name residue res# atom_type
#

chomp(my @temp_info = split(" ", $lines_edit[0]));
my $res_name = $temp_info[2];
my $first_atom_num = $temp_info[0];
my $renumber_factor;

# Determine if renumbering is necessary, based on command line flag
if (exists($arguments{"-renumber"})) {
	if ($first_atom_num != 1) {
		$renumber_factor = $first_atom_num - 1;
	}
}

# Open output file
open(HEADER_OUT, ">>section_header");
print HEADER_OUT ";\n";
print HEADER_OUT "; \tOPLS-AA topology, built by TopolGen version $v\n";
print HEADER_OUT "; \tScript written by: Justin Lemkul (jalemkul\@vt.edu)\n";
print HEADER_OUT "; \tThis is your molecule's topology\n";
print HEADER_OUT "; \tCheck it carefully for any errors. It is not necessarily perfect!\n";
print HEADER_OUT ";\n";
print HEADER_OUT "; \tTopology written on $date\n";
print HEADER_OUT ";\n";

if ($out_type eq "top") {
	print HEADER_OUT "; Include force field\n";
	print HEADER_OUT "#include \"ffoplsaa.itp\"\n";
}

print HEADER_OUT "\n";
print HEADER_OUT "[ moleculetype ]\n";
printf HEADER_OUT "; Name%18s\n", "nrexcl";
printf HEADER_OUT "%-6s%15d\n", $res_name, 3;
print HEADER_OUT "\n";
print HEADER_OUT "[ atoms ]\n";
printf HEADER_OUT ";%5s%11s%7s%8s%6s%7s%11s%11s%7s%11s%11s\n", "nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass", "typeB", "chargeB", "massB";
close(HEADER_OUT);

# Store atom number and name in a separate hash for later use
my %atom_number_name;

# Store atom number and type in a separate hash for later use
my %atom_number_type;

foreach $_ (@lines_edit) {
	chomp(my @info = split(" ", $_));

	my $atom_num;

	# Assign numbering based on -renumber
	if (exists($arguments{"-renumber"}) && ($first_atom_num != 1)) {
		$atom_num = $info[0] - $renumber_factor;
		splice(@info, 0, 1, $atom_num);
	} else {
		$atom_num = $info[0];
	}
	
	my $atom_name = $info[1];
	my $residue_name = $info[2];
	# my $atom_type = $info[4];

	# Unknowns at this point
	my $atom_type;	
	my $res_num;
	my $atom_type;

	# Conditional assignment in case chain identifier present
	if (scalar(@info) == 6) {		# chain identifier present
		$res_num = $info[4];
		$atom_type = $info[5];
	} else {				# chain identifier not present
		$res_num = $info[3];
		$atom_type = $info[4];
	}

	# add name and number to the hash
	$atom_number_name{$atom_num} = $atom_name;
	$atom_number_type{$atom_num} = $atom_type;

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
	} elsif ($atom_type eq "CL") {
		$atom_mass = 35.453;
	} elsif ($atom_type eq "BR") {
		$atom_mass = 79.904;
	} elsif ($atom_type eq "I") {
		$atom_mass = 126.9045;
	} elsif ($atom_type eq "SI") {
		$atom_mass = 28.08;
	} else {
		print "Unknown atom found.  Please check the output!\n";
		print "Check the last column of your .pdb file and make sure it contains a valid atomic element symbol.\n";
		$atom_mass = 0.000;
	}

	# Print to output
	# dummy "X" printed to charge group column for later replacement
	open(ATOMS_OUT, ">>section_atoms");
	printf ATOMS_OUT "%6d%11s%7d%7s%7s%7s%11.3f%11.5f%7s%11s%11s\n", $atom_num, "opls_XXX", $res_num, $residue_name, $atom_name, "X", $atom_charge, $atom_mass;
	close(ATOMS_OUT);

	# Print to dummy output, to be used later in improper determination
	open(ATOMS_OUT_RAW, ">>section_atoms_raw");
	print ATOMS_OUT_RAW "@info\n";
	close(ATOMS_OUT_RAW);

}

# Get atom numbers (keys) and names (values) from hash
# These will be used for writing atom names to bonded sections
my @keys_atom_numbers = keys %atom_number_name;
my @values_atom_names = values %atom_number_name;

# Generate the [ bonds ] section

open(IN, $input) or die "Problem while reading $input for bonds: $!\n";
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
# TO DO:
# 	In @lines_edit (global), element 0 is atom number, element 1 is atom name
#	Open @bond_line_edit, which has atom numbers for bonds
#	Loop:
#	If $keys_numbers_names[$key] == $bond_line_edit[0], $bonded_atom_i = $values_numbers_names[$key]
#	If $lines_edit[0] == $bond_line_edit[$n] ($n >= 1), $bonded_atom_j = $lines_edit[1]
#		counter $n has to stop at $length => for loop
#	Put this all within the } else {
#

# Define new array to hold lines
my @bond_line_edit;

# Open filehandle for bonds output
# Write the heading information for the [ bonds ] section
open(BONDS_OUT, ">>section_bonds");
printf BONDS_OUT "[ bonds ]\n";
printf BONDS_OUT ";%4s%6s%6s\n", "ai", "aj", "funct";
close(BONDS_OUT);

foreach $_ (@conect) {
	chomp(@bond_line_edit = split(" ", $_));
	
	# Remove CONECT string
	shift(@bond_line_edit);
	
	# At this point, the @bond_line_edit array should hold one CONECT entry
	# Each atom number is an element in the array

	# Need to parse thru each line and satisfy these criteria
	# First time thru: $bond_line_edit[0] = 1, $bond_line_edit[1] = 2, etc.
	# Need to write out bonds based on first atom
	# If $bond_line_edit[0] > $bond_line_edit[$i], do nothing; check next element 

	my $string;

	my $bonded_atom_i;
	my $bonded_atom_j;

	foreach $_ (@bond_line_edit) {
		if ($bond_line_edit[0] >= $_) {
			# Do nothing! Includes self check - not pretty, but effective
		} else {
			# This part finds ai and aj in the bond
			for (my $key_count=0; $key_count<scalar(@keys_atom_numbers); $key_count++) {
				if (exists($arguments{"-renumber"})) {
					if ($keys_atom_numbers[$key_count] == $bond_line_edit[0] - $renumber_factor) {
						$bonded_atom_i = $values_atom_names[$key_count];
					} elsif ($keys_atom_numbers[$key_count] == $_ - $renumber_factor) {
						$bonded_atom_j = $values_atom_names[$key_count];
					}
				} else {
					if ($keys_atom_numbers[$key_count] == $bond_line_edit[0]) {
                                                $bonded_atom_i = $values_atom_names[$key_count];
                                        } elsif ($keys_atom_numbers[$key_count] == $_) {
                                                $bonded_atom_j = $values_atom_names[$key_count];
                                        }
				}
			}

			if (exists($arguments{"-renumber"})) {
				my $renumb = $bond_line_edit[0] - $renumber_factor;
				my $default_renumb = $_ - $renumber_factor;
				open(BONDS_OUT, ">>section_bonds");
				printf BONDS_OUT "%5d%6d%6d\t;%6s%6s\n", $renumb, $default_renumb, $funct_bond, $bonded_atom_i, $bonded_atom_j;
				close(BONDS_OUT);
			} else {
				open(BONDS_OUT, ">>section_bonds");	
				printf BONDS_OUT "%5d%6d%6d\t;%6s%6s\n", $bond_line_edit[0], $_, $funct_bond, $bonded_atom_i, $bonded_atom_j;
				close(BONDS_OUT);
			}
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

open(IN, "<section_bonds") or die "Problem while reading \"section_bonds\" for writing angles: $!\n";
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

	# atom numbers need to be held in variables so they can be written to the final output
	# angle is defined as i-j-k, where j is the vertex of the angle
	# atom names are stored in elements 4 and 5 (aj, ai in this context)
	my $atom_i = $line[0];
	my $atom_j = $line[1];
	my $atom_k;

	my $angle_atom_i = $line[5];
	my $angle_atom_j = $line[4];
	my $angle_atom_k;
	
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
			for (my $key_count=0; $key_count<scalar(@keys_atom_numbers); $key_count++) {
				if ($keys_atom_numbers[$key_count] == $next_line[1]) {
					$angle_atom_k = $values_atom_names[$key_count];	
				}
			}

			# Write output	
			open(ANGLES_OUT, ">>section_angles");
                        printf ANGLES_OUT "%5d%6d%6d%6d\t;%6s%6s%6s\n", $atom_i, $atom_j, $atom_k, $funct_angle, $angle_atom_i, $angle_atom_j, $angle_atom_k; 
                        close(ANGLES_OUT);

		} elsif ($next_line[0] == $atom_i) {
			$atom_k = $next_line[1];
                        for (my $key_count=0; $key_count<scalar(@keys_atom_numbers); $key_count++) {
                                if ($keys_atom_numbers[$key_count] == $next_line[1]) {
                                        $angle_atom_k = $values_atom_names[$key_count]; 
                                }
                        }

			# Write output
			open(ANGLES_OUT, ">>section_angles");
			printf ANGLES_OUT "%5d%6d%6d%6d\t;%6s%6s%6s\n", $atom_j, $atom_i, $atom_k, $funct_angle, $angle_atom_i, $angle_atom_j, $angle_atom_k;
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
open(IN2, "<section_bonds") or die "Problem while reading \"section_bonds\" for dihedrals: $!\n";
my @dihedral_input_file_bonds = <IN2>;
close(IN2);

# shift off headers
shift(@dihedral_input_file_bonds);
shift(@dihedral_input_file_bonds);
chomp(@dihedral_input_file_bonds);

# Open file containing angles
open(IN3, "<section_angles") or die "Problem while reading \"section_angles\" for dihedrals: $!\n";
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

	# define dihedral atom names (ai, aj, ak from last three fields in angle entry)
	my $dihe_angle_i = $line[5];
	my $dihe_angle_j = $line[6];
	my $dihe_angle_k = $line[7];
	my $dihe_angle_l;

	# Search for existing bond starting with $atom_k
	my $length_bonds = scalar(@dihedral_input_file_bonds);
	for (my $j=0; $j<$length_bonds; $j++) {
		chomp(my @line2 = split(" ", $dihedral_input_file_bonds[$j]));
		pop(@line2);

		if ($line2[0] == $atom_k) {
			$atom_l = $line2[1];
                        for (my $key_count=0; $key_count<scalar(@keys_atom_numbers); $key_count++) {
                                if ($keys_atom_numbers[$key_count] == $line2[1]) {
                                        $dihe_angle_l = $values_atom_names[$key_count]; 
                                }
                        }
		
			# Write dihedral output	
			# Uses a global hash %atom_number_type to determine if dihedral should be printed
			# i.e., if $i = O, $j = C, $k = N, $l = H, no dihedral (amide)
			# also no dihedral for C-C-N-H, but only in the case of amide
			# decision guards against alcohol, where you would want the dihedral:
			#	-C-(COH)-NH2 will have dihedral for C-C-N-H

			if ($atom_number_type{$atom_i} eq "O" && $atom_number_type{$atom_j} eq "C" && $atom_number_type{$atom_k} eq "N" && $atom_number_type{$atom_l} eq "H") {
				# do nothing, there should be no dihedral for an amide!
			} elsif ($atom_number_type{$atom_i} eq "C" && $atom_number_type{$atom_j} eq "C" && ($atom_number_type{$atom_j+1} eq "O" && $atom_number_type{$atom_j+2} ne "H") && $atom_number_type{$atom_k} eq "N" && $atom_number_type{$atom_l} eq "H") {
				# do nothing, there should be no dihedral for an amide!
			} else {
				open(DIHEDRALS_OUT, ">>section_dihedrals");
                        	printf DIHEDRALS_OUT "%5d%6d%6d%6d%6d\t;%6s%6s%6s%6s\n", $atom_i, $atom_j, $atom_k, $atom_l, $funct_dihedral, $dihe_angle_i, $dihe_angle_j, $dihe_angle_k, $dihe_angle_l; 
                        	close(DIHEDRALS_OUT);
			}		

			# Write pairs output
			open(PAIRS_OUT, ">>section_pairs");
			printf PAIRS_OUT "%5d%6d%6d\t;%6s%6s\n", $atom_i, $atom_l, $funct_pairs, $dihe_angle_i, $dihe_angle_l;
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
#		* "improper_A_B_C_D" written to each line

# Open file containing bonds
open(IN4, "<section_bonds") or die "Problem while reading \"section_bonds\" for impropers: $!\n";
my @improper_input_file_bonds = <IN4>;
close(IN4);

# shift off headers
shift(@improper_input_file_bonds);
shift(@improper_input_file_bonds);
chomp(@improper_input_file_bonds);

# Open file containing angles
open(IN5, "<section_angles") or die "Problem while reading \"section_angles\" for impropers: $!\n";
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

	# define variables for atom names
	my $improper_i = $line[5];
	my $improper_j = $line[6];
	my $improper_k = $line[7];
	my $improper_l;

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
	                        for (my $key_count=0; $key_count<scalar(@keys_atom_numbers); $key_count++) {
        	                        if ($keys_atom_numbers[$key_count] == $line2[1]) {
                	                        $improper_l = $values_atom_names[$key_count];
                        	        }
                        	}

			} elsif ($line2[1] == $atom_i && ($line2[0] != $atom_j && $line2[0] != $atom_k)) {
				$atom_l = $line2[0];
				for (my $key_count=0; $key_count<scalar(@keys_atom_numbers); $key_count++) {
                                	if ($keys_atom_numbers[$key_count] == $line2[0]) {
                                        	$improper_l = $values_atom_names[$key_count];
                                	}
                        	}

			}
		}

	}

	if ($counter == 3 && $atom_l != 0) {
	
		# Write improper output	
		open(IMPROPERS_OUT_RAW, ">>section_impropers_raw");
		printf IMPROPERS_OUT_RAW "%5d%6d%6d%6d%6d\t;%6s%6s%6s%6s\n", $atom_i, $atom_j, $atom_k, $atom_l, $funct_improper, $improper_i, $improper_j, $improper_k, $improper_l; 
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

	my $atom_i = $line3[0];
	my $atom_j = $line3[1];
	my $atom_k = $line3[2];
	my $atom_l = $line3[3];

	my $improper_i = $line3[6];
	my $improper_j = $line3[7];
	my $improper_k = $line3[8];
	my $improper_l = $line3[9];

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

	# NEED TO FIX THIS IN THE CASE OF CHAIN IDENTIFIERS?
	# SET EACH OF THESE TO SCALAR(ARRAY)-1, NOT ABSOLUTE POSITION	
	foreach $_ (@atom_search) {
		chomp(my @atom_search_parse = split(" ", $_));

		my $length_atom_search_parse = scalar(@atom_search_parse);

		if ($atom_i == $atom_search_parse[0]) {
			$atom_type_i = $atom_search_parse[$length_atom_search_parse-1];
		}

		if ($atom_l == $atom_search_parse[0]) {
			$atom_type_l = $atom_search_parse[$length_atom_search_parse-1];
		}

		if ($atom_k == $atom_search_parse[0]) {
			$atom_type_k = $atom_search_parse[$length_atom_search_parse-1];
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
		# apparently, the OPLS improper format requires something like: j k i l
		# The atom name output is a total hack job, just to have it make sense
		printf IMPROPERS_OUT "%5d%6d%6d%6d%6d%4s%-16s\t;%6s%6s%6s%6s\n", $atom_j, $atom_k, $atom_i, $atom_l, $funct_improper, "    ", $improper_string, $improper_l, $improper_j, $improper_k, $improper_i;    
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

#
# Primitive atom type and charge assignment:
#	1. Loop thru "section_bonds" (opened as array) - for each atom number,
#	   count instances of that number; 1 = H/carbonyl O, 2 = sp, 3 = sp2, 4 = sp3
#	2. Based on number of bonded partners and $atom_type field, assign opls_XXX type
#	3. Use splice to replace "opls_XXX" with actual atom type, 0.000 with actual charge
#

open(BONDS_REOPEN, "<section_bonds") or die "Problem while reading \"section_bonds\" for atom typing: $!\n";
my @bonds_count_array = <BONDS_REOPEN>;
close(BONDS_REOPEN);

# shift off headers
shift(@bonds_count_array);	# [ bonds ]
shift(@bonds_count_array);	# ; ai aj ...

# get max. atom number from atoms section
open(ATOMS_REOPEN, "<section_atoms") or die "Problem while reading \"section_atoms\" for atom typing: $!\n";
my @atoms_reopen = <ATOMS_REOPEN>;
close(ATOMS_REOPEN);

# get rid of extra blank line (made it pretty earlier)
pop(@atoms_reopen);

# assign the last line of the atoms section to a string
my $atom_last_line = @atoms_reopen[scalar(@atoms_reopen - 1)];

# split the string into a new array; first element will hold the atom number
my @atom_last_line_array = split(" ", $atom_last_line);

# return the maximum atom number
my $max_atom_num = $atom_last_line_array[0];

#
# Now for some nested loops
#	For each line in the bonds array, there are two elements - atom numbers for ai and aj
#	For each element, if element matches a counter (from 1 -> $max_atom_num), add to counter
#

# initialize variable for net charge
my $net_charge = 0;

# make another atom type hash with atom numbers and $atom_type_to_write (within for-loop)
# so that better atom types can be assigned, for e.g. H atoms
my %atom_numbers_final_types;

# EXPERIMENTAL - charge group assignment
# initialize charge group, increment whenever new C or N atom type identified
# Potential problem:  if first atom is not C or N, will keep placeholder "X" as charge group
my $charge_group = 0;

# start a new loop that runs through each of the atom numbers
for (my $atom = 1; $atom <= $max_atom_num; $atom++) {
	my $count = 0;

	my $atom_type_to_write;
	my $atom_charge_to_write;
	my $atom_type_final;

	foreach $_ (@bonds_count_array) {
	        my @line = split(" ", $_);
		# remove some info from the end of the line - makes searching easier later
	        pop(@line);     # atom j (name)
	        pop(@line);     # atom i (name)
	        pop(@line);     # semicolon
	        pop(@line);     # function type

		# @line now contains two elements - atom numbers for ai and aj
		# determine number of bonds to each atom, if atom appears as atom i or j
		if ($line[0] == $atom || $line[1] == $atom) {
			$count++;
		}

		# now open global hash that contains atom numbers as keys, types as values
		my $k;		# variable for searching through keys

		# Search through hash
		# if key == atom number, then atom type = value
		foreach $k (keys (%atom_number_type)) {
			if ($k == $atom) {
				$atom_type_final = $atom_number_type{$atom};
			}
		}

	}

	# Use $atom+1/+2 to search for subsequent atom types (like N)?
	#	Example: alcohol
	#		CH2OH: $atom_number_type{$atom} = C, {$atom+1/+2} = H, {$atom+3} = O, etc.

	### SP3 Carbons ###
	if ($count == 4 && $atom_type_final eq "C") {
		
		### ALCOHOLS ###
		# covers R-COH, CH3OH, R-CH2OH
		if (($atom_number_type{$atom+1} eq "O" && $atom_number_type{$atom+2} eq "H") || ($atom_number_type{$atom+1} eq "H" && $atom_number_type{$atom+2} eq "H" && $atom_number_type{$atom+3} eq "H" && $atom_number_type{$atom+4} eq "O" && $atom_number_type{$atom+5} eq "H") || ($atom_number_type{$atom+1} eq "H" && $atom_number_type{$atom+2} eq "H" && $atom_number_type{$atom+3} eq "O" && $atom_number_type{$atom+4} eq "H")) {
			$atom_type_to_write = "opls_157";
			$atom_charge_to_write = 0.145;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;
		
		# CHOH
		} elsif ($atom_number_type{$atom+1} eq "H" && $atom_number_type{$atom+2} eq "O" && $atom_number_type{$atom+3} eq "H") {
			$atom_type_to_write = "opls_158";
			$atom_charge_to_write = 0.205;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		### AMINES ###
		# covers R-CNH, CH3N-R, R-CH2N-R
		} elsif (($atom_number_type{$atom+1} eq "N" && $atom_number_type{$atom+2} eq "H") || ($atom_number_type{$atom+1} eq "H" && $atom_number_type{$atom+2} eq "H" && $atom_number_type{$atom+3} eq "H" && $atom_number_type{$atom+4} eq "N") || ($atom_number_type{$atom+1} eq "H" && $atom_number_type{$atom+2} eq "H" && $atom_number_type{$atom+3} eq "N")) {
			$atom_type_to_write = "opls_906";
			$atom_charge_to_write = 0.060;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		### THIOLS ###
		# covers R-CSH, CH3SH, R-CH2SH
		} elsif (($atom_number_type{$atom+1} eq "S" && $atom_number_type{$atom+2} eq "H") || ($atom_number_type{$atom+1} eq "H" && $atom_number_type{$atom+2} eq "H" && $atom_number_type{$atom+3} eq "H" && $atom_number_type{$atom+4} eq "S" && $atom_number_type{$atom+5} eq "H") || ($atom_number_type{$atom+1} eq "H" && $atom_number_type{$atom+2} eq "H" && $atom_number_type{$atom+3} eq "S" && $atom_number_type{$atom+4} eq "H")) {
			$atom_type_to_write = "opls_206";
			$atom_charge_to_write = 0.060;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		### FLUORINE ###
		# covers any neighboring fluorine, could be dangerous
		} elsif (($atom_number_type{$atom-1} eq "F") || ($atom_number_type{$atom+1} eq "F")) {
			$atom_type_to_write = "opls_161";
			$atom_charge_to_write = 0.532;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		### alpha-C IN R-C-COO ###
		# additional -0.100 charge relative to alkane (i.e., -0.220 for CH2, etc.)
		# R-CH2-COO
		# R-CH2-COOH defaults to alkane CH2 (next)
		} elsif ($atom_number_type{$atom+1} eq "H" && $atom_number_type{$atom+2} eq "H" && $atom_number_type{$atom+3} eq "C" && $atom_number_type{$atom+4} eq "O" && $atom_number_type{$atom+5} eq "O" && $atom_number_type{$atom+6} ne "H") {
			$atom_type_to_write = "opls_274";
			$atom_charge_to_write = -0.220;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		### ALKANE CH2 ###
		# R-CH2-R
		} elsif ($atom_number_type{$atom+1} eq "H" && $atom_number_type{$atom+2} eq "H" && $atom_number_type{$atom+3} ne "H") {
			$atom_type_to_write = "opls_136";	# alkane CH2
			$atom_charge_to_write = -0.120;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		### ALKANE CH ###
		# R-CH-R
		} elsif ($atom_number_type{$atom+1} eq "H" && $atom_number_type{$atom+2} ne "H") {
			$atom_type_to_write = "opls_137";	# alkane CH
			$atom_charge_to_write = -0.060;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		### TERT-BUTYL ###
		} elsif ($atom_number_type{$atom+1} eq "C" && $atom_number_type{$atom+2} eq "C" && $atom_number_type{$atom+3} eq "C") {
			$atom_type_to_write = "opls_516";
			$atom_charge_to_write = 0.000;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		### DEFAULT - ALKANE ###
		} else {
			$atom_type_to_write = "opls_135";	# alkane CH3
			$atom_charge_to_write = -0.180;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;
		}

	### SP2 Carbons ###
	} elsif ($count == 3 && $atom_type_final eq "C") {

		### AROMATIC ALCOHOL (e.g. PHENOL) ###
		# covers R-COH
		if ($atom_number_type{$atom+1} eq "O" && $atom_number_type{$atom+2} eq "H") {
			$atom_type_to_write = "opls_166";
			$atom_charge_to_write = 0.150;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		### AMIDE or KETONE ###
		# covers R-CONH, R-CO-R
		} elsif (($atom_number_type{$atom+1} eq "O" && $atom_number_type{$atom+2} eq "N") || ($atom_number_type{$atom+1} eq "O" && $atom_number_type{$atom+2} eq "C" && $atom_number_type{$atom-1} eq "C")) {
			$atom_type_to_write = "opls_235";
			$atom_charge_to_write = 0.500;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		### ESTER/CARBOXYLIC ACID ###
		# covers R-COO, not R-COOH
		} elsif (($atom_number_type{$atom+1} eq "O" && $atom_number_type{$atom+2} eq "O" && $atom_number_type{$atom+3} ne "H") || ($atom_number_type{$atom-1} eq "O" && $atom_number_type{$atom+1} eq "O" && $atom_number_type{$atom+2} ne "H")) {
			$atom_type_to_write = "opls_271";
			$atom_charge_to_write = 0.700;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		# covers R-COOH (assumes H is part of carboxylate, not another moiety
		} elsif ($atom_number_type{$atom+1} eq "O" && $atom_number_type{$atom+2} eq "O" && $atom_number_type{$atom+3} eq "H") {
			$atom_type_to_write = "opls_267";
			$atom_charge_to_write = 0.520;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		### alpha-C in C(ring)-COO
		} elsif ($atom_number_type{$atom+1} eq "C" && $atom_number_type{$atom+2} eq "O" && $atom_number_type{$atom+3} eq "O" && $atom_number_type{$atom+4} ne "H") {
			$atom_type_to_write = "opls_145";
			$atom_charge_to_write = -0.100;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		### AROMATIC C WITH NO H ###
		# i.e., ring C connected to amide, amine, alkyl chain
		# extra check of atom+1 = C/N necessary so as not to over-write, e.g. phenol
		} elsif ($atom_number_type{$atom+1} ne "H" && ($atom_number_type{$atom+1} eq "C" || $atom_number_type{$atom+1} eq "N")) {
			$atom_type_to_write = "opls_145";
			$atom_charge_to_write = 0.000;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		### DEFAULT - RING CH ###
		} else {
			$atom_type_to_write = "opls_145";	# based on aromatic C
			$atom_charge_to_write = -0.115;
			$net_charge += $atom_charge_to_write;
      		        $atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;
		}

	### SP3 Nitrogen ###
	# covers only quaternary amine (NR4+)
	} elsif ($count == 4 && $atom_type_final eq "N") {
		$atom_type_to_write = "opls_287";	# based on quaternary amine
		$atom_charge_to_write = -0.300;
		$net_charge += $atom_charge_to_write;
                $atom_numbers_final_types{$atom} = $atom_type_to_write;
		$charge_group++;	

	### SP2 Nitrogen ###
	} elsif ($count == 3 && $atom_type_final eq "N") {
	
		### AMIDE ###
		# covers amide R-CONH such that carbonyl C is two atoms back, O is one atom back
		if ($atom_number_type{$atom-2} eq "C" && $atom_number_type{$atom-1} eq "O") {
			$atom_type_to_write = "opls_237";	# based on amide
			$atom_charge_to_write = -0.760;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;	

		### SECONDARY AMINE ###
		# covers R-NH-R (flanked by C)
		} elsif ($atom_number_type{$atom-1} eq "C" && $atom_number_type{$atom+1} eq "H" && $atom_number_type{$atom+2} eq "C") {
			$atom_type_to_write = "opls_901";
			$atom_charge_to_write = -0.310;		# guess
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;

		### PRIMARY AMINE ###
		} else {
			$atom_type_to_write = "opls_900";
			$atom_charge_to_write = -0.900;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
			$charge_group++;
		}

	### SP2 Oxygen w/2 bonded neighbors ###
	# covers phenol-OH and alkyl-OH oxygens
	} elsif ($count == 2 && $atom_type_final eq "O") {

		### AROMATIC GROUP (e.g., PHENOL) ###
		if ($atom_numbers_final_types{$atom-1} eq "opls_166") {
			$atom_type_to_write = "opls_167";	# based on TYR -OH
			$atom_charge_to_write = -0.585;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		# R-COOH: -OH oxygen
		} elsif ($atom_number_type{$atom-2} eq "C" && $atom_number_type{$atom-1} eq "O" && $atom_number_type{$atom+1} eq "H") {
			$atom_type_to_write = "opls_268";
			$atom_charge_to_write = -0.440;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write; 
	
		# R-COO-R (ester O)
		} elsif (($atom_number_type{$atom+1} eq "C" && $atom_number_type{$atom+2} eq "O" && $atom_number_type{$atom+3} eq "C") || ($atom_number_type{$atom-2} eq "C" && $atom_number_type{$atom-1} eq "O" && $atom_number_type{$atom+1} eq "C")) {
			$atom_type_to_write = "opls_467";
			$atom_charge_to_write = -0.350;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		### STANDARD ALCOHOL ###
		} else {
			$atom_type_to_write = "opls_154";	# based on alcohol
			$atom_charge_to_write = -0.683;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
		}

	### SP2 Oxygen w/1 bonded neighbor ###
	# covers carbonyl, COO-, guesses at phosphate
	} elsif ($count == 1 && $atom_type_final eq "O") {

		### CARBOXYLATE ###
		# R-COO- only
		if (($atom_number_type{$atom-2} eq "C" && $atom_number_type{$atom-1} eq "O" && $atom_number_type{$atom+1} ne "H") || ($atom_number_type{$atom-1} eq "C" && $atom_number_type{$atom+1} eq "O" && $atom_number_type{$atom+2} ne "H")) {
			$atom_type_to_write = "opls_272";
			$atom_charge_to_write = -0.800;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		# R-COOH: carbonyl oxygen
		} elsif ($atom_number_type{$atom-1} eq "C" && $atom_number_type{$atom+1} eq "O" && $atom_number_type{$atom+2} eq "H") {
			$atom_type_to_write = "opls_269";
			$atom_charge_to_write = -0.530;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		# R-COO-R: ester oxygen (carbonyl)
		} elsif (($atom_number_type{$atom-1} eq "C" && $atom_number_type{$atom+1} eq "O" && $atom_number_type{$atom+2} eq "C") || ($atom_number_type{$atom-1} eq "C" && $atom_number_type{$atom-2} eq "O" && $atom_number_type{$atom+1} eq "C")) {
			$atom_type_to_write = "opls_466";
			$atom_charge_to_write = -0.350;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		### PHOSPHATE ###
		} elsif (($atom_number_type{$atom+1} eq "P") || ($atom_number_type{$atom-1} eq "P") || ($atom_number_type{$atom-2} eq "P") || ($atom_number_type{$atom-3} eq "P")) {
			$atom_type_to_write = "opls_441";
			$atom_charge_to_write = -0.920;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		### DEFAULT - GENERIC CARBONYL ###
		} else {
			$atom_type_to_write = "opls_236";	# carbonyl
			$atom_charge_to_write = -0.500;
			$net_charge += $atom_charge_to_write;
               		$atom_numbers_final_types{$atom} = $atom_type_to_write;
		}

	### SULFUR ###
	# only covers R-SH
	} elsif ($count == 2 && $atom_type_final eq "S") {
		
		### THIOL ###
		if ($atom_number_type{$atom-1} eq "C" && $atom_number_type{$atom+1} eq "H") {
			$atom_type_to_write = "opls_200";	# thiol
			$atom_charge_to_write = -0.335;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
		### GENERIC S IN ALKYL GROUP, i.e. MET ###
		} else {
			$atom_type_to_write = "opls_202";	# methionine sidechain
			$atom_charge_to_write = -0.335;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
		}


	### HYDROGEN ATOMS ###
	} elsif ($count == 1 && $atom_type_final eq "H") {

		# all decisions based on preceding atom only
		# if alkane, then sequence is usually C, H, H, H - first H assigned based on
		# the C; subsequent H defined by the preceding H, etc.

		### AROMATIC H ###
		if ($atom_numbers_final_types{$atom-1} eq "opls_145") {
			$atom_type_to_write = "opls_146";	# aromatic H
			$atom_charge_to_write = 0.115;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		### ALKANE H ###
		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_135" || $atom_numbers_final_types{$atom-1} eq "opls_136" || $atom_numbers_final_types{$atom-1} eq "opls_137" || $atom_numbers_final_types{$atom-1} eq "opls_274") {
			$atom_type_to_write = "opls_140";	# alkane H
			$atom_charge_to_write = 0.060;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_140") {
			$atom_type_to_write = "opls_140";
			$atom_charge_to_write = 0.060;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		### QUATERNARY AMINE ###
		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_287") {
			$atom_type_to_write = "opls_290";	# quaternary amine H (LYSH NZ)
			$atom_charge_to_write = 0.330;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_290") {
			$atom_type_to_write = "opls_290";
			$atom_charge_to_write = 0.330;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		### TERTIARY AMINE ###
		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_900") {
			$atom_type_to_write = "opls_909";	# tertiary amine H (LYS NZ)
			$atom_charge_to_write = 0.360;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_909") {
			$atom_type_to_write = "opls_909";
			$atom_charge_to_write = 0.360;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		### SECONDARY AMINE ###
		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_901") {
			$atom_type_to_write = "opls_910";
			$atom_charge_to_write = 0.310;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		### AMIDE ###
		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_237") {
			$atom_type_to_write = "opls_240";	# amide H
			$atom_charge_to_write = 0.380;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_240") {
			$atom_type_to_write = "opls_240";
			$atom_charge_to_write = 0.380;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		### ALCOHOLS ###
		# aromatic
		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_167") {
			$atom_type_to_write = "opls_168";
			$atom_charge_to_write = 0.435;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
		
		# alkyl
		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_154") {
			$atom_type_to_write = "opls_155";	# alcohol H
			$atom_charge_to_write = 0.418;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_158") {
			$atom_type_to_write = "opls_140";
			$atom_charge_to_write = 0.060;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		### CARBOXYLIC ACID ###
		# R-COOH
		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_268") {
			$atom_type_to_write = "opls_270";
			$atom_charge_to_write = 0.450;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
	
		### PHOSPHATE ###
		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_441") {
			$atom_type_to_write = "opls_435";	# H in OH-
			$atom_charge_to_write = 0.300;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;
	
		### THIOL ###
		} elsif ($atom_numbers_final_types{$atom-1} eq "opls_200") {
			$atom_type_to_write = "opls_204";	# thiol H
			$atom_charge_to_write = 0.155;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{$atom} = $atom_type_to_write;

		### DEFAULT OUTPUT ###
		} else {
			print "Using default H type for atom $atom.\n";
			$atom_type_to_write = "opls_140";	# alkane H (generic output)
			$atom_charge_to_write = 0.060;
			$net_charge += $atom_charge_to_write;
			$atom_numbers_final_types{atom} = $atom_type_to_write;
		}

	} elsif ($count == 1 && $atom_type_final eq "F") {
		$atom_type_to_write = "opls_164";		# generic F from ffoplsaa.atp
		$atom_charge_to_write = -0.206;
		$net_charge += $atom_charge_to_write;
                $atom_numbers_final_types{$atom} = $atom_type_to_write;

	} elsif ($count == 1 && $atom_type_final eq "CL") {
		$atom_type_to_write = "opls_151";		# alkyl chloride
		$atom_charge_to_write = -0.200;
		$net_charge += $atom_charge_to_write;
		$atom_numbers_final_types{$atom} = $atom_type_to_write;

	} elsif ($count == 1 && $atom_type_final eq "BR") {
		$atom_type_to_write = "opls_722";		# bromobenzene
		$atom_charge_to_write = -0.220;
		$net_charge += $atom_charge_to_write;
		$atom_numbers_final_types{$atom} = $atom_type_to_write;

	} elsif ($count == 1 && $atom_type_final eq "I") {
		$atom_type_to_write = "opls_732";		# iodobenzene
		$atom_charge_to_write = -0.100;
		$net_charge += $atom_charge_to_write;
		$atom_numbers_final_types{$atom} = $atom_type_to_write;
	
	} elsif ($atom_type_final eq "P") {
		$atom_type_to_write = "opls_440";		# generic P from ffoplsaa.atp
		$atom_charge_to_write = 1.620;
		$net_charge += $atom_charge_to_write;
                $atom_numbers_final_types{$atom} = $atom_type_to_write;

	} elsif ($atom_type_final eq "SI") {
		$atom_type_to_write = "SI";
		$atom_charge_to_write = 0.000;
                $atom_numbers_final_types{$atom} = $atom_type_to_write;
		print "Found SI atom, cannot assign charge. Check the output and make an appropriate correction.\n";

	} else {
		print "Unknown atom type: $atom_type_final $atom with $count bonds - cannot assign atom type.\n";
	}

	# Now, open up the atoms section and make the replacements
	open(ATOMS_FINAL_MOD, "<section_atoms") or die "Problem while reading \"section_atoms\" while making atom type replacements: $!\n";
	my @atoms_final_mod = <ATOMS_FINAL_MOD>;
	close(ATOMS_FINAL_MOD);

	# Try making the substitution using regular expression pattern matching?
	foreach $_ (@atoms_final_mod) {
		if ($_ =~ /^;/ || $_ =~ /\s*\[/ || $_ =~ /\s*#/) {
			# do nothing, this is a comment line
			# probably unnecessary, now that header is written separately
		} elsif ($_ =~ /\s+$atom\s+opls_XXX/) {
			$_ =~ s/opls_XXX/$atom_type_to_write/;
			$_ =~ s/0\.000/$atom_charge_to_write/;
			$_ =~ s/X/$charge_group/;

			# split line and write to output
			my @final_line = split(" ", $_);

			open(ATOMS_FINAL_OUTPUT, ">>section_atoms_final") or die "Problem opening \"section_atoms_final\": $!\n";
			printf ATOMS_FINAL_OUTPUT "%6d%11s%7d%7s%7s%7d%11.3f%11.5f%7s%11s%11s\n", $final_line[0], $final_line[1], $final_line[2], $final_line[3], $final_line[4], $final_line[5], $final_line[6], $final_line[7];
			close(ATOMS_FINAL_OUTPUT);

		}

	}

}

# Print warning for the user if we have a net charge
if ($net_charge != 0) {
	print "\n********************************************************\n";
	print "WARNING: Net charge on the molecule is $net_charge\n";
	print "You will need to make corrections to the output topology\n";
	print "********************************************************\n";
} else {
	# Things still may not be OK...
	print "\n********************************************************\n";
	print "The net charge on your molecule is 0. HOWEVER! This does\n";
	print "not necessarily indicate that the topology is correct\n";
	print "********************************************************\n";
}

# Concatenate the final output into one single file

# Open output filehandle
open (TOPOLOGY_OUT, ">>$output") or die "Cannot open $output: $!\n";

# Open input files to be written to output

# header
open(HEAD_IN, "<section_header") or die "Cannot open \"section_header\" while writing output: $!\n";

while (my $line_out = <HEAD_IN>) {
	print TOPOLOGY_OUT $line_out;
}

# [ atoms ] section
open(FIRST_IN, "<section_atoms_final") or die "Cannot open \"section_atoms_final\" while writing output: $!\n";

while (my $line_out = <FIRST_IN>) {
	print TOPOLOGY_OUT $line_out;
}

print TOPOLOGY_OUT "\n";	# add spacing

close(FIRST_IN);

# [ bonds ] section
open(SECOND_IN, "<section_bonds") or die "Cannot open \"section_bonds\" while writing output: $!\n";

while (my $line_out = <SECOND_IN>) {
	print TOPOLOGY_OUT $line_out;
}

close(SECOND_IN);

# [ pairs ] section
open(THIRD_IN, "<section_pairs") or die "Cannot open \"section_pairs\" while writing output: $!\n";

while (my $line_out = <THIRD_IN>) {
	print TOPOLOGY_OUT $line_out;
}

close(THIRD_IN);

# [ angles ] section
open(FOURTH_IN, "<section_angles") or die "Cannot open \"section_angles\" while writing output: $!\n";

while (my $line_out = <FOURTH_IN>) {
	print TOPOLOGY_OUT $line_out;
}

close(FOURTH_IN);

# proper [ dihedrals ] section
open(FIFTH_IN, "<section_dihedrals") or die "Cannot open \"section_dihedrals\" while writing output: $!\n";

while (my $line_out = <FIFTH_IN>) {
	print TOPOLOGY_OUT $line_out;
}

close(FIFTH_IN);

# improper [ dihedrals ] section
# only write it out if it contains relevant lines
open (TEST, "<section_impropers") or die "Cannot open \"section_impropers\" while writing output: $!\n";

chomp(my @test_array = <TEST>);

my $test_length = scalar(@test_array);

if ($test_length > 3) {
	open(SIXTH_IN, "<section_impropers") or die "Cannot open \"section_impropers\" while testing for output: $!\n";

	while (my $line_out = <SIXTH_IN>) {
		print TOPOLOGY_OUT $line_out;
	}

	close(SIXTH_IN);
}

if ($out_type eq "top") {
	print TOPOLOGY_OUT "[ system ]\n";
	print TOPOLOGY_OUT "; Name\n";
	print TOPOLOGY_OUT "$res_name topology, generated by TopolGen\n";
	print TOPOLOGY_OUT "\n";
	print TOPOLOGY_OUT "[ molecules ]\n";
	printf TOPOLOGY_OUT ";%9s%13s\n", "Compound", "\#mols";
	printf TOPOLOGY_OUT "%-6s%15d\n", $res_name, 1;
}

close(TOPOLOGY_OUT);

# Clean up intermediate files
unlink "section_header";
unlink "section_atoms";
# unlink "section_atoms_raw";
unlink "section_atoms_final";
unlink "section_bonds";
unlink "section_angles";
unlink "section_pairs";
unlink "section_dihedrals";
unlink "section_impropers";

print "
TopolGen, version $v complete.

Output topology has been written. An attempt has been made to assign charges and atom types based
on existing functional groups, but they may not be correct.  No charge calculations or other 
parameterization calculations have been done.  Guesses have been made for charge groups.  Please
inspect and correct the topology before using it in any simulations.  The author of the script does 
NOT guarantee accuracy or usability of any of the content; TopolGen was written as a convenience 
for outputting a skeleton topology, and nothing more.\n\n";

exit;
