#!/usr/bin/env perl

#######################################################
## Modified version of the script from proteinortho. ##
#######################################################

#pk

##########################################################################################
#	  This file is part of proteinortho.
#	  (C) 2009 Marcus Lechner
# 
#	  proteinortho is free software; you can redistribute it and/or modify
#	  it under the terms of the GNU General Public License as published
#	  by the Free Software Foundation; either version 2, or (at your
#	  option) any later version.
#
#	  proteinortho is distributed in the hope that it will be useful, but
#	  WITHOUT ANY WARRANTY; without even the implied warranty of
#	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#	  General Public License for more details.
#
#	  You should have received a copy of the GNU General Public License
#	  along with proteinortho; see the file COPYING.  If not, write to the
#	  Free Software Foundation, Inc., 59 Temple Place - Suite 330,
#	  Boston, MA 02111-1307, USA.	
##########################################################################################

##########################################################################################
# About
##########################################################################################
# 
# @author Paul Klemm
# @email klemmp at staff.uni-marburg.de
# @company Bioinformatics, University of Leipzig
# @version 4
# @date 3-12-2020
#
##########################################################################################
#
# @contributor Andrey Rozenberg
# @email alephreish at gmail.com
# @company Technion - Israel Institute of Technology
# @version 4AR
# @date 22-11-2021
#
##########################################################################################

use POSIX;

my $usage = <<'ENDUSAGE';
proteinortho_grab_proteins.pl        greps all genes/proteins of a given fasta file
 
SYNOPSIS
 
proteinortho_grab_proteins.pl (options) QUERY FASTA1 (FASTA2 ...)

	QUERY	proteinortho.tsv FILE or search STRING or '-' for STDIN:
		a)	proteinortho output file (.tsv). This uses by default the -exact option.
		b)	string of one identifier e.g. 'tr|asd3|asd' OR multiple identifier separated by ',' (-F=)
	FASTA*	fasta file(s) (database)

	(options):
		-tofiles, -t  print everything to files instead of stdout files are called OrthoGroup**.fasta for a proteinortho.tsv file
		-E            enables regex matching otherwise the string is escaped (e.g. | -> \|)
		-exact        search patters are extended with a \b, that indicates end of word.
		-source, -s   adds the filename (FASTA1,...) to the found gene-name
                -source_only  adds the filename (FASTA1,...) without extension instead of the found gene-name
                -output_dir=s output directory
                -output_file=s save each output in a subfolder as file with this fixed name
		-F=s          char delimiter for multiple identifier if QUERY is a string input (default: ',')

DESCRIPTION
 
	This script finds and extract all given identifier of a list of fasta files. 
	The identifier can be provided as a simple string 'BDNF1', regex string 'BDNF*' 
	or in form of a proteinortho output file (myproject.proteinortho.tsv).
       
	Example:
 
 	# 1. most simple call:

	perl proteinortho_grab_proteins.pl 'BDNF1' *.faa

		STDOUT:
			>BDNF1 Brain derived neurotrophic factor OS=human(...)
			MNNGGPTEMYYQQHMQSAGQPQQPQTVTSGPMSHYPPAQPPLLQPGQPYSHGAPSPYQYG
			>BDNF15 Brain derived neurotrophic factor OS=human(...)
			MAFPLHFSREPAHAIPSMKAPFSRHEVPFGRSPSMAIPNSETHDDVPPPLPPPRHPPCTN

	    The second hit BDNF15 is reported since it also contains 'BDNF1' as a substring. 
	    To prevent such a behaviour use proteinortho_grab_proteins.pl -E 'BDNF1\b'. 
	    The \b marks the end of a word and -E enables regex expressions.

	    Or simply add -exact: perl proteinortho_grab_proteins.pl -exact 'BDNF1' *.faa

 	# 2. multiple ids:

	perl proteinortho_grab_proteins.pl 'BDNF1,BDNF2,BDNF3' *.faa

 	# 3. more complex regex search:

	perl proteinortho_grab_proteins.pl -E 'B?DNF[0-3]3+' *.faa

		This finds: BDNF13, BDNF23, DNF13, DNF033, ... 

 	# 4. proteinortho tsv file and write output to files:

	proteinortho_grab_proteins.pl -tofiles myproject.proteinortho.tsv test/*.faa

		This will produce the files: OrthoGroup0.fasta, OrthoGroup1.fasta, OrthoGroup2.fasta, ...
		Each fasta file contains all genes of one orthology group (one line in myproject.proteinortho.tsv)
 
ENDUSAGE

my $query;
my $help;
my $tofiles=0;
my $justid;
my $prefix=">";
my $doregex=0;
my $source=0;
my $source_only=1;
my $exact=0;
my $del=',';
my @ARGV_copy=@ARGV;
my $output_dir = ".";
my $output_file = "";

my @ARGV_copyiddone=(1) x (scalar(@ARGV_copy));
my $ARGV_copyiddone_counter=scalar(@ARGV_copy);
for(my $v = 0 ; $v < scalar @ARGV_copy ; $v++){
	$ARGV_copyiddone[$v]=1;
	if($ARGV_copy[$v] =~ m/--?(help|h)$/){$help=1;}
	elsif($ARGV_copy[$v] =~ m/^--?(tofiles|t)$/){$tofiles=1;}
	elsif($ARGV_copy[$v] =~ m/^--?(source|s)$/){$source=1;}
	elsif($ARGV_copy[$v] =~ m/^--?source_only$/){$source_only=1;}
	elsif($ARGV_copy[$v] =~ m/^--?F=(.*)$/){$del=$1;}
	elsif($ARGV_copy[$v] =~ m/^--?E$/){$doregex=1;}
	elsif($ARGV_copy[$v] =~ m/^--?exact$/){$exact=1;}
        elsif($ARGV_copy[$v] =~ m/^--?output_dir=(.*)$/){$output_dir=$1;}
	elsif($ARGV_copy[$v] =~ m/^--?output_file=(.*)$/){$output_file=$1;}
	elsif($ARGV_copy[$v] =~ m/^-.+/){print $usage; print STDERR "ERROR: invalid option ".$ARGV_copy[$v]."!\n\n";exit(1);}
	elsif(!defined($query)){$query = $ARGV_copy[$v];}
	else{$ARGV_copyiddone[$v]=0;$ARGV_copyiddone_counter--;}
}
if ($help){
    print $usage;
    exit(0);
}
my $fail="";
if (scalar(@ARGV_copy) == 0){
    $fail.="ERROR: no arguments found!\n";
}
if (!defined($query)){
    $fail.="ERROR: no QUERY provided!\n";
}
if ( $ARGV_copyiddone_counter==scalar(@ARGV_copy) ){
    $fail.="ERROR: no FASTA files provided!\n";
}
if($fail){
	print $usage.$fail;
	exit(1);
}
if($tofiles and $output_dir ne "." and not -d $output_dir){
	mkdir $output_dir or die "Error creating directory: $output_dir";
}

my %qdata;
# my $qdata_count = {};

our $orthogroupcounter=0;
our $genecounter=0;
my $numOfFastas=0;
my $line_i = 0;

for(my $v = 0 ; $v < scalar @ARGV_copy ; $v++){
	if($ARGV_copyiddone[$v]){next;}
	$numOfFastas++;
}
our @filenames;

my $foundHeader=0;

# print STDERR (scalar keys %qdata)." vs ".($numOfFastas);

sub processLine{
	my $line = shift;
	my $prefix = shift;

	$line=~s/[\r\n]+$//; 

	my @sp = split(/\t/,$line);
	if(substr($line,0,1) eq "#"){$foundHeader=1;@filenames=@sp; next;}
	if(scalar(@sp)>3){
		for(my $v = 3 ; $v < scalar @sp ; $v++){
			if($sp[$v] eq "*" || $sp[$v] eq ""){next;}
			my @spp = split(",",$sp[$v]);

			for(my $vv = 0 ; $vv < scalar @spp ; $vv++){

				if($spp[$vv] eq "*" || $spp[$vv] eq ""){next;}
				$spp[$vv]=~s/^\(//;$spp[$vv]=~s/\)$//;

				if(!exists $filenames[$v]){
					$filenames[$v]=""
				}

				# $qdata{$filenames[$v]}{$spp[$vv]}=$prefix.".OrthoGroup".$orthogroupcounter;
				$qdata{$filenames[$v]}{$spp[$vv]} = "OrthoGroup".$orthogroupcounter;
				$genecounter++;
			}
		}
	}
	$orthogroupcounter++;
}

unless(open(my $FH,'<',$query)) {
	if($query eq "-"){
		foreach my $line (<STDIN>) { &processLine($line, "STDIN") }
	}else{
		my @sp = split($del,$query); 
		for(my $v = 0 ; $v < scalar @sp ; $v++){
			if(length($sp[$v])==0 || $sp[$v] eq ""){next;}
			$sp[$v]=~s/^\(//;$sp[$v]=~s/\)$//;

			$qdata{"STDIN"}{$sp[$v]}=$sp[$v];
			$genecounter++;
		}
	}
}else{
	if(!$exact){print STDERR "[STDERR] WARNING The -exact option is mandatory if a proteinortho file is given. -exact is now set.\n";$exact=1;}
	if($doregex){print STDERR "[STDERR] WARNING The -E option is not allowed if a proteinortho file is given. -E is now unset.\n";$doregex=0;}

	my $query_basename=$query;
	if($query_basename =~ m/\/([^\/]+)$/){
		$query_basename=$1;
	}

	while(<$FH>){ &processLine($_, $query_basename) }
	close($FH);
	print STDERR "[STDERR] Done reading the query $query file. Now I know $orthogroupcounter groups with $genecounter genes/proteins in total.\n";
}


if( $foundHeader==0 && $numOfFastas > 3 && $genecounter > 20){
	print STDERR "\nWARNING : The header of the proteinortho file is missing, this can increase the runtime dramatically. Please include the first line (starting with '#'), to accelerate this program.\n$NC\n";
	sleep 1;
}

#if( $tofiles==1 && ($orthogroupcounter > 100) ){
#	print STDERR "\n!!!\nWARNING : This call will produce $orthogroupcounter files (one for each orthology group) !\nIn the *.html file you can individually extract single groups by clicking on the front part of a row.\n$NC";
#	print STDERR "Press 'strg+c' to prevent me from proceeding or wait 20 seconds to continue...\n!!!\n";
# 	sleep 20;
#	print STDERR "\nWell then, proceeding...\n\n";
#	sleep 1;
#}

my $cur_gene="";
my $cur_gene_filename="";
my %cur_gene_firsttime;
my $genecounterfound=0;
my $basename = "";

my $fastai=1;

my %foundIDs;

for(my $v = 0 ; $v < scalar @ARGV_copy ; $v++){

	if($ARGV_copyiddone[$v]){next;}

	print STDERR "[STDERR] ($fastai/$numOfFastas) : ";if($basename ne ""){print STDERR "Done reading $basename. "}print STDERR "Start reading the fasta file ".($ARGV_copy[$v])."\n";
	$fastai++;

	$basename = $ARGV_copy[$v];
	if($basename =~ m/\/([^\/]+)$/){
		$basename=$1;
	}

	open(my $FH,'<',$ARGV_copy[$v]);
	my $geneprintswitch = 0;
	
	while(<$FH>){
		$_=~s/[\r\n]+$//;

		if($_ eq "" || length $_ < 2 || substr($_,0,1) eq "#"){next;}
		
		my $curLine=$_;

		if(substr($curLine,0,1) eq $prefix ){
			$geneprintswitch = 0;

			if($cur_gene ne ""){
				$cur_gene_filename=~s/\\b//g;
				$cur_gene_filename=~s/[^a-zA-Z0-9.]//g;
				if($tofiles){ # print to files
                                        if($output_file ne ""){
						my $dir_name = "$output_dir/$cur_gene_filename";
						$cur_gene_filename = "$dir_name/$output_file";
						if (not -d $dir_name) {
							mkdir $dir_name or die "Error creating directory: $dir_name";
						}
					} else {
						$cur_gene_filename = "$output_dir/$cur_gene_filename.fasta";
					}
					my $writemodus=">>";
					if(!exists $cur_gene_firsttime{$cur_gene_filename}){$writemodus=">";$cur_gene_firsttime{$cur_gene_filename}=1;}
					open($FHOUT,$writemodus,$cur_gene_filename);
					print $FHOUT $cur_gene;
					close($FHOUT);
				}else{
					print $cur_gene;
				}
				$cur_gene="";
			}

			if($exact && exists $qdata{$basename}){

				my $genename=$curLine;
				my @arr=split(" ",$genename);
				if(scalar(@arr)>0){$genename=$arr[0];}
				$genename=~s/^>//;

				if(exists $qdata{$basename}{$genename}){
					
					my $headerstr=$curLine;
					if($source){$headerstr=$headerstr." ".$basename;}
					if($source_only) {
						$headerstr = ">" . ($basename =~ s/[.]\w+$//r);
					}

					$cur_gene.=$headerstr."\n";
					$cur_gene_filename=$qdata{$basename}{$genename}; # .".fasta";
					$geneprintswitch = 1;
					$genecounterfound++;

					delete $qdata{$basename}{$genename};
				}

			}else{ # fallback, if the basename (filename) does not exists, try all 

				foreach my $filename (keys %qdata) { 
					if($geneprintswitch){last;}
					foreach my $key (keys %{$qdata{$filename}}) { 
						my $regexv=$key; 
						my $curLine_test = $curLine;
						
						if(!$doregex && !$exact){$regexv=quotemeta($regexv);}

						my $test_match = 0;
						if( !$exact ){
							
							# use regular expression
							$test_match = $curLine_test =~ $regexv;

						}else{
							
							# directly compare starting with first 5 character of fasta entry as offset
							my $offset = 1; # start at 1 -> fasta entries starts with ">"
							while( !( $test_match = substr($curLine,$offset,length $key) eq $key ) && $offset < 5 ){
								$offset++;
							} 
						}

						if( $test_match ){

							if($qdata{$filename}{$key} eq ""){
								print STDERR "[STDERR] WARNI NG The input ($key) was found multiple times in the fasta files ".(!$exact ? "(maybe try --exact)." : ".")."\n";
							}

							my $headerstr=$curLine;
							if($source){$headerstr=$headerstr." ".$basename;}

							$cur_gene.=$headerstr."\n";	
							# $cur_gene_filename = $qdata{$filename}{$key}.".fasta";
							$cur_gene_filename = $qdata{$filename}{$key};
							$geneprintswitch = 1;
							$genecounterfound++;

							delete $qdata{$filename}{$key};

							last;
						}
					}
				}	
			}

		}else{
			if($geneprintswitch){
				$cur_gene.=$curLine."\n";
			}
		}
	}
	close($FH);
}

if($cur_gene ne ""){
	$cur_gene_filename=~s/\\b//g;
	$cur_gene_filename=~s/[^a-zA-Z0-9.]//g;
	if($tofiles){ # print to files
		my $writemodus=">>";
		if(!exists $cur_gene_firsttime{$cur_gene_filename}){$writemodus=">";$cur_gene_firsttime{$cur_gene_filename}=1;}
		open($FHOUT,$writemodus,$cur_gene_filename);
		print $FHOUT $cur_gene;
	}else{
		print $cur_gene;
	}
	close($FHOUT);
}

if($genecounter != $genecounterfound){
	print STDERR "[STDERR] WARNING The input ($query) contains $genecounter queries, but I extracted $genecounterfound entries out of the fasta(s).";
	if(!$exact){print STDERR " If this is not desired, please consider using the -exact option";}elsif($genecounter > $genecounterfound){print STDERR "\n-> This should not have happen, maybe some fasta files are missing as input?\n(If you cannot solve this error, please send a report to incoming+paulklemm-phd-proteinortho-7278443-issue-\@incoming.gitlab.com or visit https://gitlab.com/paulklemm_PHD/proteinortho/wikis/Error%20Codes for more help. Further more all mails to lechner\@staff.uni-marburg.de are welcome)\n";}
	print STDERR "\n";

	if($genecounterfound < $genecounter){
		print STDERR "The following ids were not found:\n";
		my $counter=0;
		foreach my $filename (keys %qdata) { 
			foreach my $key (keys %{$qdata{$filename}}) { 
				if(10 < $counter++){last}
				print STDERR $key."\n"; 
			}
			if(10 < $counter){print STDERR " ...\n";last}
		}
		print STDERR "\nPlease make sure that the upper ids are part of the given fasta files (try searching these in the given fasta files) !\n";
	}
	
}else{
	print STDERR "[STDERR] All entries of the query are found in the fasta(s).\n";
}
