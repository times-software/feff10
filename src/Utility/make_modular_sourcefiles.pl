#!/usr/bin/perl -w
# This script does the "make src" of the FEFF9 distribution.  It reads the dependency input from feff90/src/DEP/module.mk
# and then copies the source code from all listed files into one big source file feff90/mod/Seq/module_tot.f90 .
# Modules are inserted before other code.  INCLUDE statements are replaced by the referenced include-file, but only down to
# one level : that is, nested include statements WILL NOT WORK with this version of the script.  This can be remedied fairly
# simply by putting the printing statements in a subroutine that calls itself when encountering an include statement.

# The program is meant to reside in the feff90/src/Utilities folder,
# and to be called by the Makefile in feff90/src through the command 'make src'.

# INPUT AND USAGE :
# The program takes one optional input argument.  If called as :
# make_modular_sourcefiles mpi    --> it makes mpi-sources by including feff90/src/PAR/parallel.f90
# make_modular_sourcefiles   or   make_modular_sourcefiles seq  --> it makes non-mpi sourcefiles by including feff90/src/PAR/sequential.f90

# Programmed by Kevin Jorissen, kevinjorissen.pdx@gmail.com , June 2009 for FEFF90.

use strict;
#use 5.010;

# Location of the modular source files :
my $outdir = './../mod/Seq';
# Location of the dependency files :
my $depsdir = './DEP';
# List of FEFF modules :
my @progs = qw(atomic dmdw eels ff2x fms genfmt ldos mkgtr opconsat path pot rdinp screen sfconv xsph dym2feffinp compton rhorrp);

# Process input flags :
my $parallel = 0 ; # Make non-mpi source by default
if (@ARGV > 0) {
    $parallel = ($ARGV[0] eq "mpi");
}
if ($parallel) {
    $outdir = './../mod/MPI';
}
print "Making source files for FEFF modules : @progs .\n";

foreach my $i (@progs) {
#    print "now tackling FEFF module $i \n";
    my $depfile = "${depsdir}/$i.mk";
    if (! open DEP,"<","$depfile") {
        die "Cannot open $depfile : $! ";
    }
    my $filecontent = join '', <DEP>;
    close DEP;
    my @deplist = split /=/,"$filecontent";
# deplist now ought to contain 3 fields.  The first is junk.  The second contains all the source files that are not modules.
# The third contains all the source files that are modules.
    my $modulefiles = $deplist[2];
    $modulefiles =~ s%#.*$% %s ; # remove trailing comments.
    $modulefiles =~ s/^[\s]+// ; # no leading whitespaces
    $modulefiles =~ s%[\\]*[\n]*%%g ; #remove continuation characters and newlines
    $modulefiles =~ s/[\s][\s]+/ /g ; # never more than one whitespace between elements
    $modulefiles =~ s/[\s][\s]+$//g ; # no extra whitespace at end.
    my $regularfiles = $deplist[1];
    $regularfiles =~ s%#.*$% %s ;
    $regularfiles =~ s/^[\s]+// ;
    $regularfiles =~ s%[\\\n]% %g ;
    #$regularfiles =~ s%[\\]*[\n]*%%g ;
    $regularfiles =~ s/[\s][\s]+/ /g ;
    $regularfiles =~ s/[\s][\s]+$//g ; # no extra whitespace at end.
# WARNING : THE ABOVE 10 LINES AREN'T VERY ROBUST AND MAY FAIL IF THE FORMATTING OF THE .mk FILES CHANGES !!!
# Now put back in a list.
# The method of making the parallel/serial versions has changed, so that PAR/par.f90 always holds the
# correct source. Josh Kas 08/09
#    my @reglist = split / /,"$regularfiles"."PAR/par.f90";
#    if ($parallel == 1) {
#        @reglist = split / /,"$regularfiles"."PAR/par.f90";       
#    }
my @reglist = split / /,"$regularfiles" ;
my @modlist = split / /,"$modulefiles" ;
    pop(@reglist) ;
# Now start assembling source file.
    if (-e "src.tmp" ) {
        unlink "src.tmp";
    }
    if (! open FILE,">","src.tmp") {
        die "Cannot open src.tmp : $!";
    }
        select FILE;
    foreach (@modlist,@reglist) {
        if (! open SOURCE,"<","./$_") {
            die "Cannot open $_ : $!";
        }
        (my $pathway = "./$_") =~ s%(^.+[\\\/])([\w\s\.]+)$%$1%; #take the path out of the filename
        while (<SOURCE>) {
            if (! /^[\s]*[#!c]+.*/i && /include[\s]+['"`](.+)['"`]/i ) {
                # an " INCLUDE 'file' " line that is not a comment line
                my $incfile = $1; 
		my $orig_incfile = $1;
                if (! ($incfile =~ m%^\.\.[\\\/]% )) {
                    #If include file doesn't start with "../", it's in the same folder as the sourcefile - so fix the path accordingly.                    
                    $incfile = ${pathway}.${incfile};
		}
		else {
		    #If it does, we need to take out the "../"
		    $incfile =~s%^\.\.%\.% ;
                }
                if (! open INCLUDE,"<","$incfile") {
#KJ                        warn "WARNING - $incfile not found -  you will have to supply it at the time of compilation.";
                        print "      include \'${orig_incfile}\'   ! This line written by make_modular_sourcefiles.\n";
                }
                else {
                    print "!!! EXPANDING INCLUDE statement: $incfile \n";
                    while (<INCLUDE>) {

                        if (! /^[\s]*[#!c]+.*/i && /include[\s]+['"`](.+)['"`]/i ) {
                           #another include line
                           my $nested_incfile = $1; 
		           my $orig_nested_incfile = $1;
                           if (! ($nested_incfile =~ m%^\.\.[\\\/]% )) {
                              $nested_incfile = ${pathway}.${nested_incfile};
		           }
		           else {
		              $nested_incfile =~s%^\.\.%\.% ;
                           }
                           if (! open NESTED,"<","$nested_incfile") {
                              print "      include \'${orig_nested_incfile}\'   ! This line written by make_modular_sourcefiles.\n";
                           }
                           else {
                              print "!!! EXPANDING INCLUDE statement: $nested_incfile \n";
                              while (<NESTED>) {
                                 if (! /^[\s]*[#!c]+.*/i && /include[\s]+['"`](.+)['"`]/i ) {
                                    #another include line
                                    # It would be useful to have the code throw an error or warning here.  However it is part of a large automated setup ...
                                 }
                                 print;  # (Treating it as) just a regular line of code
                              }
                              close NESTED;
                              print "!!! END OF : $nested_incfile \n" ;
                           }
                        }
                        else {
                           print;  # Just a regular line of code
                        }
                    }
                    close INCLUDE;
                    print "!!! END OF : $incfile \n" ;
                }
            }
            else {
                print;   # Just a regular line of code
            }
        }
        print "\n";
	close SOURCE;
    }    

    select STDOUT;
    close FILE;
    ## Check for unexpanded INCLUDE statements (due to nesting of include files): - THIS SECTION DOESN'T WORK YET
    #my @problem_lines = `tcgrep -i include "src.tmp"`;
    #foreach (@problem_lines) {
    #    if (! /^[\s]*[#!c]+.*/i && /include[\s]+['"`](.+)['"`]/i ) {  # test for an " INCLUDE 'file' " line that is not a comment line
    #        die "Whoops - nested INCLUDE files in program $i - I can't handle that! Quitting now.";    
    #    }
    #}
# Put the file in its final location.
    if (-e "${outdir}/${i}_tot.f90") {
        unlink "${outdir}/${i}_tot.f90";
    }
    if (! rename "src.tmp", "${outdir}/${i}_tot.f90" ) {
        die "Cannot rename to ${outdir}/${i}_tot.f90 : $!";
    }
#KJ    print "Finished FEFF program $i \n";
}

#KJ print "all done!\n";
