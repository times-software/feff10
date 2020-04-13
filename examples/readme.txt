FEFF9.7 EXAMPLES
###################



The "examples" folder contains a collection of reference calculations.
E.g., in the folder "examples/ELNES/Cu" you will find a calculation of the ELNES spectrum of Cu.

Each example contains the following files :
* The input files needed to run the calculation.  Usually, this is only "feff.inp" .  For a few
   calculations, there is also a "loss.dat", a "spring.inp", or a "feff.dym" file.
* A reference spectrum.  Usually, this file is called "referencexmu.dat".  Sometimes, it is called
   "reference_eels.dat".
* An archive of the full reference calculation, called "REFERENCE.zip".  This archive extracts
   to a subdirectory "REFERENCE/" containing all files pertaining to the calculation.
   

You can run any example calculation, either by :
* opening "examples/ELNES/Cu/feff.inp" in the JFEFF GUI; or
* "cd examples/ELNES/Cu/feff.inp ; feff9" (or a similar command)  on the command line (linux/mac)

After running the calculation, you can compare the produced spectrum in "xmu.dat" to the provided
reference in "referencexmu.dat".  For example, in gnuplot:
> p 'xmu.dat' u 1:4 w lp,’referencexmu.dat' u 1:4 w lp
or, for EELS calculations,
> p 'eels.dat' u 1:2 w lp,’reference_eels.dat' u 1:2 w lp
The spectra should be identical or very close.

If you wish to compare intermediate files, you can unzip the REFERENCE.zip archive.
This will give you all the files we produced in a REFERENCE/ subfolder.  It is possible
that there will be small differences in numerical results or file formatting, because we
frequently update the FEFF9 code.

PLEASE KEEP IN MIND:  Many of the calculations here have a low cutoff radius for the FMS or SCF card.  We did that in
order to keep the tests fast.  However for real-world results the calculation needs to be converged
in terms of the SCF and FMS cutoff radio.  This may require using substantially higher and slower
values than those used in these example calculations.  E.g. if a BN calculation contains
SCF 4.0
FMS 4.0
You should run further calculation where you increase one of these parameters until the spectrum doesn’t change anymore.  You keep the lowest value that gives the converged spectrum.  Then you do the same for the other parameter.  Perhaps the final result will be something like
SCF 5.0
FMS 7.2
Typically, the SCF radius should contain between 30 and 70 atoms; the FMS radius between 100 and 300 atoms.  But it depends on the material.
Other changes to feff.inp may be necessary.  But you should at least check SCF and FMS.  For a new material, it can also be a good idea to plot the Density of States by using the LDOS, and inspect the result.

COMPTON calculations and KSPACE calculations can take a very long time.



Do not hesitate to contact us if you require assistance.


Best regards,



The FEFF team.


