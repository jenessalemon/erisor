Directory for erisor project.
QC- quality control folder:
  Fix_str.py- takes a .str file, and gets rid of the extra columns that interfere with downstream analysis.
  Duffy_QC.R- builds a matrix and pulls out the 5 closest individuals to a replicate sample. (the modified .str is the input)
  erisor_reps_≈9sout.str is the modified .str file (output from Fix_str.py)
  erisor_reps_stats.txt is the past table in the ipyrad output stats, and is used in loci_reads_regression.R
  loci_reads_regression.R reads in the stats file, plots the loci in the assembly for an individual by the raw reads for that individual, and then we regress loci on reads.
gompert- files needed for the gompert pipeline, currently on hold.
new_eri- the folder I run ipyrad from. SLURM, barcodes, and parameters files for the runs I have done. Results and output files can be found in CHPC permanent storage.
  
