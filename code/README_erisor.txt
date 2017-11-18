/uufs/chpc.utah.edu/common/home/wolf-group1/Eri_sor/DigitalCommons/code


########################### Data Assembly: iPyRAD commands  ########################
Data were run at two separate facilities. Identical barcodes were used in the different library preparations, so libraries must be demultiplexed separately.

# Create 2 new assemblies:
ipyrad -n erisor1
ipyrad -n erisor2          


# After editing "params-erisor1.txt" and "params-erisor2.txt" (provided) run the first step (demultiplex):
ipyrad -p params-erisor1.txt -s 1
ipyrad -p params-erisor2.txt -s 1


# Merge the two assemblies.
ipyrad -m erisor params-erisor1.txt params-erisor2.txt


# Edit “params-erisor.txt” (provided)
# Finish the remaining data assembly steps with the merged data.
ipyrad -p params-erisor.txt -s 234567


# To branch E. shockleyi into a separate assembly:
ipyrad -p params-erisor.txt -b shockleyi shockleyi_samples.txt


# After editing params-shockleyi.txt (provided), run the branch:
ipyrad -p params-shockleyi.txt -s 234567


######################## Data Analysis and Input Files ########################
Descriptions for each code are at the top of the file.

* Nucleotide diversity and heterozygosity: diversity_stats.R
   * -erisor.vcf
   * -stats files: “lowmsl_shock_stats.txt” and “lowmsl_sored_stats.txt”

* Isolation by distance: mantel_test_IBD.R
   * -erisor_shock.str
   * -erisor_shock_coords.csv

* Structure Statistics: structure_stats.R
   * -erisor.geno
   * -ipy_order.csv

* Structure: structure_erisor.ipynb
   * -erisor.str
   * -erisor.snps.map

* Genetic Distances Population Heatmap: population_heatmap.R
   * -erisor.str
   * -popFac.csv