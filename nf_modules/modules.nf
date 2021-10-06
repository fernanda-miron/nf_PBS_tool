#!/usr/bin/env nextflow

/*Modulos de prueba para  nextflow */
results_dir = "./test/data/results"
intermediates_dir = "./test/data/intermediates"

process fst_wrangling {

	publishDir "${results_dir}",mode:"copy"

	input:
	path fst_files
	file r_script

	output:
	file "*.csv"

	"""
	Rscript --vanilla wrangling.R ${fst_files} "pop1vsout.csv"
	"""

}

process wrangling_2 {

	publishDir "${results_dir}",mode:"copy"

	input:
	path fst_files_2
	file r_script

	output:
	file "*.csv"

	"""
	Rscript --vanilla wrangling.R ${fst_files_2} "pop1vspop2.csv"
	"""

}

process wrangling_3 {

	publishDir "${results_dir}",mode:"copy"

	input:
	path fst_files_3
	file r_script

	output:
	file "*.csv"

	"""
	Rscript --vanilla wrangling.R ${fst_files_3} "pop2vsout.csv"
	"""

}
