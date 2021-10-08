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

process af_1 {

	publishDir "${results_dir}",mode:"copy"

	input:
	file vcf_file
	file pop1

	output:
	file "*.frq"

	"""
	vcftools --vcf ${vcf_file} --keep ${pop1} --freq --out pop1
	"""

}

process af_2 {

	publishDir "${results_dir}",mode:"copy"

	input:
	file vcf_file
	file pop2

	output:
	file "*.frq"

	"""
	vcftools --vcf ${vcf_file} --keep ${pop2} --freq --out pop2
	"""

}

process af_3 {

	publishDir "${results_dir}",mode:"copy"

	input:
	file vcf_file
	file pop3

	output:
	file "*.frq"

	"""
	vcftools --vcf ${vcf_file} --keep ${pop3} --freq --out pop3
	"""

}

process moving_files {

	publishDir "${results_dir}",mode:"copy"

	input:
	file fst_files_mix
	file script_R

	output:
	file "*.png"

	errorStrategy 'ignore'
	shell:
	'''
	mkdir results
	mv !{fst_files_mix} results
	Rscript --vanilla !{script_R} ./results pbs.png
	'''
}
