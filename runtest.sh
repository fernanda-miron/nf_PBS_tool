path_data="test/data/pop1_out"
path_data_2="test/data/pop1_pop2"
path_data_3="test/data/pop2_out"
vcf_file="test/data/AF_calculation/out.recode.vcf"
pop1="test/data/AF_calculation/pop1"
pop2="test/data/AF_calculation/pop2"
pop3="test/data/AF_calculation/pop_out"
weir_files="test/data/weir_fst/"
output_directory="$(dirname $path_data)/results"

echo -e "======\n Testing NF execution \n======" \
&& rm -rf $output_directory \
&& nextflow run plot-microRNA-targets.nf \
	--path_data $path_data \
	--path_data_2 $path_data_2 \
	--path_data_3 $path_data_3 \
	--vcf_file $vcf_file \
	--pop1 $pop1 \
	--pop2 $pop2 \
	--pop3 $pop3 \
	--weir_files $weir_files \
	--output_dir $output_directory \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
