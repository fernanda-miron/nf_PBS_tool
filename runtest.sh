path_data="test/data/pop1_out"
path_data_2="test/data/pop1_pop2"
path_data_3="test/data/pop2_out"
output_directory="$(dirname $path_data)/results"

echo -e "======\n Testing NF execution \n======" \
&& rm -rf $output_directory \
&& nextflow run plot-microRNA-targets.nf \
	--path_data $path_data \
	--path_data_2 $path_data_2 \
	--path_data_3 $path_data_3 \
	--output_dir $output_directory \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
