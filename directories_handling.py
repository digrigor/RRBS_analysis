aligned_dir = main_result_dir + 'Aligned/'
if not os.path.exists(aligned_dir):
	os.makedirs(aligned_dir)

fastqc_result_dir = main_result_dir + 'FASTQC/'
if not os.path.exists(fastqc_result_dir):
	os.makedirs(fastqc_result_dir)

fastqc_before_trim_result_dir = fastqc_result_dir + 'Before_trim/'
if not os.path.exists(fastqc_before_trim_result_dir):
	os.makedirs(fastqc_before_trim_result_dir)

fastqc_after_trim_result_dir = fastqc_result_dir + 'After_trim/'
if not os.path.exists(fastqc_after_trim_result_dir):
	os.makedirs(fastqc_after_trim_result_dir)

multiqc_result_dir = main_result_dir + 'MULTIQC/'
if not os.path.exists(multiqc_result_dir):
	os.makedirs(multiqc_result_dir)

trimgalore_result_dir = aligned_dir + 'Trimgalore/'
if not os.path.exists(trimgalore_result_dir):
	os.makedirs(trimgalore_result_dir)

temp_dir = main_result_dir + 'temp/'
if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)

aligned_pe_dir = aligned_dir + 'Aligned_pe/'
if not os.path.exists(aligned_pe_dir):
	os.makedirs(aligned_pe_dir)

alinged_merged_pe_dir = aligned_pe_dir + 'Aligned_merged/'
if not os.path.exists(alinged_merged_pe_dir):
	os.makedirs(alinged_merged_pe_dir)

merged_unpaired = main_result_dir + 'Merged_all_unpaired/'
if not os.path.exists(merged_unpaired):
	os.makedirs(merged_unpaired)

merged_all = merged_unpaired + 'Merged_all/'
if not os.path.exists(merged_all):
	os.makedirs(merged_all)

aligned_se_dir = aligned_dir + 'Aligned_se/'
if not os.path.exists(aligned_se_dir):
	os.makedirs(aligned_se_dir)

merged_unpaired_all = main_result_dir + 'Merged_all/'
if not os.path.exists(merged_unpaired_all):
	os.makedirs(merged_unpaired_all)

alinged_merged_se_dir = aligned_se_dir + 'Aligned_merged/'
if not os.path.exists(alinged_merged_se_dir):
	os.makedirs(alinged_merged_se_dir)

merged_pe_se_dir = main_result_dir + 'Merged_PE_SE/'
if not os.path.exists(merged_pe_se_dir):
	os.makedirs(merged_pe_se_dir)

sorted_dir = main_result_dir + 'Sorted_merged_se_pe/'
if not os.path.exists(sorted_dir):
	os.makedirs(sorted_dir)

methcal_dir = main_result_dir +'Meth_call/'
if not os.path.exists(methcal_dir):
    os.makedirs(methcal_dir)

bedgraph_dir = main_result_dir +'Bedgraph/'
if not os.path.exists(bedgraph_dir):
	os.makedirs(bedgraph_dir)

bedgraph_cpg_dir = bedgraph_dir +'CpG/'
if not os.path.exists(bedgraph_cpg_dir):
	os.makedirs(bedgraph_cpg_dir)

bedgraph_chh_dir = bedgraph_dir +'CHH/'
if not os.path.exists(bedgraph_chh_dir):
	os.makedirs(bedgraph_chh_dir)

bedcov_dir = bedgraph_dir + 'Bedcov/'
if not os.path.exists(bedcov_dir):
	os.makedirs(bedcov_dir)