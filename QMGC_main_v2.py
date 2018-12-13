import sys
import os
import multiprocessing
import subprocess
import time
import gzip
import io
import smtplib
import re
import csv
import datetime
import smtplib
import getpass
from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from QMGC_functions import *
from optparse import OptionParser
from dependencies import *

#Information for sending the notification e-mail:
myemail='digrigor@hotmail.com'
toemail='d.grigoriadis@qmul.ac.uk'
password = '!Diosmailprotected!'
#password = getpass.getpass('Please insert your '+myemail+' email account password: ')

global outlog
outlog = '-------------\n'
outlog+= 'RUN SUMMARY\n'
outlog+= '-------------\n'
outlog+= '\n'
outlog+= 'User: Dionysios Grigoriadis (Grigor02)(dio01)\n'
outlog+= 'Time: '+str(datetime.datetime.now())+'\n'

start = time.time()

# parsing the command line arguments

parser = OptionParser()
parser.add_option("-i", "--Input_directory", dest="inputdir",
				  help="Main directory where data will be read from. This directory will have subdirectories (depth 1) that contains fastaq files")
parser.add_option("-r", "--Result_directory", dest="resultdir",
				  help="Main Result directory, where the result will be written to.")

parser.add_option("-1", "--Pair_one", dest="pair1",
				  help="suffix of paired one data")

parser.add_option("-2", "--Pair_two", dest="pair2",
				  help="suffix of paired two data")

parser.add_option("-p", "--proc",
				  dest="no_of_prcessor", default=1,
				  help="Initialize the number of processor you want to use. Default is 1")

(options, args) = parser.parse_args()

# Running the command

main_dir = options.inputdir
main_result_dir = options.resultdir
max_process = int(options.no_of_prcessor)

p1_suffix = options.pair1
p2_suffix = options.pair2

outlog+='Data input directory: '+main_dir+'\n'
outlog+='Result directory: '+main_result_dir+'\n'

print(main_dir)
print(main_result_dir)
max_process = 5

exec(open('directories_handling.py').read())
#Extract the fastq filenames
file_R1 = get_files_with_suffix(main_dir, suffix2=p1_suffix, mindepth=True)
file_R2 = get_files_with_suffix(main_dir, suffix2=p2_suffix, mindepth=True)

#Perform some QC
print('------------------------------------------')
print('RUN INITIATION:')
print(len(file_R1),'R1 files were detected')
print(len(file_R2),'R2 files were detected')
print('Total samples =',(len(file_R1)+len(file_R2))/8)
outlog+='Total samples ='+str((len(file_R1)+len(file_R2))/8)+'\n'
print('-')
print("Are R1 and R2 set in the same order?")
print("Order check:")
print([x.split('_R')[0] for x in file_R1] == [x.split('_R')[0] for x in file_R2])
outlog+='Order check: '+str([x.split('_R')[0] for x in file_R1] == [x.split('_R')[0] for x in file_R2])+'\n'
print('Number of R1s and R2s check:')
print(list(set(split_list_objects(file_R1, '_L00'))) == list(set(split_list_objects(file_R2, '_L00'))))
outlog+='Number of R1s and R2s check: '+str((list(set(split_list_objects(file_R1, '_L00'))) ==
											 list(set(split_list_objects(file_R2, '_L00')))))+'\n'
print('Number of samples check:')
print(len(list(set(split_list_objects(file_R1, '_L00'))))==(len(file_R1)+len(file_R2))/8)
outlog+='Number of samples check: '+str(len(list(set(split_list_objects(file_R1, '_L00'))))==
										(len(file_R1)+len(file_R2))/8)+'\n'
print('------------------------------------------')

#Get the unique sample names
all_rs = file_R1 + file_R2
unique_sample_names = list(set(split_list_objects(file_R1, '_L00')))
unique_sample_names_backup = list(set(split_list_objects(file_R1, '_L00')))

##FastQC Before trimming
outlog+='\n'
log = 'FastQC_before_trim_log.txt'
command_fastqc_bt = get_fastqc_command_from_list(file_R1+file_R2, fastqc_dir, fastqc_before_trim_result_dir)
parallel_command(command_fastqc_bt, 20, fastqc_before_trim_result_dir, log)
outlog+='FASTQC Pre-trimming: Done!\n'


##TRIMGALORE
log = "trimgalore_log.txt"
command_tgalore = get_trimgalore_command(file_R1, file_R2, trimgalore_dir, trimgalore_result_dir)
parallel_command(command_tgalore, 15, trimgalore_result_dir, log)


##Add the number of usable, unpaired and val reads to trimgalore reports
log = 'usable_reads_log.txt'
get_usable_reads = usable_reads(trimgalore_result_dir)
parallel_command(get_usable_reads, 40, trimgalore_result_dir, log)
outlog+='Trimgalore!: Done!\n'

#QC
check_number_of_trimmed(main_dir, trimgalore_result_dir)

##FastQC After trimming
log = 'FastQC_after_trim_log.txt'
fqsuffix = '.fq.gz'
command_fastqc_at = get_fastqc_command_from_dir(trimgalore_result_dir, fastqc_dir, fastqc_after_trim_result_dir, fqsuffix)
parallel_command(command_fastqc_at, 20, fastqc_after_trim_result_dir, log)
outlog+='FASTQC After-trimming: Done!\n'

##MULTIQC
log = 'MultiQC_log.txt'
command_mqc_before = get_multiqc_command(multiqc_dir, fastqc_result_dir, 'MultiQC_report')
parallel_command(command_mqc_before, 1, fastqc_result_dir, log)
outlog+='MULTIQC: Done!\n'

##PE Alignment - BISMARK PE directional
log = "bismark_PE_log.txt"
suffix1 = "val_1.fq.gz"
suffix2 = "val_2.fq.gz"
aligning_command_pe = get_bismark_commands_pe(bismark_dir, bowtie2_dir, bowtie2_ref, samtools_dir, trimgalore_result_dir, temp_dir, suffix1, suffix2, aligned_pe_dir)
parallel_command(aligning_command_pe, 3, aligned_pe_dir, log)
outlog+='PE QC & Alignment: Done!\n'

## Merging unpaired R1 SE after Trimgalore with unmapped R1 SE after Bismark
log = "merging_bismark_and_trimgalore_R1_SE_unpaired_log.txt"
suffix_trim = "unpaired_1.fq.gz"
suffix_bis = "unmapped_reads_1.fq.gz"
split_str = "_R1_"
suffix_out = "_merged_unpaired_1.fq.gz"
merge_unpaired_R1 = get_merge_unpaired_command(trimgalore_result_dir, aligned_pe_dir, suffix_trim, suffix_bis, split_str, merged_all, suffix_out)
parallel_command(merge_unpaired_R1, 4, merged_all, log)
outlog+='SE1 merge unpaired-unmapped: Done!\n'

## Merging unpaired R2 SE after Trimgalore with unmapped R2 SE after Bismark
log = "merging_bismark_and_trimgalore_R2_SE_unpaired_log.txt"
suffix_trim = "unpaired_2.fq.gz"
suffix_bis = "unmapped_reads_2.fq.gz"
split_str = "_R2_"
suffix_out = "_merged_unpaired_2.fq.gz"
merge_unpaired_R2 = get_merge_unpaired_command(trimgalore_result_dir, aligned_pe_dir, suffix_trim, suffix_bis, split_str, merged_all, suffix_out)
parallel_command(merge_unpaired_R2, 4, merged_all, log)
outlog+='SE2 merge unpaired-unmapped: Done!\n'

##Deleting Trimgalore fq.gz files
log="Delete_trimgalore_log.txt"
del_trim_commands = get_del_commands(trimgalore_result_dir, suffixd_2='.fq.gz')
parallel_command(del_trim_commands, 1, trimgalore_result_dir, log)
outlog+='Trimgalore files: Deleted!\n'

## Aligning R1 Single end reads in directional mode
log = "align_R1_SE_log.txt"
suffix1 = "1.fq.gz"
aligning_command_se_r1 = get_bismark_commands_se(bismark_dir, bowtie2_dir, bowtie2_ref, samtools_dir, merged_all, suffix1, aligned_se_dir, mode='')
parallel_command(aligning_command_se_r1, 4, aligned_se_dir, log)
outlog+='SE1 QC & Alignment: Done!\n'

## Aligning R2 Single end reads in pbat mode
log = "align_R2_SE_log.txt"
suffix1 = "2.fq.gz"
aligning_command_se_r2 = get_bismark_commands_se(bismark_dir, bowtie2_dir, bowtie2_ref, samtools_dir, merged_all, suffix1, aligned_se_dir, mode='--pbat')
parallel_command(aligning_command_se_r2, 4, aligned_se_dir, log)
outlog+='SE2 QC & Alignment: Done!\n'

##Merging lanes (bam files) for SE R1
log = "merge_aligned_R1_log.txt"
merging_lanes_se1_command = get_merging_bam_lanes_command(samtools_func, aligned_se_dir, alinged_merged_se_dir, unique_sample_names, read='1_bismark', mode='se')
parallel_command(merging_lanes_se1_command, 20, alinged_merged_se_dir, log)

#QC SE R1 lanes merging
merging_lanes_qc(unique_sample_names, aligned_se_dir, alinged_merged_se_dir, read_un='1_bismark', read_me='1_bismark_se')
outlog+='Merging SE1 BAM files & Performing QC: Done!\n'

##Merging lanes (bam files) for SE R2
log = "merge_aligned_R2_log.txt"
merging_lanes_se2_command = get_merging_bam_lanes_command(samtools_func, aligned_se_dir, alinged_merged_se_dir, unique_sample_names, read='2_bismark', mode='se')
parallel_command(merging_lanes_se2_command, 20, alinged_merged_se_dir, log)

#QC SE R2 lanes merging
merging_lanes_qc(unique_sample_names, aligned_se_dir, alinged_merged_se_dir, read_un='2_bismark', read_me='2_bismark_se')
outlog+='Merging SE2 BAM files & Performing QC: Done!\n'

##Merging lanes (bam files) for PE
log = "merge_aligned_PE_log.txt"
merging_lanes_pe_command = get_merging_bam_lanes_command(samtools_func, aligned_pe_dir, alinged_merged_pe_dir, unique_sample_names, read='', mode='pe')
parallel_command(merging_lanes_pe_command, 20, alinged_merged_pe_dir, log)

#QC PE lanes merging
merging_lanes_qc(unique_sample_names, aligned_pe_dir, alinged_merged_pe_dir, read_un='', read_me='')
outlog+='Merging PE BAM files & Performing QC: Done!\n'

##Deleting lane-specific bam files and no-longer needed fq files
log="Delete_lanes_bam_files_se_log.txt"
del_lanebam_commands = get_del_commands(aligned_se_dir, suffixd_2='.bam')
parallel_command(del_lanebam_commands, 1, aligned_se_dir, log)
log="Delete_lanes_bam_files_pe_log.txt"
del_lanebam_commands = get_del_commands(aligned_pe_dir, suffixd_2='.bam')
parallel_command(del_lanebam_commands, 1, aligned_pe_dir, log)
log="Delete_pe_fq_files_log.txt"
del_unmappedpe_commands = get_del_commands(aligned_pe_dir, suffixd_2='.gz')
parallel_command(del_unmappedpe_commands, 1, aligned_pe_dir, log)
outlog+='Delete Lane-specific bam files\n'

##Sorting merged bam files
log = "sort_merged_bam_files_log.txt"
command_sorting = get_sorting_command_pe_se(samtools_func, alinged_merged_se_dir, alinged_merged_pe_dir, sorted_dir)
parallel_command(command_sorting, 20, sorted_dir, log)
outlog+='Sort merged bam files: Done!\n'

#Merged to sorted bam files QC
samples_mer = get_files_with_suffix(dir=alinged_merged_pe_dir, suffix2='.bam')+get_files_with_suffix(dir=alinged_merged_se_dir, suffix2='.bam')
for x in samples_mer:
	sor = [sorted_dir+j for j in os.listdir(sorted_dir) if j.startswith(os.path.basename(x).split('_merged.bam')[0])][0]
	siz1 = os.path.getsize(x)
	siz2 = os.path.getsize(sor)
	qc = siz2/siz1
	if qc >= 0.94: print('Ok!')
	else: raise Exception('Error in bam sorting! Problem with sample: '+sor)
outlog+='Merged to sorted bam files QC: Done!\n'

##Count aligned reads in bam files for QC
log = "count_aligned_reads_QC_log.txt"
count_aligned_in_bams = get_aligned_counts_command(samtools_func, sorted_dir)
parallel_command(count_aligned_in_bams, 20, sorted_dir, log)
outlog+='Count reads in sorted bam files: Done!\n'

##Delete merged bam files
log="Delete_lanes_bam_files_se_log.txt"
del_mergedbam_commands = get_del_commands(alinged_merged_pe_dir, suffixd_2='.bam') + get_del_commands(alinged_merged_se_dir, suffixd_2='.bam')
parallel_command(del_mergedbam_commands, 1, aligned_se_dir, log)
outlog+='Delete merged bam files: Done!\n'

##Methylation Extraction for sorted bam files
log = "methcal_for_pe_log.txt"
command_mextraction = get_mExtraction(bismark_methCall_func, samtools_dir, sorted_dir, methcal_dir)
parallel_command(command_mextraction, 15, methcal_dir, log)
outlog+='Methylation Extraction: Done!\n'

##Compress sorted bam files to cram files
log='bam2cram_log.txt'
bam2cram_com = get_bam2cram_command(samtools_func, combined_reference, sorted_dir)
parallel_command(bam2cram_com, 20, sorted_dir, log)
outlog+='Sorted BAM files to CRAM compression: Done!\n'

##Delete sorted bam files
log="Delete_sorted_bam_files_log.txt"
del_sortedbam_commands = get_del_commands(sorted_dir, suffixd_2='.bam')
parallel_command(del_sortedbam_commands, 1, sorted_dir, log)
outlog+='Delete sorted BAM files: Done!\n'

#Get bedgraph files
#CpG
log = "bedcov_cpg_merge_log.txt"
command_bedcov_merge = get_b2bgraph_command(unique_sample_names, bismark_b2bgraph_func, methcal_dir, bedgraph_cpg_dir, mode='CpG')
parallel_command(command_bedcov_merge, 20, bedgraph_cpg_dir, log)
outlog+='Get bedgraph: Done!\n'

##Delete methylation extraction files
log="Delete_mem_files_log.txt"
del_me_files_commands = get_del_commands(methcal_dir, suffixd_2='.txt.gz')
parallel_command(del_me_files_commands, 1, methcal_dir, log)
outlog+='Delete methylation extraction files: Done!\n'

#Coverage to Genomewide Cytosine Report
log = "cov2cyt_log.txt"
command_cov2cyt = get_cov2cyt_command(bismark_cov2cyt_func, bowtie2_ref, bedgraph_cpg_dir)
parallel_command(command_cov2cyt, 20, bedgraph_cpg_dir, log)
outlog+='Coverage to cytosine step: Done!\n'

#Remove 0-covered CpG sites
log= "get_final_report_log.txt"
freport_command = get_final_report_command(bedgraph_cpg_dir)
parallel_command(freport_command, 20, bedgraph_cpg_dir, log)
outlog+='Remove the empty CpG sites from the report. Done!\n'

##Delete the genome-wide reports
log="Delete_genome-wide_reports_log.txt"
del_me_files_commands = get_del_commands(bedgraph_cpg_dir, suffixd_2='.CpG_report.txt.gz')
parallel_command(del_me_files_commands, 1, methcal_dir, log)
outlog+='Delete methylation extraction files: Done!\n'

outlog+='Total Input Size: '+str(sum([os.path.getsize(x) for x in file_R1+file_R2])/1000000000)+' GB'
print('dio')
end = time.time()
elapsed = end - start
outlog+=str(elapsed)+' sec. total run time'

with open(main_result_dir+'final_running_log.txt', 'w+') as fh:
	fh.write(outlog)
	fh.close()

#Create the summary/metrics files
final_report(unique_sample_names, all_rs, main_result_dir, sorted_dir, aligned_pe_dir, aligned_se_dir, trimgalore_result_dir, bedgraph_cpg_dir)

send_email(myemail, toemail, password, main_result_dir+'final_running_log.txt')