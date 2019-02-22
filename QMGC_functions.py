##QMGC_main_all

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
from optparse import OptionParser
from capturer import CaptureOutput
#Utilities


def get_files_with_suffix(dir, suffix1='', suffix2='', inside_word='', mindepth=False):
    """Function that will iterate through the files of a directory and will return a list
    of the full pathnames of the files fulfilling the requirements set with the below parameters.
    :param dir: String (Required). The directory of which the files you want to iterate through. Should always finish with an '/'
    :param suffix1: String (Optional). The prefix you want the files of the dir to START with.
    Default is None: the function will not look for files starting with a particular prefix.
    :param suffix2: String (Optional). The suffix you want the files of the dir to END with.
    Default is None: the function will not look for files ending with a particular suffix.
    :param inside_word: String (Optional). A particular string you want to be included within the
    filenames of the dir. Default is None: the function will not look for files having the given inside word.
    :param mindepth: Boolean (Optional). True: Search for files in all subdirectories of the dir.
    False: Only search for files within dir. Default to False.
    :return: A list of all the files meeting the requirements set above.
    """
    if mindepth==False:
        rfiles = [dir + x for x in os.listdir(dir) if
                  x.startswith(suffix1) and inside_word in x and x.endswith(suffix2)]
    if mindepth==True:
        temp = []
        for path, subdirs, files in os.walk(dir, followlinks=True):
            for name in files:
                temp.append(os.path.join(path, name))
        rfiles = [x for x in temp if x.startswith(suffix1) and inside_word in x and x.endswith(suffix2)]
    rfiles.sort()
    return(rfiles)

#Fuction for running linux commands in parallel and capturing the output
def parallel_command(command, n=2, wordir=os.getcwd()+'/', name='log.txt'):
    """
    Function that takes a list of linux commands as its argument, instructs the linux system to run them
    one by one in parallel using the required amount of processors and capturing the output into a log file.
    :param command: List (Required). List of Linux commands that will be running.
    :param n: Integer (Optional). Number of processors will be used. Default to 2.
    :param wordir: String (Optional). Full pathname of the place where the log file will be saved to. Default
    to the current working directory. Should always finish with an '/'.
    :param name: String (Optional). Name of the output log file. Default to log.txt
    """
    processes = set()
    max_processes = n
    with CaptureOutput() as capturer:
        for i in range(0, len(command)):
            processes.add(subprocess.Popen(command[i], shell=True))
            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update([p for p in processes if p.poll() is not None])

        for p in processes:
            if p.poll() is None:
                p.wait()

    text_file = open(wordir + name, "a+")
    text_file.write("\n" + capturer.get_text())
    text_file.close()
    return

# List all the subdirectories of the input dir
def get_sub_dir(main_dir):
    """
    Function that takes a directory name as its input and returns a list with all its subdirectories.
    :param main_dir: String. Full path name of the directory. Should always finish with an '/'
    """
    sub_dir = []
    for dir in os.listdir(main_dir):
        file_path = main_dir + dir + "/"
        sub_dir.append(file_path)
    return sub_dir


# def pair_1_files(sub_dir, suffix):
#     file_R1 = []
#     for i in range(0, len(sub_dir)):
#         file_R1+=get_files_with_suffix(dir=sub_dir[i], suffix2=suffix)
#     return file_R1
#
# if all the data is in main directory run pair_1_files_main
# def pair_1_files_main(main_dir, suffix):
#     file_R1 = []
#     for file in os.listdir(main_dir):
#         if file.endswith(suffix):
#             file_R1.append(main_dir + file)
#             file_R1.sort()
#     return file_R1

# #Check if the number of files in file_R1 and file_R2 are correct
# def check_number_of_files(main_dir, R1s):
#     command = 'find '+main_dir+' -type f | wc -l'
#     p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     output, error = p.communicate()
#     output = int(output)
#     if len(R1s)*2!=output:
#         print("WARNING!!! Error on getting fq files correctly")
#         print("Number of raw fq files %s, Number of extracted R1 files %s" % (output, len(R1s)*2))
#     else:
#         print("R1s were extracted correctly")
#         print("Number of raw R1 files %s, Number of extracted R1 files %s" % (output, len(R1s)*2))

#Extract unique sample names
def split_list_objects(inlist, splitstring):
    """
    Function that takes a list of files' full pathnames and returns their basename split by a given separator.
    :param inlist: List (Required). List of file full pathnames.
    :param splitstring: String (Required). String which will be used as split separator.
    :return: List of the basenames of the given files split by the input separator.
    """
    file_names = [os.path.basename(i) for i in inlist]
    samples = [i.split(splitstring)[0] for i in file_names]
    return samples

# def get_unique_sample_names_EIGC(file_R1):
#     file_names = [os.path.basename(i) for i in file_R1]
#     samples = ['_'.join(i.split('_')[0:4]) for i in file_names]
#     unique_sample_names = list(set(samples))
#     return unique_sample_names

#Run FastQC on the raw fq files
def get_fastqc_command_from_list(file_R1, fastqc_dir, fastqc_before_trim_result_dir):
    """
    Function that takes a list of fastq files filenames as its input and creates a linux command that runs FastQC
    on each of the fastq files.
    :param file_R1: List (Required). List of fastq files full pathnames.
    :param fastqc_dir: String (Required). Full path of FastQC software file. Should always finish with an '/'.
    :param fastqc_before_trim_result_dir: String (Required). Full path of FastQC result directory.
    :return: List of commands to perform FastQC.
    """
    command = []
    for i in range(0, len(file_R1)):
        a = fastqc_dir + ' --outdir=' + fastqc_before_trim_result_dir + ' ' + file_R1[i]
        command.append(a)
    return command

#Run FastQC on the trimmed fq files
def get_fastqc_command_from_dir(input_dir, fastqc_dir, fastqc_after_trim_result_dir, fqsuffix):
    """
    Function that takes a list of directories containing fastq files as its input and creates
    a linux command that runs FastQC on each of the fastq files in these directories.
    :param input_dir: String (Required). Fullpath of the directories containing the fastq files.
    :param fastqc_dir: String (Required). Full path of FastQC software file. Should always finish with an '/'.
    :param fastqc_after_trim_result_dir: String (Required). Full path of FastQC result directory.
    :param fqsuffix: String (Required). The suffix you want the files of the dir to END with.
    :return: List of linux commands to perform FastQC on each of the input fastq files.
    """
    command = []
    file_R1 = get_files_with_suffix(input_dir, suffix2=fqsuffix, mindepth=False)
    for i in range(0, len(file_R1)):
        a = fastqc_dir + ' --outdir=' + fastqc_after_trim_result_dir + ' ' + file_R1[i]
        command.append(a)
        return command

#Run MultiQC on the fastqc output
def get_multiqc_command(multiqc_dir, fastqc_result_dir, outname):
	"""
	Function that takes a list of directories containing fastq files as its input and creates
    a linux command that runs FastQC on each of the fastq files in these directories.
	:param multiqc_dir: String (Required). Full path of MultiQC software file. Should always finish with an '/'.
	:param fastqc_result_dir: String (Required). Full path of FastQC result directory.
	:param outname: Name of the output file.
	:return: List of linux commands to perform multiqc on the found fastqc reports.
	"""
	command = multiqc_dir+' '+fastqc_result_dir+' -n '+fastqc_result_dir+outname
	return [command]


##################################################################
def get_trimgalore_command(file_R1, file_R2, trimgalore_dir, trimgalore_result_dir):
    """
    Function that takes a list of 1st read-pairs full paths and a list of 2nd read-pairs
    full paths in the same order as its input and returns a command to perform trimgalore!
    trimming on them.
    :param file_R1: List (Required). List of R1 fastq files full pathnames.
    :param file_R2: List (Required). List of R2 fastq files full pathnames.
    :param trimgalore_dir: String (Required). Full path of Trimgalore! software file. Should always finish with an '/'.
    :param trimgalore_result_dir: String (Required). Full path of trimming result directory. Should always finish with an '/'.
    :return: List of linux commands to perform trimming on each of the input read pairs fastq files.
    """
    command = []
    for i in range(0, len(file_R1)):
        a = trimgalore_dir + ' --rrbs --paired -phred33 -q 30 --retain_unpaired -o ' + trimgalore_result_dir + ' ' + \
            file_R1[i] + ' ' + file_R2[i]
        command.append(a)
    return command

def check_number_of_trimmed(main_dir, trimgalore_result_dir):
    """
    Function that takes a the strings of path directories as its input and creates
    anc checks if the number of raw fastq files is the same as the number of trimmed
    fastq files and prints the outcome of the checking.
    :param main_dir: String (Required). The main directory of the analysis where the raw fastq files are.
    Should always finish with an '/'
    :param trimgalore_result_dir: String (Required). full path of trimming result directory.
    Should always finish with an '/'
    """
    mains = get_files_with_suffix(main_dir, suffix2='.fastq.gz', mindepth=True)
    trimmed = get_files_with_suffix(trimgalore_result_dir, suffix2='.fq.gz', mindepth=False)
    if len(mains) != len(trimmed)/2:
        print("Number of trim output seems incorrect! Check again!")
    else:
        print("Trimgalore output files are OK!")

def usable_reads(trimgalore_result_dir):
    """Function that takes the path of trimgalore results directory, searches all the trimming reports there
    and create a list of linux command which add the number of usable reads of each sample in the end of each report.
    :param trimgalore_result_dir: String (Required). full path of trimming result directory.
    Should always finish with an '/'
    :return: List of linux commands to perform the usable reads counting and reporting.
    """
    command = []
    reports = get_files_with_suffix(dir=trimgalore_result_dir, suffix2='trimming_report.txt')
    for rep in reports:
        suf = os.path.basename(rep).split('.fastq.gz_trimming_report.txt')[0]
        fs = get_files_with_suffix(dir=trimgalore_result_dir, suffix1=suf, suffix2='.gz')
        comm="A=$(zcat "+fs[0]+" | echo $((`wc -l`/4))); B=$(zcat "+fs[1]+" | echo $((`wc -l`/4))); " \
                                                                          "printf 'usable reads: '$(($A+$B))'\nunpaired_reads: '$((A))'\nval_reads: '$((B))'\n' >>"+rep
        command.append(comm)
    return(command)

def get_merge_unpaired_command(unpaired_trimgalore, unpaired_bismark, suffix_trim, suffix_bis, split_str, output_dir,
                               suffix_out):
    """
    Function that take specific software and directories paths and an output suffix as its input and returns
    a list of linux commands that concatenate the single-end fastq files which did not align in a paired-end
    mode with the fastq files of the same sample which were left unpaired after the trimming procedure.
    :param unpaired_trimgalore: String (Required). The directory where the trimgalore! unpaired fastq files are. Should always finish with an '/'.
    :param unpaired_bismark: String (Required). The directory where the pe alignment unpmapped fastq files are. Should always finish with an '/'.
    :param suffix_trim: String (Required). Suffix of the trimgalore! unpaired fastq files.
    :param suffix_bis: String (Required). Suffix of the pe alignment unmapped fastq files.
    :param split_str: String (Required). Separator between unique_sample_name and suffices in the names of the input files.
    :param output_dir: String (Required). The directory where the output will be written. Should always finish with an '/'.
    :param suffix_out: String (Required). The suffix of the output files.
    :return: List of linux commands to concatenate each pair of the input fastq files.
    """
    command = []
    file_R_trim = get_files_with_suffix(unpaired_trimgalore, suffix2=suffix_trim)
    file_R_bis = get_files_with_suffix(unpaired_bismark, suffix2=suffix_bis)
    for i in range(0, len(file_R_trim)):
        if len(file_R_bis) != len(file_R_trim) or \
                os.path.basename(file_R_trim[i]).split(split_str)[0] != os.path.basename(file_R_bis[i]).split(split_str)[0]:
            raise Exception('Unpaired/Unmapped input Error! The pair of files passed to fq files merging '
                            'were not paired reads! Check file names!')
        a = os.path.basename(file_R_trim[i]).split(split_str)[0]
        global outlog
        try: outlog+=str(file_R_trim[i])+' '+str(file_R_bis[i])+' were merged to create: '+str(output_dir)+str(a)+str(suffix_out)+"\n"
        except NameError: pass
        com = 'cat ' + file_R_trim[i] + ' ' + file_R_bis[i] + \
              ' > ' + output_dir + a + suffix_out
        command.append(com)
    return command


def get_bismark_commands_pe(bismark_dir, bowtie2_dir, bowtie2_ref, samtools_dir, input_dir, temp_dir, suffix_r1, suffix_r2, aligned_pe_dir):
    """
    Function that take specific software and directories paths and an output suffix as its input and returns
    a list of linux commands to perform paired-end Bismark alignment with bowtie2 in directional model
    using the bismark parameters: --multicore 4 --gzip --unmapped for each pair of fastq files.
    :param bismark_dir: String (Required). Full path of Bismark software file. Should always finish with an '/'.
    :param bowtie2_dir: String (Required). Full path of Bowtie2 software file. Should always finish with an '/'.
    :param bowtie2_ref: String (Required). Full path of Bowtie2 bisuplhite converted reference genome.
    :param samtools_dir: String (Required). Full path of Samtools software file. Should always finish with an '/'.
    :param input_dir: String (Required). The directory where the fastq files which are going to be aligned are.
    Should always finish with an '/'
    :param temp_dir: String (Optional). Temporary directory for bismark. Should always finish with an '/'.
    :param suffix_r1: String (Required). Suffix of the R1 files.
    :param suffix_r2: String (Required). Suffix of the R2 files
    :param aligned_pe_dir: Directory where the alignment output will be stored. Should always finish with an '/'
    :return: List of linux commands to perform PE alignment on each of the input read-pairs fastq files.
    """
    command = []
    file_R1 = get_files_with_suffix(input_dir, suffix2=suffix_r1)
    file_R2 = get_files_with_suffix(input_dir, suffix2=suffix_r2)
    for i in range(0, len(file_R1)):
        if file_R1[i].split('_R1_')[0]!=file_R2[i].split('_R2_')[0]:
            raise Exception('Alignment input Error! The pair of files passed to PE alignment '
                            'were not paired reads! Check file names!')
        else:
            a = bismark_dir + ' --multicore 4 --path_to_bowtie ' + bowtie2_dir + ' --samtools_path ' + \
                samtools_dir + ' --gzip --unmapped --temp_dir ' + temp_dir + ' --output_dir ' + \
                aligned_pe_dir + ' --genome_folder ' + bowtie2_ref + ' -1 ' + file_R1[i] + ' -2 ' + \
                file_R2[i]
            command.append(a)
    return command


def get_bismark_commands_se(bismark_dir, bowtie2_dir, bowtie2_ref, samtools_dir, input_dir,
                                        suffix_r, aligned_se_dir, mode=''):
    """
    Function that take specific software and directories paths and an output suffix as its input and returns
    a list of linux commands to perform single-end Bismark alignment with bowtie2 using the bismark parameters: --multicore 4
    --gzip --unmapped for each pair of fastq files.
    :param bismark_dir: String (Required). Full path of Bismark software file.
    :param bowtie2_dir: String (Required). Full path of Bowtie2 software file.
    :param bowtie2_ref: String (Required). Full path of Bowtie2 bisuplhite converted reference genome.
    :param samtools_dir: String (Required). Full path of Samtools software file. Should always finish with an '/'
    :param input_dir: String (Required). The directory where the fastq files which are going to be aligned are. Should always finish with an '/'.
    :param temp_dir: String (Optional). Temporary directory for bismark.
    :param suffix_r1: String (Required). Suffix of the fastq files.
    :param aligned_se_dir: Directory where the alignment output will be stored. Should always finish with an '/'.
    :param mode: String (Required). Can be either '--pbat' or ''. If --pbat the aligment will run in pbat mode.
    If '' the alignment will run in directional mode.
    :return: List of linux commands to perform PE alignment on each of the input read-pairs fastq files.
    """
    if mode=='' or mode=='--directional':
        mode=''
        gzip = ' --gzip'
    elif mode=='--pbat':
        mode='--pbat'
        gzip = ''
    command = []
    file_R1 = get_files_with_suffix(input_dir, suffix2=suffix_r)
    for i in range(0, len(file_R1)):
        a = bismark_dir +' '+mode+' --multicore 4 --path_to_bowtie ' + bowtie2_dir + ' --samtools_path ' + samtools_dir + \
            gzip+' --unmapped --output_dir ' + aligned_se_dir + ' --genome_folder ' + bowtie2_ref + ' ' + \
            file_R1[i]
        command.append(a)
    return command

def get_del_commands(path, suffixd_1='', suffixd_2=''):
    """Function that takes apecific directory path, a prefix and a suffix and returns a list of linux commands
    that delete the files in the given directory which names that start with this prefix and end with this suffix
    :param suffixd_1: String (Optional). The prefix you want the files of the dir to START with.
    Default is None: the function will not look for files starting with a particular prefix.
    :param suffixd_1: String (Optional). The suffix you want the files of the dir to END with.
    Default is None: the function will not look for files ending with a particular suffix.
    :return: List of linux commands to delete each file in the given directory with the given suffices.
    """
    command = []
    file_R1=get_files_with_suffix(path, suffix1=suffixd_1, suffix2=suffixd_2)
    command = ['rm '+x for x in file_R1]
    return(command)

########################################################################################
def get_merging_bam_lanes_command(samtools_func, aligned_dir, merged_dir, unique_sample_names, read, mode='se'):
    """
    Function that take specific software and directories paths and the required mode as its input and returns
    a list of linux commands to merge the bam files of different lanes of the same samples. Paired-end BAM files
    will be merged together in a different file than Read1 Single-end BAM files which will be merged together
    in a different file than Read2 Single-end BAM files which will be merged together.
    :param samtools_func: String (Required). Full path of Samtools software file.
    :param aligned_dir: String (Required). Full path of the directory where the aligned files are located. Should always finish with an '/'.
    :param merged_dir: String (Required). Full path of the directory where the merged aligned files will be stored. Should always finish with an '/'.
    :param unique_sample_names: List (Required). List of unique names of the samples of the analysis.
    :param read: String (Required). The sub-string inside the input alignment files which separates the R1 files from the R2 files. Example read='1_bismark' or read='2_bismark'. Leave read='' for merging paired-end files.
    :param mode: String (Required). 'se' if merging single-end alignment files, 'pe' if mergine paired-end alignment files.
    :return: List of linux commands to merge the input alignment bam files.
    """
    commands = []
    space = ' '
    if read != '': read2 = '_' + read + '_'
    else: read2='_'
    for i in range(0, len(unique_sample_names)):
        a = samtools_func + ' merge -n --threads 4 ' + merged_dir
        sample = get_files_with_suffix(dir=aligned_dir, suffix1=unique_sample_names[i], inside_word=read,
                                       suffix2='.bam')
        a += unique_sample_names[i]+read2+mode+'_merged.bam'
        for j in range(0, len(sample)):
            a += space + sample[j]
        commands.append(a)
    return commands

def merging_lanes_qc(unique_sample_names, unmerged_dir, merged_dir, read_un, read_me, perc_diff=0.005):
    """Function that takes specific directories path and some parameters, checks if the lanes merging performed correctly by comparing the bam files sizes before and after the merging and prints the result of the QC.
    :param unique_sample_names: List (Required). List of unique names of the samples of the analysis.
    :param unmerged_dir: String (Required). Full path of the directory where the unmerged alignment files are located. Should always finish with an '/'.
    :param merged_dir: String (Required). Full path of the directory where the merged alignment files are located. Should always finish with an '/'.
    :param read_un: String (Required). The sub-string inside the input unmerged alignment files which separates the R1 files from the R2 files. Example read='1_bismark' or read='2_bismark'. Leave read='' for merging paired-end files.
    :param read_me: String (Required). The sub-string inside the input merged alignment files which separates the R1 files from the R2 files. Example read='1_bismark' or read='2_bismark'. Leave read='' for merging paired-end files.
    :param perc_diff: float (Required). 0 to 1 number representing the maximum tolerated file size differences between before and after merging.
    :return: Prints a message indicating whether the merging of the BAM files of different lanes went OK or not.
    """
    for i in range(0, len(unique_sample_names)):
        sample_un = get_files_with_suffix(dir=unmerged_dir, suffix1=unique_sample_names[i], inside_word=read_un, suffix2='.bam')
        sample_me = get_files_with_suffix(dir=merged_dir, suffix1=unique_sample_names[i], inside_word=read_me, suffix2='.bam')
        samples = sample_un + sample_me
        siz = [os.path.getsize(x) for x in samples]
        qc = sum(siz[0:-1]) in range(int(siz[-1]-perc_diff*siz[-1]),int((siz[-1]+perc_diff*siz[-1])))
        if qc==True:
            print('Lanes merging QC: Ok!')
        else: raise Exception('Error in bam lanes merging! Problem with sample '+unique_sample_names[i]+' '+read_me)

def get_sorting_command_pe_se(samtools_func, merged_se_dir, merged_pe_dir, sorted_dir, cpus=4):
    """
    Function that take specific software and directories paths as its input and returns
    a list of linux commands to sort the bam files of different lanes of the same samples.
    :param samtools_func: String (Required). Full path of Samtools software file.
    :param merged_se_dir: String (Required). Full path of the directory where the single-end merged alignment files are located. Should always finish with an '/'.
    :param merged_pe_dir: String (Required). Full path of the directory where the paired-end merged alignment files are located. Should always finish with an '/'.
    :param sorted_dir: String (Required). Full path of the directory where the sorted merged alignment files will be stored. Should always finish with an '/'.
    :param cpus: Int (Required). How many cores will be used for the sorting process.
    :return: List of linux commands to sort the input alignment bam files.
    """
    files1 = get_files_with_suffix(merged_se_dir, suffix2='.bam')
    files2 = get_files_with_suffix(merged_pe_dir, suffix2='.bam')
    files=files1+files2
    commands = []
    for i in range(0, len(files)):
        comm = samtools_func + ' sort -n --threads '+cpus+' '+ files[i] + ' -o ' + \
               sorted_dir+os.path.basename(files[i])[0:-4] + "_sorted.bam"
        commands.append(comm)p
    return commands

def get_mExtraction(bismark_methCall_func, samtools_dir, sorted_dir, meth_dir):
    """
    Function that take specific software and directories paths as its input and returns
    a list of linux commands to perform the methylation calling for each sorted and merged aligned file.
    :param bismark_methCall_func: String (Required). Full path of bismark_methylation_extractor software file.
    :param samtools_dir: String (Required). Full path of Samtools software file.
    :param sorted_dir: String (Required). Full path of the directory where the sorted-merged alignment files are located. Should always finish with an '/'.
    :param meth_dir: String (Required). Full path of the directory where the methylation calling files will be stored. Should always finish with an '/'.
    :return: List of linux commands to perform the methylation calling on the input aligment files.
    """
    commands = []
    files = get_files_with_suffix(sorted_dir, suffix2='.bam')
    for i in range(0, len(files)):
        if '_pe_' in os.path.basename(files[i]): mode=' -p'
        elif '_se_' in os.path.basename(files[i]): mode=' -s'
        comm = bismark_methCall_func + mode+' --gzip --multicore 4 --buffer_size 6G --samtools_path ' + samtools_dir + \
               '  --o ' + meth_dir + ' ' + files[i]
        commands.append(comm)
    return commands

def get_b2bgraph_command(sample_names, bismark_b2bgraph_func, meth_dir, bedgraph_dir, mode='CpG'):
    """
    Function that take specific software and directories paths as its input and returns
    a list of linux commands to merge the methylation calling for each sample (each sample has three final methylation
    calling files: One for Read1 SE aligment, one for Read2 SE aligment and one for the PE aligment) and write them in a bedgraph
    file for each sample using the bismark2bedGraph function with parameters --counts --buffer_size 6G.
    :param sample_names: List (Required). List of unique names of the samples of the analysis.
    :param bismark_b2bgraph_func: String (Required). Full path of bismark2bedGraph software file.
    :param meth_dir: String (Required). Full path of the directory where the methylation calling files are located. Should always finish with an '/'.
    :param bedgraph_dir: String (Required). Full path of the directory where the bedgraph files will be stored. Should always finish with an '/'.
    :param mode: String (Required). 'CpG' to get CpG context calls. Everything elee will enable the --CX option of the bismark2bedGraph software.
    :return: List of linux commands to perform the merging and conversion of methylation calls to a bedgraph output.
    """
    commands = []
    for name in sample_names:
        files = get_files_with_suffix(meth_dir, suffix1=mode, suffix2='txt.gz', inside_word=name)
        a = 'CpG_' + name
        if mode=='CpG': b=' '
        else: b=' --CX '
        comm = bismark_b2bgraph_func + b +'--counts --buffer_size 6G -o ' + a + ' --dir ' + bedgraph_dir + ' ' + ' '.join(files)
        commands.append(comm)
    return (commands)

def get_cov2cyt_command(bismark_cov2cyt_func, bowtie2_ref, bedgraph_dir):
    """
    Function that take specific software and directories paths as its input and returns
    a list of linux commands to convert bedgraph files to cytosine-wide reports for each sample.
    :param bismark_cov2cyt_func: String (Required). Full path of bismark2bedGraph software file.
    :param bowtie2_ref: String (Required). Full path of Bowtie2 bisuplhite converted reference genome.
    :param bedgraph_dir: String (Required). Full path of the directory where the bedgraph files are located. Should always finish with an '/'.
    :return: List of linux commands to perform the conversion of the bedgraph file to a cytosine-wide file for each sample.
    """
    commands=[]
    files=get_files_with_suffix(bedgraph_dir, suffix2='.bismark.cov.gz')
    for i in range(0,len(files)):
        a = os.path.basename(files[i]).split('.gz.bismark.cov.gz')[0]
        comm = bismark_cov2cyt_func + ' -o '+a+' --dir '+bedgraph_dir+' --genome_folder '+bowtie2_ref+' --gzip '+files[i]
        commands.append(comm)
    return(commands)

def get_aligned_counts_command(samtools_func, sorted_dir):
    """
    Function that take specific software and directories paths as its input and returns
    a list of linux commands that count how many aligned reads there are in the alignment files of the given directory.
    :param bismark_cov2cyt_func: String (Required). Full path of bismark2bedGraph software file.
    :param bowtie2_ref: String (Required). Full path of Bowtie2 bisuplhite converted reference genome.
    :param bedgraph_dir: String (Required). Full path of the directory where the bedgraph files are located. Should always finish with an '/'.
    :return: List of linux commands to perform the conversion of the bedgraph file to a cytosine-wide file for each sample.
    """
    files=get_files_with_suffix(sorted_dir, suffix2='.bam')
    commands = ['A=$('+samtools_func+' flagstat -@ 2 '+x+" | head -n 1 | cut -d'+' -f1); " \
                "B=$(basename "+os.path.basename(x)+"); echo $B $A+'\n' >> "+sorted_dir+"count_aligned_reads_QC.txt" for x in files]
    return(commands)

def get_bam2cram_command(samtools_func, bowtie2_ref, sorted_dir):
    """
    Function that take specific software and directories paths as its input and returns
    a list of linux commands that create cram files for each identified bam file in the given directory..
    :param samtools_func: String (Required). Full path of Samtools software file.
    :param bowtie2_ref: String (Required). Full path of Bowtie2 bisuplhite converted reference genome.
    :param sorted_dir: String (Required). Full path of the directory where the bam files are located. Should always finish with an '/'.
    :return: List of linux commands to perform the conversion of the bam files to cram.
    """
    files=get_files_with_suffix(sorted_dir, suffix2='.bam')
    commands = [samtools_func + ' view -@2 -T '+bowtie2_ref+' -C -o '+sorted_dir+os.path.basename(x).split('.bam')[0]+'.cram '+x for x in files]
    return(commands)

def get_final_report_command(bedgraph_dir):
    """
    Function that take a directory with cytosine-wide coverage paths as its input and returns
    a list of linux commands that create a new cpg coverage file for each located file
    containing only Cs covered by at least one read.
    :param bedgraph_dir: String (Required). Full path of the directory where the bedgraph files are located. Should always finish with an '/'.
    :return: List of linux commands to perform the creation of new files.
    """
    files=get_files_with_suffix(bedgraph_dir, suffix2='.CpG_report.txt.gz')

    commands = ["gzip -dc "+x+" |  awk '{if ( ($4 + $5) > 0) {print} }' | gzip > "+bedgraph_dir+os.path.basename(x).split('.CpG_report.txt.gz')[0]+'.CpG_report_v2.txt.gz' for x in files]
    return(commands)


def final_report(unique_sample_names, all_rs, main_result_dir, sorted_dir, aligned_pe_dir, aligned_se_dir, trimgalore_result_dir, bedgraph_cpg_dir):
    """
    Function that creates two reports summarising the metrics of all samples:
    *Per_lane_stats.txt summarises information of the analysis for each pair of fastq files used Column names (2 fastq files per Lane)
     Column names: run_name, lane, usable_reads, unpaired_reads, val_reads (no. of trimgalore passed paired reads), pe_pairs_analysed,
     se_pairs_aligned, se1_reads_analysed, se2_reads_aligned, paired-end mapping efficiency. se1 mapping efficiency, se2 mapping efficiency,
     total_mapping efficiency, qc1, qc2, qc3
     Metrics of the per_lane_Stats are merged to create total metrics for all the samples of the analysis (and not specific lanes)
    *Per_sample_stats.txt summarises information for all the samples in the analysis (in our analysis: 8 fastq files as input - 4 sequencing lanes - 1 sample)
     Column names: run_name, sample_name, qc1, qc2, qc3, usable_reads, mapping_efficiency, identified_cpgs and location of cpgs
    :param unique_sample_names: List (Required). List of unique names of the samples of the analysis.
    :param main_result_dir: String (Required). Full path of the main result directory of the analysis. Should always finish with an '/'.
    :param sorted_dir: String (Required). Full path of the directory where the bam files are located. Should always finish with an '/'.
    :param aligned_pe_dir: Directory where the alignment output will be stored. Should always finish with an '/'.
    :param aligned_se_dir: Directory where the alignment output will be stored. Should always finish with an '/'.
    :param trimgalore_result_dir: String (Required). Full path of trimming result directory. Should always finish with an '/'.
    :return: The function will save the generated reports in the main result directory.
    """
    sheaders = [['run', 'sample', 'gc1', 'qc2', 'qc3', 'usable_reads', 'total_mapping_eff', 'cpgs_identified',
                 'cpg_report_location']]
    with open(main_result_dir + 'Per_sample_stats.txt', "a+", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerows(sheaders)
        fp.close()
    lheaders = [['run', 'lane', 'usable_reads', 'unpaired_reads', 'val_reads', 'pe_pairs_analysed', 'pe_pairs_aligned',
                 'se1_analysed', 'se1_aligned', 'se2_analysed', 'se2_aligned', 'pe_me', 'se1_me', 'se2_me', 'total_me',
                 'qc1', 'qc2', 'qc3']]
    with open(main_result_dir + 'Per_lane_stats.txt', "a+", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerows(lheaders)
        fp.close()
    for un in unique_sample_names:
        namesL = [x.split('/')[-1].split('_R1_')[0] for x in all_rs if '_R1_' in x]
        bam_counts = open(sorted_dir + 'count_aligned_reads_QC.txt', 'r').read().splitlines()
        bam_counts_un = [x for x in bam_counts if x.startswith(un)]
        bam_counts_un.sort()
        bam_counts_unint = [int(re.findall(' \d+ ', x)[0].strip()) for x in bam_counts_un if x != '']
        samples_results=[]
        allbams = []
        print('Processing sample: '+un)
        sample_lanes = [x for x in namesL if x.startswith(un)]
        sample_lanes.sort()
        for nameL in sample_lanes:
            lanes_results = []
            print(nameL)
            filenames = []
            counts = []
            trims = get_files_with_suffix(trimgalore_result_dir, suffix1=nameL, suffix2='trimming_report.txt')
            trims.sort()
            pes= get_files_with_suffix(aligned_pe_dir, suffix1=nameL, suffix2='PE_report.txt')
            pes.sort()
            ses = get_files_with_suffix(aligned_se_dir, suffix1=nameL, suffix2='SE_report.txt')
            ses.sort()
            comms = ["grep 'usable reads' "+trims[0], "grep 'usable reads' "+trims[1], "grep 'unpaired_reads' "+trims[0],
                     "grep 'unpaired_reads' "+trims[1], "grep 'val_reads' "+trims[0], "grep 'val_reads' "+trims[1],
                     "grep 'Sequence pairs analysed in total:' "+ pes[0], "grep 'Number of paired-end alignments with a unique best hit' "+ pes[0],
                     "grep 'Sequences analysed in total' "+ses[0], "grep 'Number of alignments with a unique best hit from the different alignments:' "+ses[0],
                     "grep 'Sequences analysed in total' " + ses[1], "grep 'Number of alignments with a unique best hit from the different alignments:' " + ses[1]]
            for comm in comms:
                p = subprocess.Popen(comm, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                output, error = p.communicate()
                output = re.findall('\d+', str(output))
                if len(output) > 1:
                    print('Error many files same sample and lane!')
                lanes_results.append(int(output[0]))
            lanes_results = [sum(lanes_results[0:2]), sum(lanes_results[2:4]), sum(lanes_results[4:6])] + lanes_results[6:]
            qc1 = lanes_results[0] in range((lanes_results[1]+lanes_results[2])-1000, (lanes_results[1]+lanes_results[2])+1000)
            qc2 = lanes_results[2]/2 in range(lanes_results[3]-1000,lanes_results[3]+1000)
            qc3 = ((lanes_results[3]-lanes_results[4])*2)+lanes_results[1] in range(lanes_results[5]+lanes_results[7]-1000,lanes_results[5]+lanes_results[7]+1000)
            # lanes_results[0] = lanes_results[0]*2
            # lanes_results[1] = lanes_results[1]*2
            unpaired = lanes_results[2]
            pe_me = (lanes_results[4]/lanes_results[3])*100
            se1_me = (lanes_results[6]/lanes_results[5])*100
            se2_me = (lanes_results[8] / lanes_results[7]) * 100
            total_me = (lanes_results[4]*2+lanes_results[6]+lanes_results[8])/lanes_results[0]*100
            lanes_results.extend([pe_me, se1_me, se2_me, total_me, qc1, qc2, qc3])
            lanes_results = [un, nameL]+lanes_results
            samples_results.append(lanes_results)
        with open(main_result_dir + 'Per_lane_stats.txt', "a+", newline='') as fp:
            wr = csv.writer(fp, dialect='excel')
            wr.writerows(samples_results)
            fp.close()
        count=0
        totals = []
        for res in zip(*samples_results):
            res = list(res)
            if  count in [2,6,8,10]:
                totals.append(sum(res))
                count+=1
            elif count==14:
                totals.append(sum(res)/len(res))
                count+=1
            else:
                count+=1
        totals[1] = totals[1] * 2
        qcb1 = totals[2] == bam_counts_unint[0]
        qcb2 = totals[3] == bam_counts_unint[1]
        qcb3 = totals[1] == bam_counts_unint[2]
        covpath = get_files_with_suffix(bedgraph_cpg_dir, suffix2='CpG_report_v2.txt.gz', inside_word=un)[0]
        sep_un = re.findall('([a-zA-Z]\d+)_',un)[0]
        p = subprocess.Popen('less '+covpath+' | wc -l', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = p.communicate()
        cpgs = int(output)
        totals = [[un] + [sep_un] + [qcb1, qcb2, qcb3] + [totals[0], totals[4]] + [cpgs] + [covpath]]
        with open(main_result_dir + 'Per_sample_stats.txt', "a+", newline='') as fp:
            wr = csv.writer(fp, dialect='excel')
            wr.writerows(totals)
            fp.close()
    return()

def send_email(from_email, to_email, password, attach):
    COMMASPACE = ', '
    sender = from_email
    gmail_password = password
    recipients = [to_email]
    # Create the enclosing (outer) message
    outer = MIMEMultipart()
    outer['Subject'] = 'Run has finished'
    outer['To'] = COMMASPACE.join(recipients)
    outer['From'] = sender
    outer.preamble = 'You will not see this in a MIME-aware mail reader.\n'
    # List of attachments
    attachments = [attach]
    # Add the attachments to the message
    for file in attachments:
        try:
            with open(file, 'rb') as fp:
                msg = MIMEBase('application', "octet-stream")
                msg.set_payload(fp.read())
            encoders.encode_base64(msg)
            msg.add_header('Content-Disposition', 'attachment', filename=os.path.basename(file))
            outer.attach(msg)
        except:
            print("Unable to open one of the attachments. Error: ", sys.exc_info()[0])
            raise
    composed = outer.as_string()
    # Send the email
    try:
        with smtplib.SMTP('smtp-mail.outlook.com', 587) as s:
            s.ehlo()
            s.starttls()
            s.ehlo()
            s.login(sender, gmail_password)
            s.sendmail(sender, recipients, composed)
            s.close()
        print("Email sent!")
    except:
        print("Unable to send the email. Error: ", sys.exc_info()[0])
        raise


