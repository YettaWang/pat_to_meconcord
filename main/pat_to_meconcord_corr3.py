import numpy as np
import pandas as pd
import math
import sys,getopt
import random
from scipy import stats
import multiprocessing
import os
import subprocess

from itertools import product
import argparse

# chrom = 'chr19'
# outpre = 'merge_9_mini'

# 定义命令
def bam_to_pat(bam_file, bam_path, out_file_path):
    bamfile_path = os.path.join(bam_path, bam_file)
    command0 =(f"wgbstools bam2pat --genome mm10 {bamfile_path} {out_file_path}")
    subprocess.run(command0, shell=True, check=True)
    return out_file_path

def process_pat(input_path, cpg_file, output_path):
    command_1 = (
        f"zcat {input_path} | "
        f"awk -F'\t' 'BEGIN {{OFS=\"\\t\"}} {{new_col = $2 + length($3) - 1; print $0 \"\\t\" new_col}}' | "
        f"awk -F'\t' 'NR==FNR {{a[$3]=$2; next}} {{print $0 \"\\t\" a[$2] \"\\t\" (a[$5]+1)}}' {cpg_file} - > {output_path}"
    )
    subprocess.run(command_1, shell=True, check=True)

def bam_to_meconcord(bam, chr, bam_path, cpg_file_path, region_path):
    ##获得对应染色体的bam文件
    bamfile_path = os.path.join(bam_path, bam)
    chr_bam = (f'{chr}_{bam}')
    chr_bamfile_path = os.path.join(bam_path, chr_bam)
    bam_to_chr =(f"samtools view -b {bamfile_path} {chr} > {chr_bamfile_path}")
    subprocess.run(bam_to_chr, shell=True, check=True)
    bam_to_pat(chr_bam, bam_path, bam_path)
    bam_name_pre = chr_bam.split(".")[0]
    pat_name = f"{bam_name_pre}.pat.gz"
    pat_path = os.path.join(bam_path, pat_name)
    pat_process_name = f"{bam_name_pre}_new.pat"
    pat_process_path = os.path.join(bam_path, pat_process_name)
    process_pat(pat_path, cpg_file_path, pat_process_path)

    pat_data = pd.read_csv(pat_process_path,sep = '\t',header = None,)
    CG_bed = pd.read_csv(cpg_file_path,sep = '\t',header = None,)
    region = pd.read_csv(region_path,sep = '\t',header = None,)
    meconcord_name = f"{bam_name_pre}_MeConcord.txt"
    meconcord_path = os.path.join(bam_path, meconcord_name)
    outdata = open(meconcord_path,'w')
    min_cpg = 4
    min_read = 4
    min_overlap = 4
    max_read = 500
    filter_reads_cpgnum = 2
    basedon0 = 0
    max_dis = 600
    reads_chr = pat_data.sort_values(by=[5,6],ascending = [True,True])
    bin_chrom = region.loc[region.iloc[:,0]==chr,:].copy()
    bin_sort = bin_chrom.sort_values(by=[2],ascending = [True])
    interval_size =5000
    i = 0
    idx = 0
    split_num = int(math.ceil((bin_sort.iloc[i,2]-bin_sort.iloc[i,1]+1)*1.0/interval_size))
    if basedon0 == 1:
        interval_name = chr+'_'+str(bin_sort.iloc[i,1]-1)+'_'+str(bin_sort.iloc[i,2])#to 0-based
    else:
        interval_name = chr+'_'+str(bin_sort.iloc[i,1])+'_'+str(bin_sort.iloc[i,2])
    for j in range(0,split_num):
        pos1 = bin_sort.iloc[i,1]+j*interval_size ##区间的位置
        # print(pos1)
        pos2 = min(bin_sort.iloc[i,2], bin_sort.iloc[i,1]+(j+1)*interval_size-1)
        # print(pos2)

        cpg_pos_wanted = CG_bed.loc[(CG_bed.iloc[:,1]>=pos1) & (CG_bed.iloc[:,1]<=pos2),:]
        # print(cpg_pos_wanted)
        cpg_pos_help = CG_bed.loc[(CG_bed.iloc[:,1]>=(pos1-2000)) & (CG_bed.iloc[:,1]<=pos2),:]
        # print(cpg_pos_help)
        if cpg_pos_wanted.shape[0]>=min_cpg:##bin has at least cpgs
            bin_reads_me = np.zeros([0,cpg_pos_wanted.shape[0]],int)
            # print(bin_reads_me)
            bin_reads_unme = np.zeros([0,cpg_pos_wanted.shape[0]],int)
            # print(bin_reads_unme)
            # print(idx,pos1,reads_sort.iloc[idx,3],reads_sort.iloc[idx,2])
            while idx > 0 and pos1 - reads_chr.iloc[idx,6]< max_dis:
                idx -= 1
            while idx < (reads_chr.shape[0]-1) and pos1 - reads_chr.iloc[idx,6] > max_dis:
                idx += 1
            while idx < (reads_chr.shape[0]-1) and reads_chr.iloc[idx,5] < pos2: ##reads的起始小于区间终止
                bin_reads_me_read = np.zeros([1,cpg_pos_wanted.shape[0]],int)
                bin_reads_unme_read = np.zeros([1,cpg_pos_wanted.shape[0]],int)
                ##原来reads的2对应现在reads的5
                if len(reads_chr.iloc[idx,2]) >=1:                        
                    if sum((cpg_pos_wanted.iloc[:,1]>=reads_chr.iloc[idx,5]) & (cpg_pos_wanted.iloc[:,1]<reads_chr.iloc[idx,6])) >= 1:
                        # cell_id = reads_chr.iloc[idx,12]
                        # cell_ids.append(cell_id)
                        # print(len(cell_ids))
                        lost_num = sum((reads_chr.iloc[idx,5]<=cpg_pos_help.iloc[:,1]) & (cpg_pos_wanted.iloc[0,1]>cpg_pos_help.iloc[:,1]))
                        # print('lost_numA:',lost_num)
                        methy_reads = reads_chr.iloc[idx,2][lost_num:]
                        reads_CG = len(methy_reads)
                        reads_start_index = reads_chr.iloc[idx,1] + lost_num #modify
                        bin_start_index = cpg_pos_wanted.iloc[0,2] #modify
                        if lost_num > 0 and reads_CG <= cpg_pos_wanted.shape[0]:
                            bin_CG = np.zeros([1, cpg_pos_wanted.shape[0]], dtype=str)
                            print(f"Shape of bin_CG_b: {bin_CG.shape}")
                            bin_CG[0, :reads_CG] = list(methy_reads)
                            print(f"Shape of bin_CG_f: {bin_CG.shape}")
                            bin_CG[0, reads_CG:] = '0'
                            # bin_CG[0, :reads_CG] = list(methy_reads[:reads_CG]) 
                        if lost_num > 0 and reads_CG > cpg_pos_wanted.shape[0]:
                            bin_CG = np.zeros([1, cpg_pos_wanted.shape[0]], dtype=str)
                            bin_CG[0, :cpg_pos_wanted.shape[0]] = list(methy_reads[:cpg_pos_wanted.shape[0]])
                        # if lost_num == 0 and reads_start_index - bin_start_index + reads_CG <= cpg_pos_wanted.shape[0]:
                        if lost_num == 0:
                            index_minus = reads_start_index - bin_start_index
                            bin_CG = np.zeros([1, cpg_pos_wanted.shape[0]], dtype=str)
                            # print(bin_CG,methy_reads,index_minus)
                            # print(bin_CG[0, index_minus:],list(methy_reads))
                            available_space = cpg_pos_wanted.shape[0] - index_minus
            
                            # Ensure methy_reads fits into the available space in bin_CG
                            if reads_CG <= available_space:
                                bin_CG[0, index_minus:index_minus + reads_CG] = list(methy_reads)
                            else:
                                bin_CG[0, index_minus:index_minus + available_space] = list(methy_reads[:available_space])
                            # bin_CG[0, index_minus:] = list(methy_reads)
                        # if lost_num == 0 and reads_start_index - bin_start_index + reads_CG > cpg_pos_wanted.shape[0]:
                        #     index_minus = reads_start_index - bin_start_index
                        #     bin_CG = np.zeros([1, cpg_pos_wanted.shape[0]], dtype=str)
                        #     bin_CG[0, index_minus:] = list(methy_reads[:cpg_pos_wanted.shape[0]-index_minus])
                        bin_reads_me_end = np.zeros([1,cpg_pos_wanted.shape[0]],int)
                        # print(f"Shape of bin_reads_me_end: {bin_reads_me_end}")
                        bin_reads_me_end[bin_CG == 'C'] = 1
                        # print(f"Shape of bin_reads_me_end: {bin_reads_me_end}")
                        # print(bin_reads_me_end,flush=True)
                        # sys.stdout.flush()
                        bin_reads_unme_end = np.zeros([1,cpg_pos_wanted.shape[0]],int)
                        bin_reads_unme_end[bin_CG == 'T'] = 1
                        # print()
                        # print(bin_reads_unme_end,flush=True)
                        # sys.stdout.flush()
                        bin_reads_me_read += bin_reads_me_end
                        bin_reads_unme_read += bin_reads_unme_end

                if bin_reads_me_read.sum() + bin_reads_unme_read.sum()>=2:
                    # print(bin_reads_me_read,flush=True)
                    # sys.stdout.flush()
                    # print(bin_reads_unme_read,flush=True)
                    # sys.stdout.flush()
                    for k in range(reads_chr.iloc[idx,3]):
                        bin_reads_me = np.vstack([bin_reads_me, bin_reads_me_read])
                        bin_reads_unme = np.vstack([bin_reads_unme, bin_reads_unme_read])
                # print(idx,bin_reads_me.shape,reads_chr.iloc[idx,:],flush=True)
                # sys.stdout.flush()
                # print(bin_reads_me_read,flush=True)
                # sys.stdout.flush()
                # print(bin_reads_unme,flush=True)
                # sys.stdout.flush()
                
                idx += 1
            # print(idx,bin_reads_me.shape,reads_chr.iloc[idx,:],flush=True)
            # sys.stdout.flush()
            # print(bin_reads_me_read,flush=True)
            # sys.stdout.flush()
            # print(bin_reads_unme,flush=True)
            # sys.stdout.flush()

            if bin_reads_me.shape[0] >= 1:
                bin_reads_total = bin_reads_me+bin_reads_unme
                methy_level = round(bin_reads_me.sum()*1.0/bin_reads_total.sum(),3)
                total_me = bin_reads_me.sum()
                total_cpg = bin_reads_total.sum()
            else:
                methy_level = np.nan
                total_me = 0
                total_cpg = 0
            if bin_reads_me.shape[0]>=min_read and bin_reads_me.shape[1]>=min_cpg and bin_reads_total.sum(axis = 0).max() >= min_overlap:
                if bin_reads_me.shape[0]>max_read:
                    f1 = random.sample(range(bin_reads_me.shape[0]),max_read)
                    bin_reads_me = bin_reads_me[f1,:]
                    bin_reads_unme = bin_reads_unme[f1,:]
                    bin_reads_total = bin_reads_total[f1,:]
                    downsample = 1
                else:
                    downsample = 0
                #concordant of reads
                reads_me = bin_reads_me.dot(bin_reads_me.T)
                reads_unme = bin_reads_unme.dot(bin_reads_unme.T)
                reads_total = bin_reads_total.dot(bin_reads_total.T)
                help_mat = np.ones([bin_reads_me.shape[0],bin_reads_me.shape[0]],int) - np.eye(bin_reads_me.shape[0],dtype = int)
                concordant_reads = round((reads_me*help_mat+reads_unme*help_mat).sum()*1.0/((reads_total*help_mat).sum()),3)
                #pvals for concordant of reads
                total_pair = (reads_total*help_mat).sum()
                wanted_pair = (reads_me*help_mat+reads_unme*help_mat).sum()
                pair_me_count = ((bin_reads_me.dot(bin_reads_total.T))*help_mat).sum()
                pair_unme_count = ((bin_reads_unme.dot(bin_reads_total.T))*help_mat).sum()
                pair_me_frac = pair_me_count*1.0/(pair_me_count+pair_unme_count)
                bino_rat = pair_me_frac*pair_me_frac+(1-pair_me_frac)*(1-pair_me_frac)
                exp_reads = round(bino_rat,3)
                if pair_me_count != 0 and pair_unme_count != 0:
                    if wanted_pair >= total_pair*bino_rat:
                        p_reads = 1-stats.binom.cdf(wanted_pair-1,total_pair,bino_rat)
                    else:
                        p_reads = stats.binom.cdf(wanted_pair,total_pair,bino_rat)
                else:
                    p_reads = 1
                #concordant of cpg site
                site_me = bin_reads_me.T.dot(bin_reads_me)
                site_unme = bin_reads_unme.T.dot(bin_reads_unme)
                site_total = bin_reads_total.T.dot(bin_reads_total)
                help_mat = np.ones([bin_reads_me.shape[1],bin_reads_me.shape[1]],int) - np.eye(bin_reads_me.shape[1],dtype = int)
                concordant_sites = round((site_me*help_mat+site_unme*help_mat).sum()*1.0/((site_total*help_mat).sum()),3)
                #pvals for concordant of CpGs
                total_pair = (site_total*help_mat).sum()
                wanted_pair = (site_me*help_mat+site_unme*help_mat).sum()
                pair_me_count = ((bin_reads_me.T.dot(bin_reads_total))*help_mat).sum()
                pair_unme_count = ((bin_reads_unme.T.dot(bin_reads_total))*help_mat).sum()
                pair_me_frac = pair_me_count*1.0/(pair_me_count+pair_unme_count)
                bino_rat = pair_me_frac*pair_me_frac+(1-pair_me_frac)*(1-pair_me_frac)
                exp_cpgs = round(bino_rat,3)
                if pair_me_count != 0 and pair_unme_count != 0:
                    if wanted_pair >= total_pair*bino_rat:
                        p_cpgs = 1-stats.binom.cdf(wanted_pair-1,total_pair,bino_rat)
                    else:
                        p_cpgs = stats.binom.cdf(wanted_pair,total_pair,bino_rat)
                else:
                    p_cpgs = 1
                
            else:
                downsample = 0
                concordant_reads = np.nan
                concordant_sites = np.nan
                p_reads = 1
                exp_reads = np.nan
                p_cpgs = 1
                exp_cpgs = np.nan
            if basedon0 == 1:
                pos1 = pos1-1 #to 0-based
            outdata.write(interval_name+ '\t'+\
                            chr+'\t'+\
                            str(pos1)+'\t'+\
                            str(pos2)+'\t'+\
                            str(bin_reads_me.shape[0])+'\t'+\
                            str(bin_reads_me.shape[1])+'\t'+\
                            str(total_me)+'\t'+\
                            str(total_cpg)+'\t'+\
                            str(methy_level)+'\t'+\
                            str(concordant_reads)+'\t'+\
                            str(concordant_sites)+'\t'+\
                            str(round(concordant_reads-exp_reads,3))+'\t'+\
                            '%.3e' % p_reads+'\t'+\
                            str(round(concordant_sites-exp_cpgs,3))+'\t'+\
                            '%.3e' % p_cpgs+'\t'+\
                            str(downsample)+'\n')
        else:
            pass
    outdata.close()
def helper_bam_to_meconcord(args):
    return bam_to_meconcord(*args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-i', dest='infile_path',  nargs='*', required=True, help='path of input bam')
    parser.add_argument('-p', dest='process', type=int, required=True, help='number of process')
    parser.add_argument('--chrlist', dest='chrlist', nargs='*', required=True, help='the list of chr')
    parser.add_argument('-r', dest='region_file', type=str, required=True, help='path of region')
    parser.add_argument('-b', dest='cpg_file_path', type=str, required=True, help='path of cpg file')
    args = parser.parse_args()
    infile_path = args.infile_path
    # outfile_path = '/home/lijia/wangyanni/project/integration/ATAC_Meth/nature_data_atac_meth/human_snmcseq3'
    process = args.process
    chr_list = args.chrlist
    region_file_raw = args.region_file
    cpg_file_path_raw = args.cpg_file_path
    print(infile_path)
    para_lists = []
    for i in infile_path:
        for chr in chr_list:
            cpg_file = f"CpG_{chr}.bed"
            cpg_file_path = cpg_file_path_raw
            cpg_file_path = os.path.join(cpg_file_path, cpg_file)
            print(infile_path)
            directory = os.path.dirname(i)
            bam_name = os.path.basename(i)
            region_file_name =f"region_{chr}.bed"
            region_file = region_file_raw
            region_file = os.path.join(region_file, region_file_name)
            para_list =[bam_name, chr, directory, cpg_file_path, region_file]
            para_lists.append(para_list)
    print(para_list)
    pool = multiprocessing.Pool(processes= process)
    pool.map(helper_bam_to_meconcord,para_lists)
    pool.close()
    








#     # create dir
#     if not os.path.exists('BAM_FILE'):
#         os.mkdir('BAM_FILE')
#     ##sra to fastq fastp cutadapt
#     sra_path = infile_path
#     sra_file = [file.split(".")[0] for file in os.listdir(sra_path) if file.endswith('.sra')]
#     fastq_file = list(set([file.split("_")[0] for file in os.listdir(sra_path) if file.endswith('fastq.gz')]))

#     if len(sra_file) > 0:
#         sra_need_pp = [file for file in sra_file if file not in fastq_file]
#         print(sra_need_pp)
#         if len(sra_need_pp) >0:
#             args_list = list(product([sra_path], sra_need_pp, [sra_path]))
#             print(args_list)
#             pool = multiprocessing.Pool(processes= process * 4)
#             pool.map(helper_sra_to_fastq,args_list)
#             pool.close()

#     align_path = 'BAM_FILE'
#     fastq_file_updata = list(set([file.split("_")[0] for file in os.listdir(sra_path) if file.endswith('fastq.gz')]))
#     bam_file = list(set([file.split("_")[0] for file in os.listdir(align_path) if file.endswith('mapq30.bam')]))
#     fastq_to_bam_file = [file for file in fastq_file_updata if file not in bam_file]
#     if len(fastq_to_bam_file) >0:
#         args_list = list(product([align_path], fastq_to_bam_file, [sra_path]))
#         pool = multiprocessing.Pool(processes=process)
#         pool.map(helper_fastp_to_bam,args_list)
#         pool.close()





    






# command0 = "wgbstools bam2pat {bamfile_path} {out_file_path}"

# command_1 = (f'zcat {input_path} | awk -F"\t" \'BEGIN {{OFS="\t"}} {{new_col = $2 + length($3) - 1; print $0 "\\t" new_col}}\' | awk -F"\t" \'NR==FNR {{a[$3]=$2; next}} {{print $0 "\\t" a[$2] "\\t" a[$5]+1}}\' CpG_chr19.bed - > {output_path}', shell=True, check=True)

# command1 = """
# zcat chr19_merge_9.pat.gz | awk -F'\\t' 'BEGIN {OFS="\\t"} {new_col = $2 + length($3) - 1; print $0 "\\t" new_col}' | \
# awk -F'\\t' 'NR==FNR {a[$3]=$2; next} {print $0 "\\t" a[$2] "\\t" a[$5]+1}' CpG_chr19.bed - > chr19_merge_9_new.pat
# """
# subprocess.run(command1, shell=True, check=True)

# pat_data = pd.read_csv('chr19_merge_9_new.pat',sep = '\t',header = None,)
# CG_bed = pd.read_csv('CpG_chr19.bed',sep = '\t',header = None,)
# region = pd.read_csv('test_region.bed',sep = '\t',header = None,)
# outdata = open('chr19_merge_9_pat_MeConcord.txt','w')
# min_cpg = 4
# min_read = 4
# min_overlap = 4
# max_read = 500
# filter_reads_cpgnum = 2
# basedon0 = 0
# max_dis = 600
# reads_chr = pat_data.sort_values(by=[5,6],ascending = [True,True])
# bin_chrom = region.loc[region.iloc[:,0]==chrom,:].copy()
# bin_sort = bin_chrom.sort_values(by=[2],ascending = [True])
# interval_size =5000
# i = 0
# idx = 0
# split_num = int(math.ceil((bin_sort.iloc[i,2]-bin_sort.iloc[i,1]+1)*1.0/interval_size))
# if basedon0 == 1:
#     interval_name = chrom+'_'+str(bin_sort.iloc[i,1]-1)+'_'+str(bin_sort.iloc[i,2])#to 0-based
# else:
#     interval_name = chrom+'_'+str(bin_sort.iloc[i,1])+'_'+str(bin_sort.iloc[i,2])
# for j in range(0,split_num):
#     pos1 = bin_sort.iloc[i,1]+j*interval_size ##区间的位置
#     # print(pos1)
#     pos2 = min(bin_sort.iloc[i,2], bin_sort.iloc[i,1]+(j+1)*interval_size-1)
#     # print(pos2)

#     cpg_pos_wanted = CG_bed.loc[(CG_bed.iloc[:,1]>=pos1) & (CG_bed.iloc[:,1]<=pos2),:]
#     # print(cpg_pos_wanted)
#     cpg_pos_help = CG_bed.loc[(CG_bed.iloc[:,1]>=(pos1-2000)) & (CG_bed.iloc[:,1]<=pos2),:]
#     # print(cpg_pos_help)
#     if cpg_pos_wanted.shape[0]>=min_cpg:##bin has at least cpgs
#         bin_reads_me = np.zeros([0,cpg_pos_wanted.shape[0]],int)
#         # print(bin_reads_me)
#         bin_reads_unme = np.zeros([0,cpg_pos_wanted.shape[0]],int)
#         # print(bin_reads_unme)
#         # print(idx,pos1,reads_sort.iloc[idx,3],reads_sort.iloc[idx,2])
#         while idx > 0 and pos1 - reads_chr.iloc[idx,6]< max_dis:
#             idx -= 1
#         while idx < (reads_chr.shape[0]-1) and pos1 - reads_chr.iloc[idx,6] > max_dis:
#             idx += 1
#         while idx < (reads_chr.shape[0]-1) and reads_chr.iloc[idx,5] < pos2: ##reads的起始小于区间终止
#             bin_reads_me_read = np.zeros([1,cpg_pos_wanted.shape[0]],int)
#             bin_reads_unme_read = np.zeros([1,cpg_pos_wanted.shape[0]],int)
#             ##原来reads的2对应现在reads的5
#             if len(reads_chr.iloc[idx,2]) >=1:                        
#                 if sum((cpg_pos_wanted.iloc[:,1]>=reads_chr.iloc[idx,5]) & (cpg_pos_wanted.iloc[:,1]<reads_chr.iloc[idx,6])) >= 1:
#                     # cell_id = reads_chr.iloc[idx,12]
#                     # cell_ids.append(cell_id)
#                     # print(len(cell_ids))
#                     lost_num = sum((reads_chr.iloc[idx,5]<=cpg_pos_help.iloc[:,1]) & (cpg_pos_wanted.iloc[0,1]>cpg_pos_help.iloc[:,1]))
#                     # print('lost_numA:',lost_num)
#                     methy_reads = reads_chr.iloc[idx,2][lost_num:]
#                     reads_CG = len(methy_reads)
#                     reads_start_index = reads_chr.iloc[idx,1]
#                     bin_start_index = cpg_pos_help.iloc[0,2]
#                     if lost_num > 0 and reads_CG <= cpg_pos_wanted.shape[0]:
#                         bin_CG = np.zeros([1, cpg_pos_wanted.shape[0]], dtype=str)
#                         bin_CG[0, :reads_CG] = list(methy_reads)
#                     if lost_num > 0 and reads_CG > cpg_pos_wanted.shape[0]:
#                         bin_CG = np.zeros([1, cpg_pos_wanted.shape[0]], dtype=str)
#                         bin_CG[0, :cpg_pos_wanted.shape[0]] = list(methy_reads[:cpg_pos_wanted.shape[0]])
#                     if lost_num == 0 and reads_start_index - bin_start_index + reads_CG <= cpg_pos_wanted.shape[0]:
#                         index_minus = reads_start_index - bin_start_index
#                         bin_CG = np.zeros([1, cpg_pos_wanted.shape[0]], dtype=str)
#                         print(bin_CG,methy_reads,index_minus)
#                         print(bin_CG[0, index_minus:],list(methy_reads))
#                         available_space = cpg_pos_wanted.shape[0] - index_minus
        
#                         # Ensure methy_reads fits into the available space in bin_CG
#                         if reads_CG <= available_space:
#                             bin_CG[0, index_minus:index_minus + reads_CG] = list(methy_reads)
#                         else:
#                             bin_CG[0, index_minus:index_minus + available_space] = list(methy_reads[:available_space])
#                         # bin_CG[0, index_minus:] = list(methy_reads)
#                     # if lost_num == 0 and reads_start_index - bin_start_index + reads_CG > cpg_pos_wanted.shape[0]:
#                     #     index_minus = reads_start_index - bin_start_index
#                     #     bin_CG = np.zeros([1, cpg_pos_wanted.shape[0]], dtype=str)
#                     #     bin_CG[0, index_minus:] = list(methy_reads[:cpg_pos_wanted.shape[0]-index_minus])
#                     bin_reads_me_end = np.zeros([1,cpg_pos_wanted.shape[0]],int)
#                     bin_reads_me_end[bin_CG == 'C'] = 1
#                     bin_reads_unme_end = np.zeros([1,cpg_pos_wanted.shape[0]],int)
#                     bin_reads_unme_end[bin_CG == 'T'] = 1
#                     bin_reads_me_read += bin_reads_me_end
#                     bin_reads_unme_read += bin_reads_unme_end

#             if bin_reads_me_read.sum() + bin_reads_unme_read.sum()>=2:
#                 bin_reads_me = np.vstack([bin_reads_me, bin_reads_me_read])
#                 bin_reads_unme = np.vstack([bin_reads_unme, bin_reads_unme_read])
            
#             idx += 1

#         if bin_reads_me.shape[0] >= 1:
#             bin_reads_total = bin_reads_me+bin_reads_unme
#             methy_level = round(bin_reads_me.sum()*1.0/bin_reads_total.sum(),3)
#             total_me = bin_reads_me.sum()
#             total_cpg = bin_reads_total.sum()
#         else:
#             methy_level = np.nan
#             total_me = 0
#             total_cpg = 0
#         if bin_reads_me.shape[0]>=min_read and bin_reads_me.shape[1]>=min_cpg and bin_reads_total.sum(axis = 0).max() >= min_overlap:
#             if bin_reads_me.shape[0]>max_read:
#                 f1 = random.sample(range(bin_reads_me.shape[0]),max_read)
#                 bin_reads_me = bin_reads_me[f1,:]
#                 bin_reads_unme = bin_reads_unme[f1,:]
#                 bin_reads_total = bin_reads_total[f1,:]
#                 downsample = 1
#             else:
#                 downsample = 0
#             #concordant of reads
#             reads_me = bin_reads_me.dot(bin_reads_me.T)
#             reads_unme = bin_reads_unme.dot(bin_reads_unme.T)
#             reads_total = bin_reads_total.dot(bin_reads_total.T)
#             help_mat = np.ones([bin_reads_me.shape[0],bin_reads_me.shape[0]],int) - np.eye(bin_reads_me.shape[0],dtype = int)
#             concordant_reads = round((reads_me*help_mat+reads_unme*help_mat).sum()*1.0/((reads_total*help_mat).sum()),3)
#             #pvals for concordant of reads
#             total_pair = (reads_total*help_mat).sum()
#             wanted_pair = (reads_me*help_mat+reads_unme*help_mat).sum()
#             pair_me_count = ((bin_reads_me.dot(bin_reads_total.T))*help_mat).sum()
#             pair_unme_count = ((bin_reads_unme.dot(bin_reads_total.T))*help_mat).sum()
#             pair_me_frac = pair_me_count*1.0/(pair_me_count+pair_unme_count)
#             bino_rat = pair_me_frac*pair_me_frac+(1-pair_me_frac)*(1-pair_me_frac)
#             exp_reads = round(bino_rat,3)
#             if pair_me_count != 0 and pair_unme_count != 0:
#                 if wanted_pair >= total_pair*bino_rat:
#                     p_reads = 1-stats.binom.cdf(wanted_pair-1,total_pair,bino_rat)
#                 else:
#                     p_reads = stats.binom.cdf(wanted_pair,total_pair,bino_rat)
#             else:
#                 p_reads = 1
#             #concordant of cpg site
#             site_me = bin_reads_me.T.dot(bin_reads_me)
#             site_unme = bin_reads_unme.T.dot(bin_reads_unme)
#             site_total = bin_reads_total.T.dot(bin_reads_total)
#             help_mat = np.ones([bin_reads_me.shape[1],bin_reads_me.shape[1]],int) - np.eye(bin_reads_me.shape[1],dtype = int)
#             concordant_sites = round((site_me*help_mat+site_unme*help_mat).sum()*1.0/((site_total*help_mat).sum()),3)
#             #pvals for concordant of CpGs
#             total_pair = (site_total*help_mat).sum()
#             wanted_pair = (site_me*help_mat+site_unme*help_mat).sum()
#             pair_me_count = ((bin_reads_me.T.dot(bin_reads_total))*help_mat).sum()
#             pair_unme_count = ((bin_reads_unme.T.dot(bin_reads_total))*help_mat).sum()
#             pair_me_frac = pair_me_count*1.0/(pair_me_count+pair_unme_count)
#             bino_rat = pair_me_frac*pair_me_frac+(1-pair_me_frac)*(1-pair_me_frac)
#             exp_cpgs = round(bino_rat,3)
#             if pair_me_count != 0 and pair_unme_count != 0:
#                 if wanted_pair >= total_pair*bino_rat:
#                     p_cpgs = 1-stats.binom.cdf(wanted_pair-1,total_pair,bino_rat)
#                 else:
#                     p_cpgs = stats.binom.cdf(wanted_pair,total_pair,bino_rat)
#             else:
#                 p_cpgs = 1
            
#         else:
#             downsample = 0
#             concordant_reads = np.nan
#             concordant_sites = np.nan
#             p_reads = 1
#             exp_reads = np.nan
#             p_cpgs = 1
#             exp_cpgs = np.nan
#         if basedon0 == 1:
#             pos1 = pos1-1 #to 0-based
#         outdata.write(interval_name+ '\t'+\
#                         chrom+'\t'+\
#                         str(pos1)+'\t'+\
#                         str(pos2)+'\t'+\
#                         str(bin_reads_me.shape[0])+'\t'+\
#                         str(bin_reads_me.shape[1])+'\t'+\
#                         str(total_me)+'\t'+\
#                         str(total_cpg)+'\t'+\
#                         str(methy_level)+'\t'+\
#                         str(concordant_reads)+'\t'+\
#                         str(concordant_sites)+'\t'+\
#                         str(round(concordant_reads-exp_reads,3))+'\t'+\
#                         '%.3e' % p_reads+'\t'+\
#                         str(round(concordant_sites-exp_cpgs,3))+'\t'+\
#                         '%.3e' % p_cpgs+'\t'+\
#                         str(downsample)+'\n')
#     else:
#         pass
# outdata.close()
            


            
