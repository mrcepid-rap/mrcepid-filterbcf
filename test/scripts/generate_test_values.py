import sys
import sh
from pathlib import Path

import pandas as pd

def set_null(gt, ALT, dp, gq, binom_p, is_indel, dp_snp_thresh, dp_indel_thresh, binom_thresh):

    if (ALT != "*" and
            ((is_indel == False and dp >= dp_snp_thresh and (
                    (gt == "0/0" and gq >= 20) or
                    (gt == "0/1" and gq >= 20 and binom_p > binom_thresh) or
                    (gt == "1/1"))) or
             (is_indel == True and dp >= dp_indel_thresh and gq >= 20))):
        return gt
    else:
        return './.'

def calc_missingness(vec):

    tot_miss = 0
    for gt_n, gt in enumerate(vec):
        if gt == './.':
            tot_miss += 1
    return tot_miss / gt_n

def calc_ac(vec):

    tot_ac_zero = 0
    for gt_n, gt in enumerate(vec):
        if gt == '0/1':
            tot_ac_zero += 1
        elif gt == '1/1':
            tot_ac_zero += 2
        else:
            pass
    return tot_ac_zero

# sh.Rscript(plot_report, output_path, created_after_date, created_before_date)
vcf_location = Path(sys.argv[1])
binom = float(sys.argv[2]) # 0.001
dp_snp = int(sys.argv[3]) # 10
dp_indel = int(sys.argv[4]) # 15
miss = float(sys.argv[5]) # 0.1

sh.bcftools.query('-HH', '-f', "[%POS\\t%GT\\t%DP\\t%GQ\\t%REF\\t%ALT\\t%AD\\t%PBINOM(AD)\\n]", '-o', 'genotypes.txt', vcf_location)
gt_file = Path('genotypes.txt')

gt_frame = pd.read_csv(gt_file, sep='\t', na_values=".")
gt_frame['binom_p'] = gt_frame['AD.1'].apply(lambda x: 10**(-x / 10))
gt_frame['is_indel'] = gt_frame.apply(lambda x: len(x['REF']) != len(x['ALT']), axis=1)
gt_frame['GT'] = gt_frame['GT'].apply(lambda x: '0/1' if x == '1/0' else x)
gt_frame['ID'] = gt_frame.apply(lambda x: f'{x["#POS"]}_{x["REF"]}_{x["ALT"]}', axis=1)

gt_frame['adj_GT'] = gt_frame.apply(lambda x: set_null(x['GT'], x['ALT'], x['DP'], x['GQ'], x['binom_p'], x['is_indel'], dp_snp, dp_indel, binom), axis=1)
miss_vals = gt_frame.groupby('ID').aggregate(miss=pd.NamedAgg(column="adj_GT", aggfunc=calc_missingness),
                                             ac=pd.NamedAgg(column="adj_GT", aggfunc=calc_ac))
miss_n = len(miss_vals) - len(miss_vals.query('miss <= @miss & ac != 0'))

print(f'Number of ./. genotypes pre-filtering: {gt_frame["GT"].value_counts()["./."]}')
print(f'Number of ./. genotypes post-filtering: {gt_frame["adj_GT"].value_counts()["./."]}')
print(f'Number of missing sites with missingness <= {miss} and AC != 0: {miss_n}')