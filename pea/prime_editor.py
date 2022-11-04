import logging
import sys
import pkg_resources
import yaml
from logging.config import fileConfig
from pathlib import Path

import pandas as pd
from absl import app, flags
from termcolor import colored, cprint

from .util import NamedArgs, TargetRegion, UserRegion, get_aligner
from .util import matchUpto1, get_user_region_query, revertedSeq


fileConfig(pkg_resources.resource_filename('pea', 'log.ini'))
logger = logging.getLogger('align_mutations')


FLAGS = flags.FLAGS

flags.DEFINE_string('amplicon_seq', None, "Amplicon sequence")
flags.DEFINE_string('target_seq', None, "Target RGEN sequence")
flags.DEFINE_string('user_target_mutation', None, "User sequence with desired mutations")

flags.DEFINE_string('output_nametag', 'out', "Output filename tag")
flags.DEFINE_string('input', None, "Input fastqjoin file")

flags.DEFINE_integer('comparison_radius', 60, "Radius for a comparison range")
flags.DEFINE_integer('indicator_seq_length', 15, "Length of indicator sequences")
flags.DEFINE_integer('user_region_length', 30, "Length for a comparison range")
flags.DEFINE_integer('user_region_beg_offset', 3, "Starting offset of user region, from a PAM start position.")
flags.DEFINE_integer('pam_length', 3, "Length of PAM seq; e.g. 3 for NGG")

flags.mark_flag_as_required('amplicon_seq')
flags.mark_flag_as_required('target_seq')
flags.mark_flag_as_required('input')


def main(argv):
    del argv

    aseq_arg = NamedArgs.build(FLAGS.amplicon_seq).upper()
    aseq = aseq_arg.str
    target_seq_arg = NamedArgs.build(FLAGS.target_seq).upper()
    target_seq = target_seq_arg.str
    user_target_mutation_arg = NamedArgs.build(FLAGS.user_target_mutation).upper()
    user_target_mutation = user_target_mutation_arg.str
    
    
    output_nametag = FLAGS.output_nametag
    input_fqj = FLAGS.input    
    fastqjoin_path = Path(input_fqj)
    
    comparison_radius = FLAGS.comparison_radius
    user_region_length = FLAGS.user_region_length
    indicator_seq_length = FLAGS.indicator_seq_length
    user_region_beg_offset = FLAGS.user_region_beg_offset
    pam_len = FLAGS.pam_length
    
    tr = TargetRegion.build(aseq, target_seq, comparison_radius, pam_len)
    if tr is None:
        logger.info("Cannot find target sequence in the amplicon sequence.")
        logger.info("Find target sequence in the reverted amplicon sequence.")
        aseq = revertedSeq(aseq)
        tr = TargetRegion.build(aseq, target_seq, comparison_radius, pam_len)
    if tr is None:
#         raise RuntimeError("Cannot find the target sequence in the amplicon_seq.")
        logger.error("cannot find the target sequence in the amplicon_seq")
        logger.error(f"Target sequence: {target_seq}")
        logger.error(f"Amplicon sequence: {aseq}")
        sys.exit(1)
    ur = UserRegion.build(tr.pam_beg, user_region_beg_offset, user_region_length)
    idc_l = tr.left_indicator_seq(indicator_seq_length)
    idc_r = tr.right_indicator_seq(indicator_seq_length)

    user_region_seq = tr.get_slice(ur.beg, ur.end)
    
    
    logger.info(f"입력데이터: {fastqjoin_path}")
    
    colored_target_seq=colored(aseq[tr.beg:tr.pam_beg], 'green'
                                        )+colored(aseq[tr.pam_beg:tr.end], 'green', attrs=['bold', 'reverse'])
    colored_aseq = aseq[:tr.beg]+colored_target_seq+aseq[tr.end:]
    colored_aseq2 = aseq[:ur.beg]+colored(aseq[ur.beg:ur.end], 'yellow', attrs=['underline'])+aseq[ur.end:]
    
    
    logger.info(f"Green color: {colored('target sequence', 'green')} with {colored('PAM', 'green', attrs=['bold', 'reverse'])}.")    
    logger.info(f"Target seq: {colored_target_seq}")
    logger.info(f"Amplicon seq: {colored_aseq}")
    logger.info(f"Amplicon seq: {colored_aseq2}")
    logger.info(f"Yellow color: {colored('User defined region', 'yellow', attrs=['underline'])} to align with fastqjoin reads.")
    


    # print(tr)
    # print(ur)
    # print(idc_l)
    # print(idc_r)

    
    seq_lines_fwd = open(fastqjoin_path).read().split('\n')[1::4]
    seq_lines_rev = [revertedSeq(x) for x in seq_lines_fwd]
    seq_lines = seq_lines_fwd + seq_lines_rev
    
    left = [matchUpto1(idc_l, x) for x in seq_lines]
    right = [matchUpto1(idc_r, x) + len(idc_r) for x in seq_lines]
    
    query_seqs = [x[b:e] for x, b, e in zip(seq_lines, left, right) if 0 < b < e]
    df_query_seqs = pd.DataFrame(query_seqs, columns=['query_seq'])

    counts = df_query_seqs.query_seq.value_counts().to_frame('n')  # .reset_index()
    counts.index.name = 'query_seq'
    seq_counts = counts.reset_index()

    aligner = get_aligner()

    seq_counts['user_region_query'] = seq_counts.query_seq.apply(lambda x: get_user_region_query(aligner, x, tr, ur))
    user_region_query_counts = seq_counts.groupby('user_region_query').n.sum().to_frame('n_reads').sort_values('n_reads', ascending=False)

    df_mut_all = user_region_query_counts.query('index!=@user_region_seq').copy()
    df_mut_all['seq_len']=df_mut_all.index.str.len()
    n_all_mutation = df_mut_all.n_reads.sum()
    n_target_mutation=df_mut_all.query('index==@user_target_mutation').n_reads.sum()
    n_incorrect_indels=df_mut_all.query(f'index!=@user_target_mutation and seq_len!={len(user_region_seq)}').n_reads.sum()
    n_mut_others = n_all_mutation - n_target_mutation - n_incorrect_indels
                                                                                             
    n_aligned = user_region_query_counts.n_reads.sum()
    n_all = len(seq_lines) // 2 #divide by 2 since seq_lines = seq_lines_fwd + seq_lines_rev

    fqname = fastqjoin_path.name
    if fastqjoin_path.parent != Path('.'):
        fqname = f'{fastqjoin_path.parent.name}.{fqname}'
    prefix = f'{fqname}.{aseq_arg}.{target_seq_arg}.{user_target_mutation_arg}.prime_editor.{output_nametag}'
    with open(f'{prefix}.summary.txt', 'w') as f:        
        f.write(f'{fqname}\t{output_nametag}\t{aseq_arg}\t{target_seq}\t{user_target_mutation}\t'
                f'{n_target_mutation}\t{n_target_mutation/n_aligned:.5f}\t'
                f'{n_incorrect_indels}\t{n_incorrect_indels/n_aligned:.5f}\t'
                f'{n_mut_others}\t{n_mut_others/n_aligned:.5f}\t'
                # f'{n_all_mutation}\t{n_all_mutation/n_aligned:.5f}\t'
                f'{n_aligned}\t{n_all}\n')
    user_region_query_counts.to_csv(f'{prefix}.count.csv')


if __name__ == '__main__':
    app.run(main)
