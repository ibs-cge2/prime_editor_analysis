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
from .util import matchUpto1, revertedSeq, align_read_to_user_region


fileConfig(pkg_resources.resource_filename('pea', 'log.ini'))
logger = logging.getLogger('align_mutations')


FLAGS = flags.FLAGS

flags.DEFINE_string('amplicon_seq', None, "Amplicon sequence")
flags.DEFINE_string('target_seq', None, "Target RGEN sequence")

flags.DEFINE_string('output_nametag', 'out', "Output filename tag")
flags.DEFINE_string('input', None, "Input fastqjoin file")

flags.DEFINE_integer('comparison_radius', 60, "Radius for a comparison range")
flags.DEFINE_integer('indicator_seq_length', 15, "Length of indicator sequences")
flags.DEFINE_integer('user_region_length', 30, "Length for a comparison range")
flags.DEFINE_integer('user_region_beg_offset', 3, "Starting offset of user region based against a PAM position.")
flags.DEFINE_integer('pam_length', 3, "Length of PAM seq; e.g. 3 for NGG")

flags.DEFINE_bool('indel_in_alignment', True, "Flag for allowing indels during sequence alignment")


flags.register_validator('comparison_radius',
                         lambda x: x > 20,
                         message='--comparison_radius should be larger than 20')

flags.mark_flag_as_required('amplicon_seq')
flags.mark_flag_as_required('target_seq')
flags.mark_flag_as_required('input')


def main(argv):
    del argv
    
    logger.info('Program start.')

    aseq_arg = NamedArgs.build(FLAGS.amplicon_seq).upper()
    aseq = aseq_arg.str
    target_seq_arg = NamedArgs.build(FLAGS.target_seq).upper()
    target_seq = target_seq_arg.str
    
    output_nametag = FLAGS.output_nametag
    input_fqj = FLAGS.input
    
    comparison_radius = FLAGS.comparison_radius
    user_region_length = FLAGS.user_region_length
    indicator_seq_length = FLAGS.indicator_seq_length
    user_region_beg_offset = FLAGS.user_region_beg_offset    
    pam_len = FLAGS.pam_length

    open_gap_score= -2.01 if FLAGS.indel_in_alignment else -20
    
    tr = TargetRegion.build(aseq, target_seq, comparison_radius, pam_len)    
    if tr is None:
        logger.info("Cannot find target sequence in the amplicon sequence.")
        logger.info("Find target sequence in the reverted amplicon sequence.")
        aseq = revertedSeq(aseq)
        tr = TargetRegion.build(aseq, target_seq, comparison_radius, pam_len)
    if tr is None:
        logger.error("cannot find the target sequence in the amplicon_seq")
        logger.error(f"Target sequence: {target_seq}")
        logger.error(f"Amplicon sequence: {aseq}")
        sys.exit(1)
    ur = UserRegion.build(tr.pam_beg, user_region_beg_offset, user_region_length)
    idc_l = tr.left_indicator_seq(indicator_seq_length)
    idc_r = tr.right_indicator_seq(indicator_seq_length)

    user_region_seq = tr.get_slice(ur.beg, ur.end)
    fastqjoin_path = Path(input_fqj)
    
    if not fastqjoin_path.exists():        
        logger.error(f"No such file or directory: {fastqjoin_path}")
        sys.exit(0)
    
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
        
    
    seq_lines_fwd = open(fastqjoin_path).read().split('\n')[1::4]
    seq_lines_rev = [revertedSeq(x) for x in seq_lines_fwd]
    seq_lines = seq_lines_fwd + seq_lines_rev
    
    left = [matchUpto1(idc_l, x) for x in seq_lines]
    right = [matchUpto1(idc_r, x, reverse=True) for x in seq_lines]
    
    query_seqs = [x[b:e+len(idc_r)] for x, b, e in zip(seq_lines, left, right) if 0 <= b < e]
    df_query_seqs = pd.DataFrame(query_seqs, columns=['query_seq'])
    
    logger.info(f"Total {len(seq_lines_fwd)} reads in the fastqjoin file.")
    logger.info(f"{len(query_seqs)} reads are aligned to the amplicon sequence.")

    counts = df_query_seqs.query_seq.value_counts().to_frame('n')  # .reset_index()
    counts.index.name = 'query_seq'
    seq_counts = counts.reset_index()

    aligner = get_aligner(open_gap_score)
    
    df=seq_counts.apply(lambda x: align_read_to_user_region(aligner, x, tr, ur), axis=1, result_type="expand")
    align_count = df.groupby(['ref','alignment','read']).n.sum().to_frame("n_reads").sort_values("n_reads", ascending=False).reset_index()
    
    align_count['ratio'] = align_count.n_reads / align_count.n_reads.sum()
    
    prefix = f'{fastqjoin_path.name}.{aseq_arg}.{target_seq_arg}.align_mutations.{output_nametag}'
    if fastqjoin_path.parent != Path('.'):
        prefix = f'{fastqjoin_path.parent.name}.'+prefix
        
    align_count.to_csv(f'{prefix}.align.csv', index=False)
    logger.info(f"Save output to {f'{prefix}.align.csv'}")
    logger.info("Program end.")
    


if __name__ == '__main__':
    app.run(main)
