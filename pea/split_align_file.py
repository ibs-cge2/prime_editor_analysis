import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from absl import app, flags


FLAGS = flags.FLAGS

flags.DEFINE_multi_integer('pos', None, "1-based position to use as flags")
flags.DEFINE_integer('beg', None, "1-based start position of the subsequence of interest")
flags.DEFINE_integer('end', None, "1-based end position of the subsequence of interest")
flags.DEFINE_string('subdir', 'trimmed', "Sub-directory name to store trimmed align files")

flags.mark_flag_as_required('pos')
flags.mark_flag_as_required('beg')
flags.mark_flag_as_required('end')



def trim_seq(df, idx_beg, idx_end):
    trimmed = pd.concat([df.ref.str[idx_beg:idx_end],
                         df.alignment.str[idx_beg:idx_end],
                         df.read.str[idx_beg:idx_end],
                         df.drop(['ref','alignment','read'], axis=1)
                        ], axis=1)
    return trimmed


def get_tag_element(row, one_based_idx):
    aligner = list(np.cumsum(np.array(list(row.ref))!='-')-1)
    i_aligned =  aligner.index(one_based_idx-1)
    return row.ref[i_aligned] + str(one_based_idx) + row.read[i_aligned]


def main(argv):
    del argv
    
    poss = FLAGS.pos
    subdir = FLAGS.subdir
    
    idx_beg = FLAGS.beg
    idx_end = FLAGS.end
        
    dfs={x:pd.read_csv(x) for x in sorted(glob.glob('*.align.csv'))}
    
    Path(subdir).mkdir(parents=True, exist_ok=True)
    

    for name, df in dfs.items():
        df['tag']=df.apply(lambda x: '_'.join([get_tag_element(x, i) for i in poss]), axis=1)
        for tag,dg in df.groupby('tag'):
            new_name = f'{subdir}/{tag}.{name}'
            if '-' in tag:
                continue
            trim_seq(dg, idx_beg, idx_end).to_csv(new_name, index=False)

if __name__=='__main__':
    app.run(main)

