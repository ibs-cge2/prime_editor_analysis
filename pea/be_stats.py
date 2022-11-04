import glob
import numpy as np
import pandas as pd
from absl import app, flags


FLAGS = flags.FLAGS

flags.DEFINE_multi_string('be_mut', [], "Base editing mutation of interest. Can be specified multiple times.")
flags.DEFINE_multi_string('nonX_mut', [], "X to non-X mutation of interest. Can be specified multiple times.")
flags.DEFINE_string('input_file_glob', '*.fastqjoin.*.align.csv', "Glob pattern for input files")
flags.DEFINE_string('be_output', 'summary.base_editing.csv', "Output filename")
flags.DEFINE_string('be_output_overall', 'summary.be_overall.csv', "Output filename")
flags.DEFINE_string('nonX_output', 'summary.mutations.csv', "Output filename")
flags.DEFINE_string('nonX_detail_output', 'summary.nonX_per_pos.mutations.csv', "Output filename")
flags.DEFINE_bool('exclude_indel_reads', True, "Do not use reads with indels for base editing stats")

# flags.mark_flag_as_required('be_mut')

be_mut_values = [f'{x}to{y}' for x in list("ACGT") for y in list("ACGT") if x!=y]
flags.register_validator('be_mut',
                         lambda xs: all([x in be_mut_values for x in xs]),
                         message=f'--be_mut must be one of {be_mut_values}')

nonX_mut_values = list("ACGT")
flags.register_validator('nonX_mut',
                         lambda xs: all([x in nonX_mut_values for x in xs]),
                         message=f'--nonX_output must be one of {nonX_mut_values}')



def base_editing_mask(df):
    df_ref=pd.DataFrame([[c for c in x] for x in df.ref], index=df.index)
    df_alignment=pd.DataFrame([[c for c in x] for x in df.alignment], index=df.index)
    df_read=pd.DataFrame([[c for c in x] for x in df.read], index=df.index)
    
    ref_cols = [f'{c}{i+1}' for i, c in enumerate(df.ref.iloc[0].replace('-',''))]
    df_ref.columns = ref_cols
    df_alignment.columns = ref_cols
    df_read.columns = ref_cols
    
    df_ref=df_ref[(df_alignment=='X') | (df_alignment=='.')]
    df_read=df_read[(df_alignment=='X') | (df_alignment=='.')]
    mutations=df_ref+"to"+df_read
    return mutations


def base_editing_per_pos(df, mut):
    tmp=base_editing_mask(df)

    mut_match=(tmp==mut)
    mut_match['any']=mut_match.any(axis=1)
    mut_match *=df.n_reads.values.reshape((-1,1))
    mut_match['filename']=df.filename
    mut_match['mut']=mut
    return mut_match


def base_editing_summary(df, target_muts):
    mutations = base_editing_mask(df)
    
    tmp=pd.concat([(mutations == t).any(axis=1).to_frame(t) for t in target_muts], axis=1)
    tmp['any']=tmp.any(axis=1)
    tmp*=df.n_reads.values.reshape((-1,1))
    tmp['filename']=df.filename
    return tmp.groupby('filename').sum()


def x_to_nonX_mutation_count(df, bps):
    df_ref=pd.DataFrame([[c for c in x] for x in df.ref], index=df.index)
    df_alignment=pd.DataFrame([[c for c in x] for x in df.alignment], index=df.index)
    
    mutations=df_ref[(df_alignment=='X') | (df_alignment=='.')]
    count = pd.concat([(mutations==bp).any(axis=1).to_frame(bp) for bp in bps], axis=1)
    count['any']=count.any(axis=1)
    count *= df[['n_reads']].values
    return count


def x_to_nonX_mutation_per_pos_count(df):
    df_ref=pd.DataFrame([[c for c in x] for x in df.ref], index=df.index)
    df_alignment=pd.DataFrame([[c for c in x] for x in df.alignment], index=df.index)
    
    ref_cols = [f'{c}{i+1}' for i, c in enumerate(df.ref.iloc[0].replace('-',''))]
    
    mutations=df_ref[(df_alignment=='X') | (df_alignment=='.')].notnull()
    mutations['mut_any']=mutations.any(axis=1)
    all_muts=mutations*df.n_reads.values.reshape(-1,1)    
    all_muts_per_pos = all_muts.sum(axis=0)
    all_muts_per_pos.index = ref_cols+['n_mut_any']
    return all_muts_per_pos


def mask_insertion(row):
    mask=np.array(list(row.ref))!='-'
    f=lambda x:''.join(np.array(list(x))[mask])
    return [f(row.ref), f(row.alignment), f(row.read), row.n_reads]#, names=['ref','alignment','read','count']


def main(argv):
    del argv
    input_file_glob=FLAGS.input_file_glob
    be_muts=FLAGS.be_mut
    nonX_muts=FLAGS.nonX_mut
    be_output=FLAGS.be_output
    be_output_overall = FLAGS.be_output_overall
    nonX_output=FLAGS.nonX_output
    nonX_detail_output = FLAGS.nonX_detail_output
    exclude_indel_reads=FLAGS.exclude_indel_reads
    dfs_all={x:pd.read_csv(x) for x in sorted(glob.glob(input_file_glob))}
    for k, df in dfs_all.items():
        if df.empty:
            print(f"No reads to analyze: {k}")
    
    dfs = {k:df.dropna()  for k,df in dfs_all.items() if not df.empty} #Drop empty data
    
    if exclude_indel_reads:
        dfs = {k:df[~df.alignment.str.contains("-").astype(bool)] for k,df in dfs_all.items()} #No-indels. For Python 3.7
    else:
        f_mask_insertion=lambda df : pd.DataFrame(df.apply(mask_insertion, axis=1).to_list(),
                                                  columns=['ref','alignment','read','n_reads'])
        dfs = {k:f_mask_insertion(df) for k,df in dfs.items() if not df.empty}
        
    dfs = {k:df for k,df in dfs.items() if not df.empty}
    
    def append_col(df, name, col):
        df=df.copy()
        df[name]=col
        return df
    df_all_mutation_raw = pd.concat([append_col(df, "filename", k) for k,df in dfs.items()])
    df_all_mutation_raw.to_csv("all_mutation_raw.csv", index=False)
       
    df_read_count = pd.DataFrame({name: df['n_reads'].sum() for name,df in dfs.items()}, index=['n_aligned_noindel']).T    
    df_read_count_all = pd.DataFrame({name: df['n_reads'].sum() for name,df in dfs_all.items()}, index=['n_aligned_total']).T
    df_read_counts = pd.concat([df_read_count, df_read_count_all], axis=1, sort=False).dropna().astype(np.int32)
    df_read_counts.index.name='filename'    
    df_read_counts.to_csv('read_counts.csv')   
    
    
    if be_muts:
        df_be = pd.concat([base_editing_per_pos(df_all_mutation_raw, m) for m in be_muts])
        be_per_pos = df_be.groupby(['mut','filename']).sum().reset_index()
        be_per_pos = be_per_pos.merge(df_read_counts.reset_index(), on='filename')
        be_per_pos = be_per_pos.set_index(['mut','filename','n_aligned_noindel','n_aligned_total']).sort_index()
        be_per_pos.to_csv(be_output, sep='\t')
        
        df_be_summary = base_editing_summary(df_all_mutation_raw, be_muts)
        df_be_summary = df_be_summary.merge(df_read_counts.reset_index(), on='filename')
        df_be_summary.to_csv(be_output_overall, sep='\t', index=False)
    
    if nonX_muts:
        nonXmuts={x:x_to_nonX_mutation_count(df, nonX_muts).sum(axis=0) for x, df in dfs.items()}
        df_nonX=pd.DataFrame.from_dict(nonXmuts).T
        df_nonX.index.name = 'filename'
        df_nonX=df_nonX.merge(df_read_counts, left_index=True, right_index=True)
        df_nonX.to_csv(nonX_output, sep='\t')
        
        nonX_per_pos = pd.DataFrame.from_dict({x:x_to_nonX_mutation_per_pos_count(df) for x,df in dfs.items()}).T
        nonX_per_pos=nonX_per_pos.merge(df_read_counts, left_index=True, right_index=True)
        nonX_per_pos.index.name = 'nonX_per_pos'
        nonX_per_pos.to_csv(nonX_detail_output, sep='\t')
        

if __name__=='__main__':
    app.run(main)

