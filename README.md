## Prime editor analysis
Prime editor data analysis program

### Installation
```
#After downloading this repo,
pip install .
```

### Usage
#### Prime editor analysis
```
USAGE: prime_editor.py [flags]
flags:

  --amplicon_seq: Amplicon sequence
  --comparison_radius: Radius for a comparison range
    (default: '60')
    (an integer)
  --indicator_seq_length: Length of indicator sequences
    (default: '15')
    (an integer)
  --input: Input fastqjoin file
  --output_nametag: Output filename tag
    (default: 'out')
  --pam_length: Length of PAM seq; e.g. 3 for NGG
    (default: '3')
    (an integer)
  --target_seq: Target RGEN sequence
  --user_region_beg_offset: Starting offset of user region, from a PAM start position.
    (default: '3')
    (an integer)
  --user_region_length: Length for a comparison range
    (default: '30')
    (an integer)
  --user_target_mutation: User sequence with desired mutations



USAGE: align_mutations.py [flags]
flags:
  --amplicon_seq: Amplicon sequence
  --comparison_radius: Radius for a comparison range
    (default: '60')
    (an integer)
  --[no]indel_in_alignment: Flag for allowing indels during sequence alignment
    (default: 'true')
  --indicator_seq_length: Length of indicator sequences
    (default: '15')
    (an integer)
  --input: Input fastqjoin file
  --output_nametag: Output filename tag
    (default: 'out')
  --pam_length: Length of PAM seq; e.g. 3 for NGG
    (default: '3')
    (an integer)
  --target_seq: Target RGEN sequence
  --user_region_beg_offset: Starting offset of user region based against a PAM position.
    (default: '3')
    (an integer)
  --user_region_length: Length for a comparison range
    (default: '30')
    (an integer)



USAGE: be_stats.py [flags]
flags:
  --be_mut: Base editing mutation of interest. Can be specified multiple times.;
    repeat this option to specify a list of values
    (default: '[]')
  --be_output: Output filename
    (default: 'summary.base_editing.csv')
  --be_output_overall: Output filename
    (default: 'summary.be_overall.csv')
  --[no]exclude_indel_reads: Do not use reads with indels for base editing stats
    (default: 'true')
  --input_file_glob: Glob pattern for input files
    (default: '*.fastqjoin.*.align.csv')
  --nonX_detail_output: Output filename
    (default: 'summary.nonX_per_pos.mutations.csv')
  --nonX_mut: X to non-X mutation of interest. Can be specified multiple times.;
    repeat this option to specify a list of values
    (default: '[]')
  --nonX_output: Output filename
    (default: 'summary.mutations.csv')
```

#### MAUND program
```
USAGE: maund_default.py [-h] [-c COMPARISON_RANGE] [-b WINDOW_BEG]
                        [-e WINDOW_END] [-ib IDXSEQ_BEG] [-ie IDXSEQ_END]
                        [-t {A,C,G,T}] [-mcut MISMATCH_CUTOFF]
                        aseq rgen [files [files ...]]

positional arguments:
  aseq
  rgen
  files

optional arguments:
  -h, --help            show this help message and exit
  -c COMPARISON_RANGE, --comparison_range COMPARISON_RANGE
  -b WINDOW_BEG, --window_beg WINDOW_BEG
                        The 1-based start index of the window sequence.
                        (default: 4)
  -e WINDOW_END, --window_end WINDOW_END
                        The 1-based end index of the window sequence; last
                        inclusive. (default: 7)
  -ib IDXSEQ_BEG, --idxseq_beg IDXSEQ_BEG
                        The 1-based start index of the index sequence.
                        (default: 13)
  -ie IDXSEQ_END, --idxseq_end IDXSEQ_END
                        The 1-based end index of the index sequence; last
                        inclusive. (default: 22)
  -t {A,C,G,T}, --target_nt {A,C,G,T}
                        Nucleobase to watch mutations in the window sequence.
                        (default: A)
  -mcut MISMATCH_CUTOFF, --mismatch_cutoff MISMATCH_CUTOFF

```

### Run example

```
#Analysis 1
maund_default.py -t T -b 18 -e 18 CTCCCTAGGTGCTGGCTTCcagcccagccaaacttgtcaaccagtatcccggtgcaggagctgcacatactagcccctgtctaggaaaagctgtcctgcgacgccctctggaggaagcagggcttcctttcctctgccatcacgtgctcagtctgggccccaaggattgacccaggccagggctggagaagcagaaaaaaagcaTCAAGCCTACAAATGCATGC	GGCCCAGACTGAGCACGTGATGG	 1.fastqjoin


#Analysis 2
python -m gea.prime_editor  --user_region_length 15 --user_region_beg_offset 6 --amplicon_seq aseq1:cttggagagttttaagcaagggctgatgtgggctgcctagaaaggcatggatgagagaagcctggagacagggatcccagggaaacgcccatgcaattagtctatttctgctgcaagtaagcatgcatttgtaggcttgatgctttttttctgcttctccagccctggcctgggtcaatccttggggcccagactgagcacgtgatggcagaggaaaggaagccctgcttcctccagagggcgtcgcaggacagcttttcctagacaggggctagtatgtgc --target_seq trg1:ggcccagactgagcacgtga --user_target_mutation mut1:agcacgGATTACAAGGATGACGACGATAAGtgatggcag --input 73.fastqjoin --output_nametag opts1

#Analysis 3
python -m pea.align_mutations  --user_region_length 149 --user_region_beg_offset 79 --amplicon_seq aseq1:ccctggtcaacctcaacctaggcctcctatttattctagccacctctagcctagccgtttactcaatcctctgatcagggtgagcatcaaactcaaactacgccctgatcggcgcactgcgagcagtagcccaaacaatctcatatgaagtcaccctagccatcattctactatcaacattactaataagtggctcctttaacctctccaccctt --target_seq site1:actcaatcctctgatc --input 2.fastqjoin --output_nametag opts1

python -m pea.be_stats --nonX_mut C --nonX_mut G
```

### Output files
- `{input}.{target_seq}.out.Miseq_summary.txt` : result summary
- `{input}.{target_seq}.out._window.txt` : window-read counts
- `{input}.{target_seq}.out._aligned.txt`: based on alignment of target sequence in a comparison range.


