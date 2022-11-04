from typing import NamedTuple

import numpy as np

from Bio import Align


class NamedArgs(NamedTuple):
    val: str
    name: str
        
    @staticmethod
    def build(arg):
        try:
            name, val = arg.split(':')
        except:
            return NamedArgs(arg,arg)
        return NamedArgs(val, name)
    
    def __str__(self):
        return f'{self.name}'
    def upper(self):
        return NamedArgs(self.val.upper(), self.name)
    @property
    def str(self):
        return self.val


class TargetRegion(NamedTuple):
    aseq: str
    target_seq: str
    beg: int
    end: int
    pam_beg: int
    pam_end: int
    cmp_beg: int
    cmp_end: int

    @staticmethod
    def build(aseq, target_seq, cmp_region_radius, pem_len):
        i = aseq.find(target_seq)
        if i == -1:
            return None
        beg = i
        end = i + len(target_seq)
        pam_beg = end - pem_len
        pam_end = end
        cmp_beg = max(0, pam_beg - cmp_region_radius)
        cmp_end = min(len(aseq), cmp_beg + cmp_region_radius * 2)
        return TargetRegion(aseq, target_seq,
                            beg, end, pam_beg, pam_end, cmp_beg, cmp_end)

    def comparison_region(self):
        return self.aseq[self.cmp_beg:self.cmp_end]

    def left_indicator_seq(self, length):
        return self.comparison_region()[:length]

    def right_indicator_seq(self, length):
        return self.comparison_region()[-length:]

    def get_slice(self, beg, end):
        return self.aseq[beg:end]

    
class UserRegion(NamedTuple):
    beg: int
    end: int
    @staticmethod
    def build(pam_beg, user_region_beg_offset, length):
        beg = max(pam_beg-user_region_beg_offset, 0)
        return UserRegion(beg, beg+length)


seq_trans=str.maketrans("ATGC", "TACG")
def revertedSeq(seq):
    return seq.translate(seq_trans)[::-1]

def get_aligner(open_gap_score=-2.01):
    aligner = Align.PairwiseAligner()
    # aligner.mode='local'
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.1
    # aligner.gap_score = -1.00
    aligner.open_gap_score = open_gap_score
    # dir(aligner)
    return aligner


toCharArray = lambda seq: np.array([seq]).astype('S').view('S1')


def mismatch(seq1, seq2):
    n = min(len(seq1), len(seq2))
    m = max(len(seq1), len(seq2))
    mismatch = np.sum(toCharArray(seq1[:n]) != toCharArray(seq2[:n]))
    return mismatch + m - n


def matchUpto1(target, seq, reverse=False):
    finder=seq.find
    if reverse:
        finder=seq.rfind
    
    half_len = len(target) // 2
    fst = target[:half_len]
    snd = target[half_len:]

    i_fst = finder(fst)
    while (i_fst != -1):
        beg = i_fst
        end = beg + len(target)
        alignedSeq = seq[beg:end]
        if mismatch(alignedSeq, target) < 2:
            return beg
        i_fst = finder(fst, end)

    i_snd = finder(snd)
    while (i_snd != -1):
        beg = i_snd - half_len
        end = beg + len(target)
        alignedSeq = seq[beg:end]
        if mismatch(alignedSeq, target) < 2:
            return beg
        i_snd = finder(snd, end)
    return -1


def get_user_region_query(aligner, query_seq, tr, ur):
    user_seq_beg = max(0, ur.beg - tr.cmp_beg)
    user_seq_end = max(0, ur.end - tr.cmp_beg)

    alignment = aligner.align(tr.comparison_region(), query_seq)[0]
#     beg_mapper = get_beg_mapper(alignment)
#     end_mapper = get_end_mapper(alignment)
#     assert alignment.target == tr.comparison_region()

    [str_ref, str_align, str_read, _] = str(alignment).split("\n")
    if len(str_ref)==tr.cmp_end-tr.cmp_beg :
        beg=user_seq_beg
        end=user_seq_end
    else:
        seq=list(str_ref)
        mapper=np.cumsum(np.array(seq)!='-')-1
        beg=np.argmax(mapper==user_seq_beg)
        end=np.argmax(mapper==user_seq_end)
    
    return str_read[beg:end].replace('-','')

#     return alignment.query[beg_mapper[user_seq_beg]:end_mapper[user_seq_end]]


def align_read_to_user_region(aligner, x, tr, ur):
    query_seq = x.query_seq
    user_seq_beg = max(0, ur.beg - tr.cmp_beg)
    user_seq_end = max(0, ur.end - tr.cmp_beg)

    
    alignment = aligner.align(tr.comparison_region(), query_seq)[0]

    [str_ref, str_align, str_read, _] = str(alignment).split("\n")
    if len(str_ref)==tr.cmp_end-tr.cmp_beg :
        f = lambda x: x[user_seq_beg:user_seq_end]
    else:
        seq=list(str_ref)
        mapper=np.cumsum(np.array(seq)!='-')-1
        beg=np.argmax(mapper==user_seq_beg)
        end=np.argmax(mapper==user_seq_end)
        f = lambda x: x[beg:end]
    # return f(str_ref), f(str_align), f(str_read)
    return {"ref": f(str_ref), "alignment": f(str_align), "read":f(str_read), "n": x.n}

