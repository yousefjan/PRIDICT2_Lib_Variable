import re
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
from reportlab.pdfgen import canvas

SF = 'GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'
fivep_homo = 'ATCTTGTGGAAAGGACGAGGTACCG'
threep_homo = 'CGCGGTTCTATCTAGTTACGCGTTA'
codon_table = {
    'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'],
    'S': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
    'Q': ['CAA', 'CAG'],
    'M': ['ATG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCT', 'CCG', 'CCA', 'CCC'],
    'K': ['AAG', 'AAA'],
    'STOP': ['TAG', 'TGA', 'TAA'],
    'T': ['ACC', 'ACA', 'ACG', 'ACT'],
    'F': ['TTT', 'TTC'],
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'G': ['GGT', 'GGG', 'GGA', 'GGC'],
    'I': ['ATC', 'ATA', 'ATT'],
    'L': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
    'H': ['CAT', 'CAC'],
    'R': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
    'W': ['TGG'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'E': ['GAG', 'GAA'],
    'Y': ['TAT', 'TAC']
}

# 5' - N(NG G)NN - 3'
# First codon (N(NG destroy last G)
pam_snv = {'AAG': 'AAA', 'CAG': 'CAA', 'GAG': 'GAA', 'TAG': 'TAG',
           'ACG': 'ACC', 'CCG': 'CCT', 'GCG': 'GCC', 'TCG': 'TCT',
           'ATG': 'ATG', 'CTG': 'CTT', 'GTG': 'GTT', 'TTG': 'TTA',
           'AGG': 'AGA', 'CGG': 'CGT', 'GGG': 'GGC', 'TGG': 'TGG'}
  

def _random_filler() -> str:
    """Returns a 9bp random DNA sequence for filler
    """
    return ''.join(random.choice("ATGC") for _ in range(9))


def _c(seq: str) -> str:
    """Returns the COMPLEMENT of <seq>
    """
    _seq = seq.upper()
    _c_ = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(_c_[bp] for bp in _seq)


def _r(seq: str) -> str:
    """Returns the REVERSE of <seq>
    """
    return seq[::-1]

def _overlap(x: str, y: str, z: str) -> str:
    """Returns <z> (RTT) contained in <y> (saturation area) which is substring of <x> (gene)
    """
    # seq,sseq,rtt
    s1, s2 = x.find(y), x.find(z)
    e1, e2 = s1 + len(y) - 1, s2 + len(z) - 1
    if s1 < s2:
        if e1 < s2:
            # no overlap, y ends before z begins
            print(z)
            raise ValueError("Error in saturation area")
        else:
            if e2 < e1:
                # z is completely contained in y
                return z
            else:
                return x[s2:e1 + 1]
    else:
        if e2 < s1:
            # no overlap, y ends before z begins
            # print(z)
            # raise ValueError("Error in saturation area")
            return z
        else:
            if e1 < e2:
                # z is completely contained in y
                return y
            else:
                return x[s1:e2 + 1]


def trim_string(seq, sseq) -> str:
    """Returns a subsequence of <seq> that accomodates RTTs with at least 1 bp in <sseq>
    
    MAINTAINS FRAME
    """
    # returns substring of seq in same frame
    index = seq.find(sseq)
    start = 0

    while start < index-41:
        start += 3

    end = min(len(seq), index + len(sseq) + 8*3)
    trimmed_string = seq[start:end]

    return trimmed_string


def _find_rtt(seq, sseq, strand) -> list:
    """Returns a list of RTTs in <seq> for the <strand> strand
    """
    pam = r'(?=(.{20}GG.{19}))'
    rtts = []

    if strand == "+":
        seq = trim_string(seq, sseq)
        pam_matches = re.findall(pam, seq)
        for match in pam_matches:
            rtt = _r(_c(match[16:]))
            rtts.append(rtt)

    if strand == "-":
        seq, sseq = _r(_c(seq)), _r(_c(sseq))
        seq = trim_string(seq, sseq)

        pam_matches = re.findall(pam, seq)
        for match in pam_matches:
            rtt = _r(_c(match[16:]))
            rtts.append(rtt)

    if strand is None:
        seq = trim_string(seq, sseq)
        pam_matches = re.findall(pam, seq)
        for match in pam_matches:
            rtt = _r(_c(match[16:]))
            rtts.append(rtt)
            
        seq, sseq = _r(_c(seq)), _r(_c(sseq))
        seq = trim_string(seq, sseq)
        
        pam_matches = re.findall(pam, seq)
        for match in pam_matches:
            rtt = _r(_c(match[16:]))
            rtts.append(rtt)

    return rtts


def _get_control_rtt(seq, sseq, rtt, frame, strand, syn, splice) -> str:
    frame -= 1

    if strand == "+":
        rtt_ = _r(_c(rtt))
        idx = (seq.index(rtt_) + 1) % 3
        idx = 1 if idx == 0 else 0 if idx == 1 else 2
        rtt_start = idx + frame

        if rtt_start == 4:
            rtt_start = 1

        if syn:
            rtt_ = _get_synony_rtt(seq, sseq, rtt, frame+1, strand, splice) ####
            if rtt_:
                rtt_ = _r(_c(rtt_[0]))

        codons_rtt = split_into_codons(rtt_, rtt_start)
        codons_rtt[3 if len(codons_rtt[0]) < 3 else 2] = 'TGA'

        ctl_rtt = _r(_c(''.join(codons_rtt)))

        return ctl_rtt

    else: # strand=="-"

        # if syn and strand=='-': # Get seqs_r handles rc of seq
        #     seq = _r(_c(seq))
        #     rtt = _r(_c(rtt))

        idx = (seq.index(rtt) + 1) % 3
        idx = 1 if idx == 0 else 0 if idx == 1 else 2
        rtt_start = idx + frame

        if rtt_start == 4:
            rtt_start = 1

        if syn:
            rtt_ = _get_synony_rtt(seq, sseq, rtt, frame+1, strand, splice) ####
            if rtt_:
                rtt = rtt_[0]

        codons_rtt = split_into_codons(rtt, rtt_start)
        codons_rtt[-4 if len(codons_rtt[-1]) < 3 else -3] = 'TGA'

        return ''.join(codons_rtt)


def _find_rtts(seq: str, rtts: list, sseq: str, frame, syn, strand, splice) -> dict:
    """Returns a dict where keys are to WT RTTs in <rtts> and values are list of every correponding RTT for
    nucleotide level saturation
    """
    seq = trim_string(seq, sseq)
    sseq_end = seq.index(sseq) + len(sseq)

    bps = ['A', 'C', 'T', 'G']
    all_rtts = []
    for rtt in rtts:
        r = ''
        if strand == '+':
            overlap = _r(_c(_overlap(seq, sseq, _r(_c(rtt)))))  # RTT in saturation area

            # Get portion of RTT not in saturation area
            if seq.index(_r(_c(rtt))) + 25 > sseq_end:  # RTT hangs over end of sseq
                overlap_index = rtt.rfind(overlap)
                r = rtt[:overlap_index]
            
            if seq.index(_r(_c(rtt))) < seq.index(sseq):  # RTT hangs over start of sseq
                r = rtt.replace(overlap, '', 1)
        if strand == '-':
            overlap = _overlap(seq, sseq, rtt)

            # RTT not in saturation area
            if seq.index(rtt) + 25 > sseq_end:  # RTT hangs over end of sseq
                r = rtt.replace(overlap, '', 1)
            
            if seq.index(rtt) < seq.index(sseq):  # RTT hangs over start of sseq
                overlap_index = rtt.rfind(overlap)
                r = rtt[:overlap_index]

        one_rtt = []
        for pos in range(len(overlap)):
            for bp in bps:
                if r + overlap in rtts:
                    one_rtt.append(r + overlap[:pos] + bp + overlap[pos + 1:])
                else:
                    one_rtt.append(overlap[:pos] + bp + overlap[pos + 1:] + r)


        one_rtt.append(rtt[:9] + rtt[10:])
        one_rtt.append(_get_control_rtt(seq, sseq, rtt, frame, strand, syn, splice)) ###
        all_rtts.append(one_rtt)

    for i in range(len(rtts)):
        all_rtts[i] = [item for item in all_rtts[i] if item not in rtts]

    d = {}
    for i, rtt in enumerate(rtts):
        d[rtt] = all_rtts[i]

    return d


def _find_pams(seq: str, sseq, strand, n=0) -> list:
    """Returns list of PAMs in <seq> for the <strand> strand

    (-) strand PAM numbering starts at last (+) strand PAM

    ex: [('CGG', 1, '+'), ('TGG', 2, '+') ... ]
    """
    if strand == "+":
        seq = trim_string(seq, sseq)
    
    if strand == "-":
        seq, sseq = _r(_c(seq)), _r(_c(sseq))
        seq = trim_string(seq, sseq)

    pattern = r'(?=(.{19})([ATGC])GG.{19})'
    matches = re.finditer(pattern, seq)

    sequences = []

    for match in matches:
        sequence = match.group(2) + 'GG'
        n += 1
        sequences.append((sequence, n, strand))

    return sequences


def _find_seqs(seq: str, p: list, sseq: str, frame: int, splice, strand: str) -> dict:
    """Returns entire peg combination: a dict with keys as PAMs and values as list of corresponding peg sequence

    ex: {('CGG', 1, '+'): ['spacer', 'pbs', ['rtt1', 'rtt2', ..., 'rtt75']], ('TGG', 2, '+'): ... }
    """
    if strand == "+":
        seq = trim_string(seq, sseq)
    else:
        seq, sseq = _r(_c(seq)), _r(_c(sseq))
        seq = trim_string(seq, sseq)

    pam = r'(?=(.{20}GG.{19}))'
    pam_matches = re.findall(pam, seq)

    spacers = []
    pbss = []
    for match in pam_matches:
        if strand == "+":
            spacer = 'G' + match[:-22]
        else:
            spacer = 'G' + match[0:19]
        spacers.append(spacer)
        pbs = _r(_c(match[3:16]))
        pbss.append(pbs)

    rtts = _find_rtt(seq, sseq, strand=strand)
    all_rtts = _find_rtts(seq, rtts, sseq, frame, syn=False, strand=strand, splice=splice)

    seqs = {}
    for i in range(len(p)):
        if spacers[i]:
            seqs[p[i]] = [spacers[i], pbss[i], all_rtts[rtts[i]]]

    return seqs


# TODO: Change to a single _make_df()
def _make_df(seqs: dict) -> pd.DataFrame:

    data = [(key[0], key[1], key[2], value[0], value[1], value[2]) for
            key, value in seqs.items()]

    df = pd.DataFrame(data, columns=['PAM', 'PAM No.', 'Strand', 'Spacer',
                                     'PBS', 'RTTs'])
    df_ = df.explode('RTTs').reset_index(drop=True)
    df_['Scaffold'] = pd.Series([SF] * len(df_))
    df_['Complete epegRNA'] = \
        df_['Spacer'] + df_['Scaffold'] + df_['RTTs'] + df_['PBS']
    df_['Length (bp)'] = df_['Complete epegRNA'].apply(lambda x: len(x))

    return df_
def _make_df_homo(seqs: dict) -> pd.DataFrame:
    df_ = _make_df(seqs)
    df_['LHA'] = fivep_homo
    df_['RHA'] = threep_homo
    df_['Complete epegRNA'] = \
        df_['LHA'] + df_['Spacer'] + df_['Scaffold'] + df_['RTTs'] + df_[
            'PBS'] + df_['RHA']
    df_['Length (bp)'] = df_['Complete epegRNA'].apply(lambda x: len(x))

    return df_
def _make_df_cloning(seq, sseq, seqs: dict) -> pd.DataFrame:

    wt_rtts_f = _find_rtt(seq, sseq, '+')
    wt_rtts_r = _find_rtt(seq, sseq, '-')
    wt_rtts = wt_rtts_f + wt_rtts_r

    df_ = _make_df_homo(seqs)
    df_.drop(columns='Scaffold', inplace=True)
    df_.drop(columns='Complete epegRNA', inplace=True)

    df_['epeg No.'] = df_.index + 1
    df_['Filler'] = df_.apply(lambda row: 'GTTTCGAGACG' + _random_filler() +
                                          'CGTCTCGGTGC', axis=1)
    df_['Complete epegRNA'] = \
        df_['LHA'] + df_['Spacer'] + df_['Filler'] + df_['RTTs'] + df_['PBS'] + \
        df_['RHA']
    df_['Length (bp)'] = df_['Complete epegRNA'].apply(lambda x: len(x))
    df_['Complete epegRNA (SF)'] = \
        df_['LHA'] + df_['Spacer'] + SF + df_['RTTs'] + df_['PBS'] + df_['RHA']
    df_['Length (bp) (SF)'] = df_['Complete epegRNA (SF)'].apply(lambda x: len(x))

    df_['epeg No. within PAM'] = df_.groupby('PAM No.').cumcount() + 1

    df_['Edited DNA Sequence'] = df_.apply(lambda row: process_row(row, seq, wt_rtts), axis=1)

    return df_


def _make_df_freq(seq: str, rtts: list) -> pd.DataFrame:
    """Makes dataframe for frequency table and plot
    """
    counts = [0] * len(seq)

    for rtt in rtts:
        rtt_c = _r(_c(rtt))
        for current_rtt in [rtt, rtt_c]:
            start = 0
            while start < len(seq):
                index = seq.find(current_rtt, start)
                if index == -1:
                    break

                for i in range(len(current_rtt)):
                    counts[index + i] += 1

                start = index + 1

    data = [(index + 1, seq[index], counts[index]) for index in range(len(seq))]

    df = pd.DataFrame(data, columns=['Position', 'Nucleotide', 'Frequency'])
    return df


# For getting reference sequence
# TODO: Combine these functions
def highlight_differences(seq, rtt, start, strand):
    diff = []
    diff.append(seq[:start])
    if len(rtt) == 25:
        for i in range(len(rtt)):
            if seq[start + i] != rtt[i]:
                diff.append(rtt[i].lower())
            else:
                diff.append(seq[start + i])

        diff.append(seq[start+len(rtt):])
        return ''.join(diff)
    else:
        if strand == '+':
            return seq[:start+15] + '-' + seq[start+16:]
        else:
            return seq[:start+9] + '-' + seq[start+10:]
def process_row(row, seq, wt_rtts):
    if row['Strand'] == '+':
        rtt = _r(_c((row['RTTs'])))
        wt_rtt = _r(_c(wt_rtts[row['PAM No.']-1]))
    else:
        rtt = row['RTTs']
        wt_rtt = wt_rtts[row['PAM No.']-1]

    start = seq.find(wt_rtt)
    return highlight_differences(seq, rtt, start, row['Strand'])


def run_cloning(seq: str, sseq: str, frame, splice) -> tuple:
    """Function to run in web to generate epegRNA saturation libraries.

    Homology arms included, epegRNA with SF and filler included, edited sequence included.

    Returns library, library without controls, library with only controls
    """

    seq, sseq = seq.upper(), sseq.upper()

    seqs_f = _find_seqs(seq, p=_find_pams(seq, sseq, '+'), sseq=sseq, frame=frame, splice=splice, strand='+')
    seqs_r = _find_seqs(seq, p=_find_pams(seq, sseq, '-', len(_find_pams(seq, sseq, '+'))), sseq=sseq, frame=frame, splice=splice, strand='-')
    
    df_f = _make_df_cloning(seq, sseq, seqs_f)
    df_r = _make_df_cloning(seq, sseq, seqs_r)
    
    df = pd.concat([df_f, df_r], ignore_index=True)

    df['epeg No.'] = df.index + 1

    df = df[['epeg No.', 'PAM', 'PAM No.', 'Strand', 'LHA', 'Spacer', 'Filler',
             'RTTs', 'PBS', 'RHA', 'Complete epegRNA', 'Length (bp)', 'Complete epegRNA (SF)',
             'Length (bp) (SF)', 'Edited DNA Sequence', 'epeg No. within PAM']]

    df_only_ctl = df.groupby('PAM No.').tail(2)
    df_no_ctl = df.drop(df_only_ctl.index)

    df_only_ctl = df_only_ctl.set_index('epeg No.').reset_index(drop=False)
    df_no_ctl = df_no_ctl.set_index('epeg No.').reset_index(drop=False)

    return df, df_no_ctl, df_only_ctl


def run_cloning_(seq: str, sseq: str, frame, splice) -> tuple:
    """Function to run in web to generate epegRNA saturation libraries.

    Homology arms NOT included, epegRNA with SF and filler included, edited sequence included.

    Returns library, library without controls, library with only controls
    """

    df, df_no_ctl, df_only_ctl = run_cloning(seq, sseq, frame, splice)

    df.drop(columns=['LHA', 'RHA'], inplace=True)
    df_only_ctl.drop(columns=['LHA', 'RHA'], inplace=True)
    df_no_ctl.drop(columns=['LHA', 'RHA'], inplace=True)

    df['Complete epegRNA'] = df['Spacer'] + df['Filler'] + df['RTTs'] + df['PBS']
    df['Length (bp)'] = df['Complete epegRNA'].apply(lambda x: len(x))
    df['Complete epegRNA (SF)'] = df['Spacer'] + SF + df['RTTs'] + df['PBS']
    df['Length (bp) (SF)'] = df['Complete epegRNA (SF)'].apply(lambda x: len(x))
    df['Length (bp)'] = df['Complete epegRNA (SF)'].apply(lambda x: len(x))

    df_only_ctl['Complete epegRNA'] = df_only_ctl['Spacer'] + df_only_ctl['Filler'] + df_only_ctl['RTTs'] + df_only_ctl['PBS']
    df_only_ctl['Length (bp)'] = df_only_ctl['Complete epegRNA'].apply(lambda x: len(x))
    df_only_ctl['Complete epegRNA (SF)'] = df_only_ctl['Spacer'] + SF + df_only_ctl['RTTs'] + df_only_ctl['PBS']
    df_only_ctl['Length (bp) (SF)'] = df_only_ctl['Complete epegRNA (SF)'].apply(lambda x: len(x))
    df_only_ctl['Length (bp)'] = df_only_ctl['Complete epegRNA (SF)'].apply(lambda x: len(x))

    df_no_ctl['Complete epegRNA'] = df_no_ctl['Spacer'] + df_no_ctl['Filler'] + df_no_ctl['RTTs'] + df_no_ctl['PBS']
    df_no_ctl['Length (bp)'] = df_no_ctl['Complete epegRNA'].apply(lambda x: len(x))
    df_no_ctl['Complete epegRNA (SF)'] = df_no_ctl['Spacer'] + SF + df_no_ctl['RTTs'] + df_no_ctl['PBS']
    df_no_ctl['Length (bp) (SF)'] = df_no_ctl['Complete epegRNA (SF)'].apply(lambda x: len(x))
    df_no_ctl['Length (bp)'] = df_no_ctl['Complete epegRNA (SF)'].apply(lambda x: len(x))

    return df, df_no_ctl, df_only_ctl


def run_freq_table(seq: str, sseq: str) -> pd.DataFrame:
    """Function to run in web to generate frequency table
    """
    seq = seq.upper()
    sseq = sseq.upper()
    seq = trim_string(seq, sseq)

    df = _make_df_freq(seq, _find_rtt(seq, sseq=sseq, strand=None))
    start = seq.index(sseq)
    end = start + len(sseq) - 1
    df_ = df[start:end + 1]
    df_['Position'] = np.arange(len(df_)) + 1
    return df_


def run_freq_plot(seq: str, sseq: str) -> None:
    """Function to run in web to generate frequency plot
    """
    seq = seq.upper()
    sseq = sseq.upper()
    seq = trim_string(seq, sseq)

    df = _make_df_freq(seq, _find_rtt(seq, sseq, strand=None))

    start = seq.index(sseq)
    end = start + len(sseq) - 1
    df_ = df[start:end + 1]
    df_['Position'] = np.arange(len(df_)) + 1

    fig, ax = plt.subplots()

    ax.bar(df_['Position'],
           [freq if freq != 0 else -1 for freq in df_['Frequency']])
    ax.set_xlabel('Position')
    ax.set_ylabel('Frequency')
    ax.set_title('Frequency Plot')
    ax.margins(x=0, y=0)
    ax.axhline(y=0, color='r', linestyle='-')

    fig.canvas.draw()
    y_labels = [item.get_text() for item in ax.get_yticklabels()]

    y_labels[0] = '0'
    y_labels[1] = ''
    ax.set_yticklabels(y_labels)

    num_minus_ones = len([freq for freq in df_['Frequency'] if freq == 0])
    minus_one_positions = [pos for freq, pos in
                           zip(df_['Frequency'], df_['Position']) if freq == 0]
    minus_one_text = f'Total 0 Count: {num_minus_ones}\nPositions: {", ".join(map(str, minus_one_positions))}'

    ax.text(1.02, 0.5, minus_one_text, transform=ax.transAxes,
            verticalalignment='center',
            bbox=dict(facecolor='lightgray', alpha=0.5))

    plt.savefig('freq_plot.pdf', bbox_inches='tight')

# TODO: This is incorrect, (-) PAMs not all accounted for
def run_figure(seq, sseq):
    """Function to run in web to generate frequency figure
    """
    seq, sseq = seq.upper(), sseq.upper()
    seq = trim_string(seq, sseq)
    cseq = _c(seq)

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    max_lines = (len(seq) - 1) // 50 + 1
    width, height = 50 + len(seq) * 10, 100 + max_lines * 20

    pams = ['AGG', 'TGG', 'CGG', 'GGG']

    c = canvas.Canvas('pam_figure.pdf', pagesize=(width, height))

    x, y = 20, height - 50
    df = _make_df_freq(seq, _find_rtt(seq, sseq, None))
    positions = df['Position'].unique()

    for i, position in enumerate(positions):
        freq = df.loc[df['Position'] == position, 'Frequency']
        freq_str = freq.to_string(index=False)
        c.drawString(20 + i * 10, y + 24, freq_str)

    sseq_start = seq.index(sseq)
    sseq_end = sseq_start + len(sseq)

    # (+) strand boxes
    count = 1
    for i, base in enumerate(seq):
        complement_base = complement[base]
        c.drawString(x, y, base)
        c.drawString(x, y - 12, complement_base)
        x += 10

        counting = False
        for pam in pams:
            if seq[
               i:i + 3] == pam and sseq_start - 21 <= i <= sseq_end + 2 and i >= 19 and len(
                seq[i:None]) >= 22:
                c.rect(20 + i * 10, y - 2, 30, 12)
                c.drawString(20 + i * 10, y + 12, str(count))
                counting = True

        if counting is True:
            count += 1

    # (-) strand boxes
    for i in range(len(cseq), -1, -1):
        counting = False
        for pam in pams:
            if cseq[i - 3:i] == _r(
                    pam) and sseq_start - 2 <= i <= sseq_end + 21 and i <= len(
                seq) - 19 and len(seq[None:i]) >= 22:
                c.rect(i * 10 - 10, y - 14, 30, 12)
                c.drawString(i * 10 - 10, y - 24, str(count))
                counting = True

        if counting is True:
            count += 1

    line_x_start = 20 + sseq_start * 10
    line_x_end = 20 + sseq_end * 10
    line_y = y - 30
    c.line(line_x_start, line_y, line_x_end, line_y)
    c.drawString((line_x_start + line_x_end) / 2, y - 45, 'Saturation Area')

    c.save()


def split_into_codons(seq: str, start_frame: int) -> list:
    """Split <seq> into list codons based on <start_frame>
    """
    if start_frame != 0:
        codons = [seq[:start_frame]]
        codons += [seq[j:j + 3] for j in range(start_frame, len(seq), 3)]
    else:
        codons = [seq[j:j + 3] for j in range(0, len(seq), 3)]

    return codons


def is_one_char_different(codon1, codon2):
    """Returns True iff <codon1> and <codon2> differ by one base
    """
    diff_count = 0
    for c1, c2 in zip(codon1, codon2):
        if c1 != c2:
            diff_count += 1
        if diff_count > 1:
            return False
    return diff_count == 1


def get_edit_position(str1, str2, j=0):
    """Returns index of edit position in <str1>"""
    for i in range(len(str2)):
        if i+j >= len(str1):
            return None # assuming splice sites always in saturation area.
        
        if str1[i+j] != str2[i]:
            return i+j
    return 0


#TODO: Do this in a non retarded way
def _get_synony_rtt(seq: str, sseq: str, rtt: str, frame: int, strand: str, splice: list) -> list:
    """Returns list of RTTs installing silent mutations based off WT <rtt>. Avoids splice sites.
    """
    
    def find_synony_codon(codon: str, left_restrict=0, right_restrict=0, reverse_order=0):
        """Returns priority silent mutation installing codon
        """
        for _, codons in codon_table.items():
            if codon in codons:
                filtered_codons = [c for c in codons if c != codon]

                if right_restrict > 0:
                    filtered_codons = [c for c in filtered_codons if c[-right_restrict:] == codon[-right_restrict:]]
                elif left_restrict > 0:
                    filtered_codons = [c for c in filtered_codons if c[:left_restrict] == codon[:left_restrict]]

                one_char_diff_codons = [c for c in filtered_codons if is_one_char_different(c, codon)]

                sorted_codons = sorted(one_char_diff_codons, key=lambda x: [x[i] != codon[i] for i in range(len(x))], reverse=bool(reverse_order))

                if sorted_codons:
                    return sorted_codons[0]

        return None


    def find_upstream_codons(seq_codons, rtt_codons, downstream=True):
        """Finds codons in <seq_codons> that are upstream (in seed region) of PAM corresponding to RTT <rtt_codons>
        
        For reverse complement RTT (a (+) strand PAM) go downstream.
        """
        seq_str = ''.join(seq_codons)
        rtt_str = ''.join(rtt_codons)

        start_index = seq_str.find(rtt_str)
        if start_index == -1:
            return []

        codon_start_index = start_index // 3
        codon_end_index = (start_index + len(rtt_str)) // 3

        if not downstream:  # (+)
            return seq_codons[max(0, codon_start_index):codon_start_index+2]
        else:  # (-)
            return seq_codons[codon_end_index-1:codon_end_index+1]
    
    seq = seq.upper()
    frame -= 1
    synony_rtts=[]

    if strand == "+":
        rtt_ = _r(_c(rtt))
        idx = (seq.index(rtt_) + 1) % 3
        idx = 1 if idx == 0 else 0 if idx == 1 else 2
        rtt_start = idx + frame
        if rtt_start == 4:
            rtt_start = 1

        codons_rtt = split_into_codons(rtt_, rtt_start)
        codons_seq = split_into_codons(seq, frame)
        upstream = find_upstream_codons(codons_seq, codons_rtt, downstream=False)

        # PAM is +3, 5'-NNN N(NG G)NN-3', 'G' in G)NN codon never destroyed
        if len(codons_rtt[0])==2: 
            if codons_rtt[1] != pam_snv[codons_rtt[1]]:
                new_rtt = rtt[:-5] + _r(_c(pam_snv[codons_rtt[1]])) + rtt[-2:]
                synony_rtts.append(new_rtt)

            # check for synonymous mutations in upstream bases
            codon = find_synony_codon(codons_rtt[1], right_restrict=1)
            if codon is not None:
                new_rtt = rtt[:-5] + _r(_c(codon)) + rtt[-2:]
                synony_rtts.append(new_rtt)

            codon = find_synony_codon(upstream[1], left_restrict=2)
            if codon is not None:
                new_rtt = rtt[:-2] + _r(_c(codon))[:2]
                synony_rtts.append(new_rtt)

            # Remove synony RTTs that edit splice site
            count = 0
            rtt_ = _overlap(seq, sseq, _r(_c(rtt))) # revese complement of overlap
            for _ in range(len(synony_rtts)):
                synony = _r(_c(synony_rtts[count]))
                if seq.index(_r(_c(rtt))) + len(rtt) > (seq.index(sseq) + len(sseq)):  # RTT hangs over end of sseq
                    index_rtt_in_sseq = len(rtt_)
                    synony_overlap = synony[:index_rtt_in_sseq]
                    pos = get_edit_position(sseq, synony_overlap, sseq.rfind(rtt_))
                else:
                    index_rtt_in_sseq = _r(_c(rtt)).rfind(rtt_) 
                    synony_overlap = synony[index_rtt_in_sseq:]
                    pos = get_edit_position(sseq, synony_overlap, sseq.index(rtt_))
                
                if pos in splice:
                    synony_rtts.pop(count)
                else:
                    count += 1

            # downstream of PAM.
            i = 0
            while len(synony_rtts) < 4:
                codon = find_synony_codon(codons_rtt[2+i], reverse_order=1)
                if codon is not None:
                    new_rtt = rtt[:-8-3*i] + _r(_c(codon)) + rtt[-5-3*i:]
                    synony_rtts.append(new_rtt)
                i += 1
            
            # Remove synony RTTs that edit splice site
            count = 0
            rtt_ = _overlap(seq, sseq, _r(_c(rtt))) # revese complement of overlap
            for _ in range(len(synony_rtts)):
                synony = _r(_c(synony_rtts[count]))
                if seq.index(_r(_c(rtt))) + len(rtt) > (seq.index(sseq) + len(sseq)):  # RTT hangs over end of sseq
                    index_rtt_in_sseq = len(rtt_)
                    synony_overlap = synony[:index_rtt_in_sseq]
                    pos = get_edit_position(sseq, synony_overlap, sseq.rfind(rtt_))
                else:
                    index_rtt_in_sseq = _r(_c(rtt)).rfind(rtt_) 
                    synony_overlap = synony[index_rtt_in_sseq:]
                    pos = get_edit_position(sseq, synony_overlap, sseq.index(rtt_))
                
                if pos in splice:
                    synony_rtts.pop(count)
                else:
                    count += 1


            return synony_rtts[:min(len(synony_rtts), 2)]

        # if PAM is +1 5'-NNN (NGG)-3'
        elif len(codons_rtt[0])==3:
            if codons_rtt[1] == 'AGG':
                new_rtt = rtt[:-6] + _r(_c('AGA')) + rtt[-3:]
                synony_rtts.append(new_rtt)
            elif codons_rtt[1] == 'CGG':
                new_rtt = rtt[:-6] + _r(_c('AGA')) + rtt[-3:]
                synony_rtts.append(new_rtt)
            elif codons_rtt[1] == 'GGG':
                new_rtt = rtt[:-6] + _r(_c('AGA')) + rtt[-3:]
                synony_rtts.append(new_rtt)

            # check for synonymous mutations in upstream bases - change NNN in NNN TGG
            codon = find_synony_codon(codons_rtt[0])
            if codon is not None:
                new_rtt = rtt[:-3] + _r(_c(codon))
                synony_rtts.append(new_rtt)

            codon = find_synony_codon(codons_rtt[0], left_restrict=2)
            if codon is not None:
                new_rtt = rtt[:-3] + _r(_c(codon))
                if new_rtt not in synony_rtts:
                    synony_rtts.append(new_rtt)

            codon = find_synony_codon(codons_rtt[0], left_restrict=1)
            if codon is not None:
                new_rtt = rtt[:-3] + _r(_c(codon))
                if new_rtt not in synony_rtts:
                    synony_rtts.append(new_rtt)

            codon = find_synony_codon(codons_rtt[0])
            if codon is not None:
                new_rtt = rtt[:-3] + _r(_c(codon))
                if new_rtt not in synony_rtts:
                    synony_rtts.append(new_rtt)

            # Remove synony RTTs that edit splice site
            count = 0
            rtt_ = _overlap(seq, sseq, _r(_c(rtt))) # revese complement of overlap
            for _ in range(len(synony_rtts)):
                synony = _r(_c(synony_rtts[count]))
                if seq.index(_r(_c(rtt))) + len(rtt) > (seq.index(sseq) + len(sseq)):  # RTT hangs over end of sseq
                    index_rtt_in_sseq = len(rtt_)
                    synony_overlap = synony[:index_rtt_in_sseq]
                    pos = get_edit_position(sseq, synony_overlap, sseq.rfind(rtt_))
                else:
                    index_rtt_in_sseq = _r(_c(rtt)).rfind(rtt_) 
                    synony_overlap = synony[index_rtt_in_sseq:]
                    pos = get_edit_position(sseq, synony_overlap, sseq.index(rtt_))
                
                if pos in splice:
                    synony_rtts.pop(count)
                else:
                    count += 1

            # downstream of PAM
            i = 0
            while len(synony_rtts) < 4:
                codon = find_synony_codon(codons_rtt[2+i], reverse_order=1)
                if codon is not None:
                    new_rtt = rtt[:-9-3*i] + _r(_c(codon)) + rtt[-6-3*i:]
                    if new_rtt not in synony_rtts:
                        synony_rtts.append(new_rtt)
                i += 1

            # Remove synony RTTs that edit splice site
            count = 0
            rtt_ = _overlap(seq, sseq, _r(_c(rtt))) # revese complement of overlap
            for _ in range(len(synony_rtts)):
                synony = _r(_c(synony_rtts[count]))
                if seq.index(_r(_c(rtt))) + len(rtt) > (seq.index(sseq) + len(sseq)):  # RTT hangs over end of sseq
                    index_rtt_in_sseq = len(rtt_)
                    synony_overlap = synony[:index_rtt_in_sseq]
                    pos = get_edit_position(sseq, synony_overlap, sseq.rfind(rtt_))
                else:
                    index_rtt_in_sseq = _r(_c(rtt)).rfind(rtt_) 
                    synony_overlap = synony[index_rtt_in_sseq:]
                    pos = get_edit_position(sseq, synony_overlap, sseq.index(rtt_))
                
                if pos in splice:
                    synony_rtts.pop(count)
                else:
                    count += 1

            return synony_rtts[:min(len(synony_rtts), 2)]

        # PAM is +2, 5'-N NN(N GG)N-3'
        else:
            codon = find_synony_codon(upstream[1], left_restrict=2)
            if codon is not None:
                new_rtt = rtt[:-4] + _r(_c(codon)) + rtt[-1]
                synony_rtts.append(new_rtt)

            codon = find_synony_codon(upstream[1], left_restrict=1)
            if codon is not None:
                new_rtt = rtt[:-4] + _r(_c(codon)) + rtt[-1]
                if new_rtt not in synony_rtts:
                    synony_rtts.append(new_rtt)

            codon = find_synony_codon(upstream[1])
            if codon is not None:
                new_rtt = rtt[:-4] + _r(_c(codon)) + rtt[-1]
                if new_rtt not in synony_rtts:
                    synony_rtts.append(new_rtt)

            codon = find_synony_codon(upstream[0], left_restrict=2)
            if codon is not None:
                new_rtt = rtt[:-1] + _r(_c(codon))[0]
                synony_rtts.append(new_rtt)

            
            # Remove synony RTTs that edit splice site
            count = 0
            rtt_ = _overlap(seq, sseq, _r(_c(rtt))) # revese complement of overlap
            for _ in range(len(synony_rtts)):
                synony = _r(_c(synony_rtts[count]))
                if seq.index(_r(_c(rtt))) + len(rtt) > (seq.index(sseq) + len(sseq)):  # RTT hangs over end of sseq
                    index_rtt_in_sseq = len(rtt_)
                    synony_overlap = synony[:index_rtt_in_sseq]
                    pos = get_edit_position(sseq, synony_overlap, sseq.rfind(rtt_))
                else:
                    index_rtt_in_sseq = _r(_c(rtt)).rfind(rtt_) 
                    synony_overlap = synony[index_rtt_in_sseq:]
                    pos = get_edit_position(sseq, synony_overlap, sseq.index(rtt_))
                
                if pos in splice:
                    synony_rtts.pop(count)
                else:
                    count += 1

            # downstream of PAM
            i = 0
            while len(synony_rtts) < 4:
                codon = find_synony_codon(codons_rtt[2+i])
                if codon is not None:
                    new_rtt = rtt[:-7-3*i] + _r(_c(codon)) + rtt[-4-3*i:]
                    synony_rtts.append(new_rtt)
                i += 1

            # Remove synony RTTs that edit splice site
            count = 0
            rtt_ = _overlap(seq, sseq, _r(_c(rtt))) # revese complement of overlap
            for _ in range(len(synony_rtts)):
                synony = _r(_c(synony_rtts[count]))
                if seq.index(_r(_c(rtt))) + len(rtt) > (seq.index(sseq) + len(sseq)):  # RTT hangs over end of sseq
                    index_rtt_in_sseq = len(rtt_)
                    synony_overlap = synony[:index_rtt_in_sseq]
                    pos = get_edit_position(sseq, synony_overlap, sseq.rfind(rtt_))
                else:
                    index_rtt_in_sseq = _r(_c(rtt)).rfind(rtt_) 
                    synony_overlap = synony[index_rtt_in_sseq:]
                    pos = get_edit_position(sseq, synony_overlap, sseq.index(rtt_))
                
                if pos in splice:
                    synony_rtts.pop(count)
                else:
                    count += 1

            return synony_rtts[:min(len(synony_rtts), 2)]


    if strand == "-":  # Only (+) strand coding (gene for sat.)
        idx = (seq.index(rtt) + 1) % 3
        idx = 1 if idx == 0 else 0 if idx == 1 else 2
        rtt_start = idx + frame

        if rtt_start == 4:
            rtt_start = 1

        codons_rtt = split_into_codons(rtt, rtt_start)
        codons_seq = split_into_codons(seq, frame)
        upstream = find_upstream_codons(codons_seq, codons_rtt)

        # PAM is +3, 5'-NN(C CN)N NN-3'
        if len(codons_rtt[-1])==2:

            # NN(C cN)N NN
            codon = find_synony_codon(codons_rtt[-2], reverse_order=1, right_restrict=2)
            if codon is not None:
                new_rtt = rtt[:-5] + codon + rtt[-2:]
                synony_rtts.append(new_rtt)

            # NN(c CN)N NN
            codon = find_synony_codon(codons_rtt[-3], reverse_order=1, left_restrict=2)
            if codon is not None:
                new_rtt = rtt[:-8] + codon + rtt[-5:]
                synony_rtts.append(new_rtt)

            # NN(C Cn)n NN
            codon = find_synony_codon(codons_rtt[-2], reverse_order=1, left_restrict=1)
            if codon is not None:
                new_rtt = rtt[:-5] + codon + rtt[-2:]
                synony_rtts.append(new_rtt)

            # NN(C CN)N nn
            codon = find_synony_codon(upstream[0], reverse_order=1, right_restrict=1)
            if codon is not None:
                new_rtt = rtt[:-2] + codon[:2]
                synony_rtts.append(new_rtt)

            # Remove synony RTTs that edit splice site
            count = 0
            rtt_ = _overlap(seq, sseq, rtt)
            for _ in range(len(synony_rtts)):
                synony = synony_rtts[count]   
                if seq.index(rtt) + len(rtt) > (seq.index(sseq) + len(sseq)):  # RTT hangs over end of sseq
                    index_rtt_in_sseq = len(rtt_)
                    pos = get_edit_position(sseq, synony[:index_rtt_in_sseq], sseq.index(rtt_))
                else:
                    index_rtt_in_sseq = rtt.rfind(rtt_)
                    pos = get_edit_position(sseq, synony[index_rtt_in_sseq:], sseq.index(rtt_))
                
                if pos in splice:
                    synony_rtts.pop(count)
                else:
                    count += 1

            # downstream of PAM
            i = 0
            while len(synony_rtts) < 2:
                codon = find_synony_codon(codons_rtt[-4-i], reverse_order=1)
                if codon is not None:
                    new_rtt = rtt[:-11-3*i] + codon + rtt[-8-3*i:]
                    synony_rtts.append(new_rtt)
                i += 1

            # Remove synony RTTs that edit splice site
            count = 0
            rtt_ = _overlap(seq, sseq, rtt)
            for _ in range(len(synony_rtts)):
                synony = synony_rtts[count]   
                if seq.index(rtt) + len(rtt) > (seq.index(sseq) + len(sseq)):  # RTT hangs over end of sseq
                    index_rtt_in_sseq = len(rtt_)
                    pos = get_edit_position(sseq, synony[:index_rtt_in_sseq], sseq.index(rtt_))
                else:
                    index_rtt_in_sseq = rtt.rfind(rtt_)
                    pos = get_edit_position(sseq, synony[index_rtt_in_sseq:], sseq.index(rtt_))
                
                if pos in splice:
                    synony_rtts.pop(count)
                else:
                    count += 1

            return synony_rtts[:min(len(synony_rtts), 2)]

        # if PAM is +1, 5'-(CCN) NNN-3'
        elif len(codons_rtt[-1]) == 3 :
            # no direct possible (CCN always P)
            # so check for synonymous mutations in upstream bases

            codon = find_synony_codon(codons_rtt[-1], reverse_order=1, right_restrict=2)
            if codon is not None:
                new_rtt = rtt[:-3] + codon
                if new_rtt not in synony_rtts:
                    synony_rtts.append(new_rtt)

            codon = find_synony_codon(codons_rtt[-1], reverse_order=1, right_restrict=1)
            if codon is not None:
                new_rtt = rtt[:-3] + codon
                if new_rtt not in synony_rtts:
                    synony_rtts.append(new_rtt)

            codon = find_synony_codon(codons_rtt[-1], reverse_order=1)
            if codon is not None:
                new_rtt = rtt[:-3] + codon
                if new_rtt not in synony_rtts:
                    synony_rtts.append(new_rtt)

            # Remove synony RTTs that edit splice site
            count = 0
            rtt_ = _overlap(seq, sseq, rtt)
            for _ in range(len(synony_rtts)):
                synony = synony_rtts[count]   
                if seq.index(rtt) + len(rtt) > (seq.index(sseq) + len(sseq)):  # RTT hangs over end of sseq
                    index_rtt_in_sseq = len(rtt_)
                    pos = get_edit_position(sseq, synony[:index_rtt_in_sseq], sseq.index(rtt_))
                else:
                    index_rtt_in_sseq = rtt.rfind(rtt_)
                    pos = get_edit_position(sseq, synony[index_rtt_in_sseq:], sseq.index(rtt_))
                
                if pos in splice:
                    synony_rtts.pop(count)
                else:
                    count += 1
            
            # Downstream RTTS
            i = 0
            while len(synony_rtts) < 4:
                codon = find_synony_codon(codons_rtt[-3-i])
                if codon is not None:
                    new_rtt = rtt[:-9-3*i] + codon + rtt[-6-3*i:]
                    if new_rtt not in synony_rtts:
                        synony_rtts.append(new_rtt)
                i += 1
            
            # Remove synony RTTs that edit splice site
            count = 0
            rtt_ = _overlap(seq, sseq, rtt)
            for _ in range(len(synony_rtts)):
                synony = synony_rtts[count]   
                if seq.index(rtt) + len(rtt) > (seq.index(sseq) + len(sseq)):  # RTT hangs over end of sseq
                    index_rtt_in_sseq = len(rtt_)
                    pos = get_edit_position(sseq, synony[:index_rtt_in_sseq], sseq.index(rtt_))
                else:
                    index_rtt_in_sseq = rtt.rfind(rtt_)
                    pos = get_edit_position(sseq, synony[index_rtt_in_sseq:], sseq.index(rtt_))
                
                if pos in splice:
                    synony_rtts.pop(count)
                else:
                    count += 1

            return synony_rtts[:min(len(synony_rtts), 2)]

        # PAM is +2, 5'-N(CC N)NN N-3'
        else:
            # 5'-N(Cc N)NN N-3'
            codon = find_synony_codon(codons_rtt[-3], left_restrict=2)
            if codon is not None:
                new_rtt = rtt[:-7] + codon + rtt[-4:]
                synony_rtts.append(new_rtt)

            # 5'-N(cc N)NN N-3'
            codon = find_synony_codon(codons_rtt[-3], left_restrict=1)
            if codon is not None:
                new_rtt = rtt[:-7] + codon + rtt[-4:]
                if new_rtt not in synony_rtts:
                    synony_rtts.append(new_rtt)

            # 5'-N(CC n)NN N-3'
            codon = find_synony_codon(codons_rtt[-2], right_restrict=2)
            if codon is not None:
                new_rtt = rtt[:-4] + codon + rtt[-1]
                synony_rtts.append(new_rtt)

            # 5'-N(CC N)Nn N-3'
            codon = find_synony_codon(codons_rtt[-2], left_restrict=2)
            if codon is not None:
                new_rtt = rtt[:-4] + codon + rtt[-1]
                if new_rtt not in synony_rtts:
                    synony_rtts.append(new_rtt)

            # 5'-N(CC N)nn N-3'
            codon = find_synony_codon(codons_rtt[-2], left_restrict=1)
            if codon is not None:
                new_rtt = rtt[:-4] + codon + rtt[-1]
                if new_rtt not in synony_rtts:
                    synony_rtts.append(new_rtt)


            # 5'-N(CC N)NN n-3'
            codon = find_synony_codon(upstream[1], right_restrict=2)
            if codon is not None:
                new_rtt = rtt[:-1] + codon[0]
                synony_rtts.append(new_rtt)

            # downstream of PAM
            codon = find_synony_codon(codons_rtt[-3])
            if codon is not None:
                new_rtt = rtt[:-7] + codon + rtt[-4:]
                if new_rtt not in synony_rtts:
                    synony_rtts.append(new_rtt)

            
            # Remove synony RTTs that edit splice site
            count = 0
            rtt_ = _overlap(seq, sseq, rtt)
            index_rtt_in_sseq = rtt.index(rtt_)
            for _ in range(len(synony_rtts)):
                synony = synony_rtts[count]   
                if seq.index(rtt) + len(rtt) > (seq.index(sseq) + len(sseq)):  # RTT hangs over end of sseq
                    pos = get_edit_position(sseq, synony[:index_rtt_in_sseq], sseq.index(rtt_))
                else:
                    pos = get_edit_position(sseq, synony[index_rtt_in_sseq:], sseq.index(rtt_))
                
                if pos in splice:
                    synony_rtts.pop(count)
                else:
                    count += 1

            # Downstream RTTs
            i = 0
            while len(synony_rtts) < 4:
                codon = find_synony_codon(codons_rtt[-4-i])
                if codon is not None:
                    new_rtt = rtt[:-10-3*i] + codon + rtt[-7-3*i:]
                    if new_rtt not in synony_rtts:
                        synony_rtts.append(new_rtt)
                i += 1

            # Remove synony RTTs that edit splice site
            count = 0
            rtt_ = _overlap(seq, sseq, rtt)
            index_rtt_in_sseq = rtt.index(rtt_)
            for _ in range(len(synony_rtts)):
                synony = synony_rtts[count]   
                if seq.index(rtt) + len(rtt) > (seq.index(sseq) + len(sseq)):  # RTT hangs over end of sseq
                    pos = get_edit_position(sseq, synony[:index_rtt_in_sseq], sseq.index(rtt_))
                else:
                    pos = get_edit_position(sseq, synony[index_rtt_in_sseq:], sseq.index(rtt_))
                
                if pos in splice:
                    synony_rtts.pop(count)
                else:
                    count += 1

            return synony_rtts

    # No PAM destruction nor upstream base synonymous edit possible
    return None

def get_preserving_rtt(rtt_lst, sat_rtt, wt_rtt):
    for synony_rtt in rtt_lst:
        synony_pos = get_edit_position(synony_rtt, wt_rtt)  # Position of synony edit in the current synony rtt from str_lst
        sat_pos = get_edit_position(sat_rtt, wt_rtt)

        if synony_pos == sat_pos:
            continue  # Try the next synony_rtt

        if isinstance(synony_pos, int):
            sat_rtt = sat_rtt[:synony_pos] + synony_rtt[synony_pos] + sat_rtt[synony_pos+1:]
            return sat_rtt
        else:
            raise ValueError(f"Invalid edit position found {synony_pos}")

    return 'No Alternate Silent Mutations Available'

def run_synony(seq: str, sseq: str, frame: int, df, HA, splice):
    """Run in web to generate complete synonymous epegRNA saturation cloning library

    <df> is no_ctl library from run_cloning or run_cloning_
    """

    seq, sseq = seq.upper(), sseq.upper()

    wt_rtts_f = _find_rtt(seq, sseq, '+')
    wt_rtts_r = _find_rtt(seq, sseq, '-')
    wt_rtts = wt_rtts_f + wt_rtts_r

    # Makes list of strand of each PAM used e.g. ['+', '+', '+', '-', '-']
    strands = [''] * df['PAM No.'].max()
    for _, row in df.iterrows():
        strands[row['PAM No.'] - 1] = row['Strand']

    # Dict of PAM No. to its corresponding synony list
    pams = df['PAM No.'].unique()
    pam_num_to_lst = {}
    for i in range(len(pams)):
        lst = _get_synony_rtt(seq, sseq, wt_rtts[i], frame, strands[i], splice)
        pam_num_to_lst[i+1] = lst

    # Row for each epeg using the dict
    new_rows = []
    for _, row in df.iterrows():
        new_row = row.copy()
        synony_rtts = pam_num_to_lst[row['PAM No.']]
        j = row['PAM No.'] - 1
        pres_rtt = get_preserving_rtt(synony_rtts, row['RTTs'], wt_rtts[j])
        new_row['RTTs'] = pres_rtt
        new_row['Syn. Mutation Position'] = 42-get_edit_position(new_row['RTTs'], row['RTTs'])
        new_row['Edited DNA Sequence (syn)'] = process_row(new_row, seq, wt_rtts)
        new_row['Syn. Mutation Position (seq)'] = get_edit_position(new_row['Edited DNA Sequence'], new_row['Edited DNA Sequence (syn)']) + 1
        new_rows.append(new_row)
    
    new_rows_df = pd.DataFrame(new_rows)

    # Add controls
    dfs = []
    for i, group in new_rows_df.groupby('PAM No.'):
        new_row_stop = group.iloc[0].copy()
        new_row_del = group.iloc[0].copy()

        new_row_stop['RTTs'] = _get_control_rtt(seq, sseq, wt_rtts[i-1], frame, strand=strands[i-1], syn=True, splice=splice)
        new_row_stop['Edited DNA Sequence (syn)'] = process_row(new_row_stop, seq, wt_rtts)
        new_row_stop['Filler'] = 'GTTTCGAGACG' + _random_filler() + 'CGTCTCGGTGC'
        
        new_row_del['RTTs'] = wt_rtts[i-1][:9] + wt_rtts[i-1][10:]
        new_row_del['Syn. Mutation Position'] = None
        new_row_del['Edited DNA Sequence (syn)'] = process_row(new_row_del, seq, wt_rtts)
        new_row_del['Filler'] = 'GTTTCGAGACG' + _random_filler() + 'CGTCTCGGTGC'

        group_ = pd.concat([group, pd.DataFrame([new_row_stop, new_row_del])], ignore_index=True)
        dfs.append(group_)
    
    new_rows_df = pd.concat(dfs, ignore_index=True)

    new_rows_df.drop('Edited DNA Sequence', axis=1, inplace=True)

    
    if HA:
        new_rows_df['Complete epegRNA'] = new_rows_df['LHA'] + new_rows_df['Spacer'] + new_rows_df['Filler'] + new_rows_df['RTTs'] + new_rows_df['PBS'] + new_rows_df['RHA']
        new_rows_df['Complete epegRNA (SF)'] = new_rows_df['LHA'] + new_rows_df['Spacer'] + SF + new_rows_df['RTTs'] + new_rows_df['PBS'] + new_rows_df['RHA']
    else:
        new_rows_df['Complete epegRNA'] = new_rows_df['Spacer'] + new_rows_df['Filler'] + new_rows_df['RTTs'] + new_rows_df['PBS']
        new_rows_df['Complete epegRNA (SF)'] = new_rows_df['Spacer'] + SF + new_rows_df['RTTs'] + new_rows_df['PBS']

    new_rows_df['Length (bp)'] = new_rows_df['Complete epegRNA'].apply(lambda x: len(x))
    new_rows_df['Length (bp) (SF)'] = new_rows_df['Complete epegRNA (SF)'].apply(lambda x: len(x))

    new_rows_df['epeg No.'] = new_rows_df.index+1
    new_rows_df['epeg No. within PAM'] = new_rows_df.groupby('PAM No.').cumcount() + 1

    return new_rows_df
