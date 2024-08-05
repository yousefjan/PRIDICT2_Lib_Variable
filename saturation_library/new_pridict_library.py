import random
import pandas as pd
from Bio.Seq import Seq
from model_pred import pred
from main import SF, fivep_homo, threep_homo, _get_control_rtt, _get_synony_rtt, _random_filler, get_preserving_rtt, _c, get_edit_position

def rc(dna):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    return ''.join(complement[base] for base in reversed(dna))

def get_wt_rtt(seq, rtt): # pass rc of rtt for (+) PAM
    for i, base in enumerate(rtt):
        bases = ['A', 'T', 'G', 'C']
        if base.islower():
            bases.remove(base.upper())
            rtts = [(rtt[:i] + char + rtt[i+1:None]) for char in bases]
            for rtt in rtts:
                if rtt in seq:
                    return rtt

    return None


def generate_strings(seq, satArea):
    '''Returns all PRIDICT inputs for all bases in satArea
    '''
    allSeq = []
    startIndex = seq.index(satArea)
    endIndex = startIndex+len(satArea)
 
    for i in range(len(satArea)):
        bases = ['A','T','C','G']
        bases.remove(satArea[i])
        for j in range(len(bases)):
            editSeq = f'{seq[:startIndex]}{satArea[:i]}({satArea[i]}/{bases[j]}){satArea[i+1:]}{seq[endIndex:]}'
            allSeq.append(editSeq)
 
    return allSeq
 


def run_pridict_lib(seq, sseq):
    """Returns PRIDICT2.0 based saturation library (with variable 3' extension structure)

    Sorted by PAM number (spacer)
    """

    def sort_result(unsorted_rs):
        # unsorted_rs = pd.read_csv('./npc_result.csv', index_col=False)
        pos_df = unsorted_rs[unsorted_rs['Strand']=='(+)'].reset_index(drop=True)
        neg_df = unsorted_rs[unsorted_rs['Strand']=='(-)'].reset_index(drop=True)

        pos_indx_lst = [seq.index(row['Spacer'][1:]) for _, row in pos_df.iterrows()]
        pos_df['Spacer Index'] = pos_indx_lst
        pos_df = pos_df.sort_values('Spacer Index', ignore_index=True)

        neg_indx_lst = [seq.index(rc(row['Spacer'][1:])) for _, row in neg_df.iterrows()]
        neg_df['Spacer Index'] = neg_indx_lst
        neg_df = neg_df.sort_values('Spacer Index', ascending=False, ignore_index=True)
        
        sorted_rs = pd.concat([pos_df, neg_df], ignore_index=True).reset_index(drop=True)

        PAMnum = 1
        old_indx = sorted_rs.iloc[0]['Spacer Index']
        for i, row in sorted_rs.iterrows():
            new_indx = row['Spacer Index']
            if old_indx!=new_indx:
                PAMnum+=1
            sorted_rs.at[i, 'PAM No.'] = PAMnum
            old_indx = new_indx
        
        # sorted_rs.insert(1, 'Spacer Index', sorted_rs.pop('Spacer Index'))
        sorted_rs.drop('Spacer Index', axis=1, inplace=True)

        sorted_rs['PAM No.'] = sorted_rs['PAM No.'].astype(int)
        sorted_rs.insert(1, 'PAM No.', sorted_rs.pop('PAM No.'))
        
        return sorted_rs
    
    pridict_input_sequences = generate_strings(seq, sseq)
    dfs = []

    for i, input in enumerate(pridict_input_sequences):
        pred(input)
        output = pd.read_csv('./predictions/seq_pegRNA_Pridict_full.csv')
        output = output[(output['RTrevcomp'].str.len() < 40) & (output['PBSrevcomp'].str.len() < 15) & (output['PBSrevcomp'].str.len() > 7)]
        max_scores = output.loc[output.groupby('Spacer-Sequence')['PRIDICT2_0_editing_Score_deep_HEK'].idxmax()]
        top_spacers = max_scores.sort_values(by='PRIDICT2_0_editing_Score_deep_HEK', ascending=False).head(4)

        for j in range(len(top_spacers.index)):
            df = pd.DataFrame()
            df['peg No. (within edit)'] = [j+1]
            df['Edit Position (sat. area)'] = [i // 3 +1]
            df['PAM'] = [_c(top_spacers.iloc[j]['RTrevcomp'][-4]) + 'GG']
            df['Strand'] = ['(+)' if top_spacers.iloc[j]['Target-Strand'] == 'Fw' else '(-)']
            df['Edit'] = [f'{top_spacers.iloc[j]['OriginalAllele']}>{top_spacers.iloc[j]['EditedAllele']}']
            df['LHA'] = [fivep_homo]
            df['Spacer'] = [top_spacers.iloc[j]['Spacer-Sequence']]
            df['RTTs'] = [top_spacers.iloc[j]['RTrevcomp']]
            df['PBS'] = [top_spacers.iloc[j]['PBSrevcomp']]
            df['RHA'] = [threep_homo]
            df['Filler'] = 'GTTTCGAGACG' + _random_filler() + 'CGTCTCGGTGC'

            df['Complete epegRNA'] = [fivep_homo + top_spacers.iloc[j]['Spacer-Sequence'] + df['Filler'] + top_spacers.iloc[j]['RTrevcomp'] + top_spacers.iloc[j]['PBSrevcomp'] + threep_homo]
            df['Length (bp)'] = df['Complete epegRNA'].str.len()

            df['Complete epegRNA (SF)'] = [fivep_homo + top_spacers.iloc[j]['epegRNA'] + threep_homo] # Uses scaffold
            df['Length (bp) (SF)'] = df['Complete epegRNA (SF)'].str.len()

            df['Reference Sequence'] = [top_spacers.iloc[j]['wide_mutated_target']]
            df['PRIDICT2.0 Score'] = [top_spacers.iloc[j]['PRIDICT2_0_editing_Score_deep_HEK']]
            dfs.append(df)

    unsorted_lib = pd.concat(dfs, ignore_index=True)
    sorted_lib = sort_result(unsorted_lib)

    return sorted_lib


def run_pridict_library_synony(seq, sseq, frame, HA, splice):
    """Returns PRIDICT2.0 based saturation library (with variable 3' extension structure) containing 
    silent mutation installing RTTs

    Sorted by PAM number (spacer)
    """

    def get_synony_reference(rtt, wt_rtt, strand):
        """Returns reference sequence for PRIDICT2.0 based synony saturation library
        """
        if strand == '+':
            rtt = rc(rtt)

        start = seq.find(wt_rtt)

        diff = []
        diff.append(seq[:start])
        if len(rtt) == len(wt_rtt): # not 1bp del ctl
            for i in range(len(rtt)):
                if seq[start + i] != rtt[i]:
                    diff.append(rtt[i].lower())
                else:
                    diff.append(seq[start + i])

            diff.append(seq[start+len(rtt):])
            return ''.join(diff)
        else:
            del_idx = int(len(wt_rtt)/2)
            if strand == '+':
                return seq[:start+len(wt_rtt)-9-1] + '-' + seq[start+len(wt_rtt)-9:]
            else:
                return seq[:start+del_idx] + '-' + seq[start+del_idx+1:]


    # df = run_pridict_lib(seq, sseq)
    df = pd.read_csv('./saturation_library/npc_result.csv')

    # Get WT RTT for controls
    wt_rtts = []
    for _, group in df.groupby('PAM No.'):
        ctl_row = group.loc[group['PRIDICT2.0 Score'].idxmax()].copy()
        if ctl_row['Strand']=='(-)':
            wt_rtt = get_wt_rtt(seq, ctl_row['RTTs'])
        else:
            wt_rtt= get_wt_rtt(seq, rc(ctl_row['RTTs']))

        wt_rtts.append(wt_rtt)

    rows = []
    for _, row in df.iterrows():
        row_syn = row.copy()
        if row['Strand']=='(-)':
            wt_rtt = get_wt_rtt(seq, row['RTTs'])
        else:
            wt_rtt= get_wt_rtt(seq, rc(row['RTTs']))

        synony = _get_synony_rtt(seq=seq, sseq=sseq, rtt=wt_rtt, frame=frame, strand=row['Strand'][1], splice=splice)
        row_syn['RTTs'] = get_preserving_rtt(synony, row['RTTs'].upper(), wt_rtt)

        row_syn['Complete epegRNA'] = [fivep_homo + row['Spacer'] + 'GTTTCGAGACG' + _random_filler() + 'CGTCTCGGTGC' + row['RTTs'] + row['PBS'] + threep_homo]
        row_syn['Complete epegRNA (SF)'] = [fivep_homo + row['Spacer'] + SF + row['RTTs'] + row['PBS'] + threep_homo]
        row_syn['Syn. Mutation Position'] = 42-get_edit_position(row_syn['RTTs'], row['RTTs'].upper())
        row_syn['Reference Sequence'] = get_synony_reference(row_syn['RTTs'].upper(), wt_rtt, row['Strand'][1])

        rows.append(row_syn)

    df = pd.DataFrame(rows)

    # Getting control RTTs
    groups = []
    for i, group in df.groupby('PAM No.'):
        new_row_stop = group.loc[group['PRIDICT2.0 Score'].idxmax()].copy()
        new_row_stop['RTTs'] = _get_control_rtt(seq, sseq, wt_rtts[i-1], frame, new_row_stop['Strand'][1], True, splice)
        new_row_stop['Filler'] = 'GTTTCGAGACG' + _random_filler() + 'CGTCTCGGTGC'
        new_row_stop['Reference Sequence'] = get_synony_reference(new_row_stop, seq, wt_rtt)

        new_row_del = group.loc[group['PRIDICT2.0 Score'].idxmax()].copy()
        del_idx = int(len(wt_rtts[i-1])/2)
        new_row_del['RTTs'] = wt_rtts[i-1][:del_idx] + wt_rtts[i-1][del_idx+1:]
        new_row_del['Filler'] = 'GTTTCGAGACG' + _random_filler() + 'CGTCTCGGTGC'
        new_row_del['Reference Sequence'] = get_synony_reference(new_row_del, seq, wt_rtt)

        group_ = pd.concat([group, pd.DataFrame([new_row_stop, new_row_del])], ignore_index=True)
        groups.append(group_)

    new_rows_df = pd.concat(groups, ignore_index=True)

    if HA:
        new_rows_df['Complete epegRNA'] = new_rows_df['LHA'] + new_rows_df['Spacer'] + new_rows_df['Filler'] + new_rows_df['RTTs'] + new_rows_df['PBS'] + new_rows_df['RHA']
        new_rows_df['Complete epegRNA (SF)'] = new_rows_df['LHA'] + new_rows_df['Spacer'] + SF + new_rows_df['RTTs'] + new_rows_df['PBS'] + new_rows_df['RHA']
    else:
        new_rows_df['Complete epegRNA'] = new_rows_df['Spacer'] + new_rows_df['Filler'] + new_rows_df['RTTs'] + new_rows_df['PBS']
        new_rows_df['Complete epegRNA (SF)'] = new_rows_df['Spacer'] + SF + new_rows_df['RTTs'] + new_rows_df['PBS']

    new_rows_df['Length (bp)'] = new_rows_df['Complete epegRNA'].apply(lambda x: len(x))
    new_rows_df['Length (bp) (SF)'] = new_rows_df['Complete epegRNA (SF)'].apply(lambda x: len(x))

    return new_rows_df


def verify_spacer_in_seq(seq):
    df = pd.read_csv('./saturation_library/npc_result.csv')
    no_match = {}

    for _, row in df.iterrows():
        spacer_without_g = row['Spacer'][1:]
        if spacer_without_g not in seq:
            if str(Seq(spacer_without_g).reverse_complement()) not in seq:
                if row['Edit Position (sat. area)'] not in no_match:
                    no_match[row['Edit Position (sat. area)']] = [(row['Edit'], row['peg No. (within edit)'])]
                else:
                    no_match[row['Edit Position (sat. area)']].append((row['Edit'], row['peg No. (within edit)']))
    
    if not no_match:
        print('All spacers from PRIDICT2.0 output are in seq!')
    else:
        print(no_match)
        

if __name__=='__main__':
    npc_seq = 'TACAGCTGGGTCTGACCTCTGAGTCCAGGGTCAGGTGATTTTGCTTAGCCTCAAGTGCTCAGATTCTGCTGATATTTTGCAAGACCTGGACTCTCTTGACACCCAGGATTCTTTCCTCAGGGGACATGCTGCCTATAGTTCTGCAGTTAACATCCTCCTTGGCCATGGCACCAGGGTCGGAGCCACGTACTTCATGACCTACCACACCGTGCTGCAGACCTCTGCTGACTTTATTGACGCTCTGAAGAAAGCCCGACTTATAGCCAGTAATGTCACCGAAACCATGGGCATTAACGGCAGTGCCTACCGAGTATTTCCTTACAGGTAAAGCCTGCCCTTTTTCAATGGGGTTTACCCAGCAAAGGGCCTACACTGGGTGGGAGTGGGGAGGGTTCCCTTGGCAAGATGCTGATTTTCAGGTTGGGTTCTGGCCCCTGCTCCATT'
    npc_sseq = 'ACTTA'

    # run_pridict_library_synony(npc_seq, npc_sseq).to_csv('./saturation_library/npc_result.csv', index=False)
    run_pridict_library_synony(npc_seq, npc_sseq, frame=2, HA=False, splice=[]).to_csv('./saturation_library/synony_npc_result.csv', index=False)
    
    