import random
import pandas as pd
from Bio.Seq import Seq
from model_pred import pred
from main import process_row, _get_control_rtt, SF, fivep_homo, threep_homo, _get_synony_rtt, _random_filler, get_preserving_rtt, _c

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

    filename = './saturation_library/inputs.txt'

    # with open(filename, 'w') as outfile:
    #     outfile.write('\n'.join(str(i) for i in allSeq))
 
    return allSeq
 


def run_pridict_lib(seq, sseq):

    def sort_result(seq, unsorted_rs):
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
    unsorted_result = pd.DataFrame()
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

            df['Complete pegRNA'] = [fivep_homo + top_spacers.iloc[j]['Spacer-Sequence'] + 'GTTTCGAGACG' + _random_filler() + 'CGTCTCGGTGC' + top_spacers.iloc[j]['RTrevcomp'] + top_spacers.iloc[j]['PBSrevcomp'] + threep_homo]
            df['Length (bp)'] = df['Complete pegRNA'].str.len()

            df['Complete pegRNA (SF)'] = [fivep_homo + top_spacers.iloc[j]['pegRNA'] + threep_homo] # Uses scaffold
            df['Length (bp) (SF)'] = df['Complete pegRNA (SF)'].str.len()

            df['Reference Sequence'] = [top_spacers.iloc[j]['wide_mutated_target']]
            df['PRIDICT2.0 Score'] = [top_spacers.iloc[j]['PRIDICT2_0_editing_Score_deep_HEK']]
            dfs.append(df)

    unsorted_result = pd.concat(dfs, ignore_index=True)
    sorted_result = sort_result(seq, unsorted_result)
    
    return sorted_result


def run_pridict_library_synony(seq, sseq, frame, HA, splice):
    # df = run_pridict_lib(seq, sseq)
    df = pd.read_csv('./saturation_library/npc_result.csv')

    rows = []
    for i, row in df.iterrows():
        row_ = row.copy()
        if row['Strand']=='(-)':
            wt_rtt = get_wt_rtt(seq, row['RTTs'])
        else:
            wt_rtt= get_wt_rtt(seq, row['RTT'])

        synony = _get_synony_rtt(seq=seq, sseq=sseq, rtt=wt_rtt, frame=frame, strand=row['Strand'][1], splice=splice)
        row_['RTTs'] = get_preserving_rtt(synony, row['RTTs'], wt_rtt)

        row['Complete pegRNA'] = [fivep_homo + row['Spacer'] + row['Filler'] + row['RTTs'] + row['PBS'] + threep_homo]
        row['Complete pegRNA (SF)'] = [fivep_homo + row['Spacer'] + SF + row['RTTs'] + row['PBS'] + threep_homo]
        # row['Reference Sequence'] = process_row(row, seq, wt_rtt)
        rows.append(row_)

    pd.DataFrame(rows).to_csv('./saturation_library/synony.csv')


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
    # TODO: How to get ref seq in synony lib, add ctls - sort by PAM (TGA + 1bp del per PAM)
    npc_seq = 'TACAGCTGGGTCTGACCTCTGAGTCCAGGGTCAGGTGATTTTGCTTAGCCTCAAGTGCTCAGATTCTGCTGATATTTTGCAAGACCTGGACTCTCTTGACACCCAGGATTCTTTCCTCAGGGGACATGCTGCCTATAGTTCTGCAGTTAACATCCTCCTTGGCCATGGCACCAGGGTCGGAGCCACGTACTTCATGACCTACCACACCGTGCTGCAGACCTCTGCTGACTTTATTGACGCTCTGAAGAAAGCCCGACTTATAGCCAGTAATGTCACCGAAACCATGGGCATTAACGGCAGTGCCTACCGAGTATTTCCTTACAGGTAAAGCCTGCCCTTTTTCAATGGGGTTTACCCAGCAAAGGGCCTACACTGGGTGGGAGTGGGGAGGGTTCCCTTGGCAAGATGCTGATTTTCAGGTTGGGTTCTGGCCCCTGCTCCATT'
    npc_sseq = 'ACTTA'

    run_pridict_lib(npc_seq, npc_sseq).to_csv('./saturation_library/npc_result.csv', index=False)
    verify_spacer_in_seq(npc_seq)
    # run_pridict_library_synony(npc_seq, npc_sseq, frame=2, HA=True, splice=[])
    # print(get_wt_rtt(npc_seq, 'GAAGAAAGCCCGtCTTATAGCCAGTAATGTCACCGAAA'))
    