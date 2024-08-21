import subprocess
import os
import time
import pandas as pd

U6_homology_arm = "ttatatatcttgtggaaaggacgaaacacc"
pe_scaffold = "gtttcagagctatgctggaaacagcatagcaagttgaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc"
tevopreQ1 = "CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA"
tevopreQ1_ori = 'cgcggttctatctagttacgcgttaaaccaactagaatttttttaagcttgggccgctcgaggtacctctctacatatgacatgtgagcaaaaggccagcaaaaggccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggt'


def pred(wts='', batch=False):
    if batch is True:
        command = [
            'python', '/Users/Dong-Kyu Kim/PRIDICT2_library/pridict2_pegRNA_design.py', 'batch',  # CHANGE DIRECTORY NAME HERE
            '--input-fname', 'batch_template.csv',
            '--output-fname', 'batchseqs'
        ]

        csv_file_path = 'batchseqs'

    else:
        if len(wts) < 203:
            return 'not available - wts+edit too short'

        command = [
            'python', '/Users/Dong-Kyu Kim/PRIDICT2_library/pridict2_pegRNA_design.py', 'manual',  # CHANGE DIRECTORY NAME HERE
            '--sequence-name', 'seq',
            '--sequence', wts
        ]

        csv_file_path = '/Users/Dong-Kyu Kim/PRIDICT2_library/predictions/seq_pegRNA_Pridict_full.csv'  # CHANGE DIRECTORY NAME HERE

    subprocess.run(command, check=True)

    while not os.path.exists(csv_file_path):
        time.sleep(1)
    
    return None



def manual_pred(wts, rtt, pbs):

    wts, rtt, pbs = wts.upper(), rtt.upper(), pbs.upper()

    if len(wts) < 203:
        return 'not available - wts+edit too short'

    command = [
        'python', '/home/yjsk/mysite/PRIDICT2/pridict2_pegRNA_design.py', 'manual',
        '--sequence-name', 'seq',
        '--sequence', wts
    ]

    csv_file_path = '/home/yjsk/mysite/PRIDICT2/predictions/seq_pegRNA_Pridict_full.csv'

    try:
        subprocess.run(command, check=True)
    except Exception as e:
        print(f"Error processing: {e}")
        return 'not available'

    while not os.path.exists(csv_file_path):
        time.sleep(1)

    pridict2_output = pd.read_csv(csv_file_path)

    os.remove(csv_file_path)

    pridict2_output['RTrevcomp'] = pridict2_output['RTrevcomp'].str.upper()
    pridict2_output['PBSrevcomp'] = pridict2_output['PBSrevcomp'].str.upper()

    score_row = pridict2_output.loc[
        (pridict2_output['RTrevcomp'] == rtt) &
        (pridict2_output['PBSrevcomp'] == pbs),
        'PRIDICT2_0_editing_Score_deep_HEK'
    ]

    if score_row.empty:
        return 'RTT not in model output'
    else:
        return score_row.iloc[0]


