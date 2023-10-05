# PRIDICT2.0: PRIme editing guide RNA preDICTion 2.0

![PRIDICT logo](pridict_logo.png)

This is the repository to run the `PRIDICT2.0` model on your local machine. `PRIDICT2.0` builds on previous work from [PRIDICT](https://rdcu.be/c3IM5) and further details on changes/improvements can be accessed [here](xxx).
Also check out our `ePRIDICT` prediction model ([link to repository](https://github.com/Schwank-Lab/epridict)), which is complementary to the sequence-context based `PRIDICT2.0` as it predicts how local chromatin-context influences prime editing efficiencies.

For accessing Supplementary Files, click [here](xxx).

To run `PRIDICT2.0` online, see our [webapp](https://pridict.it/).

--------------------------

### Installation using Anaconda (Linux, Mac OS or Windows) 🐍

The easiest way to install and manage Python packages on various OS platforms is through [Anaconda](https://docs.anaconda.com/anaconda/install/). Once installed, any package (even if not available on Anaconda channel) could be installed using pip. 

* Install [Anaconda](https://docs.anaconda.com/anaconda/install/).
* Start a terminal and run:
    ```shell
    # clone PRIDICT2.0 repository
    git clone https://github.com/uzh-dqbm-cmi/PRIDICT2.git
    # navigate into repository
    cd PRIDICT2
    # create conda environment and install dependencies for PRIDICT2 (only has to be done before first run/install)
    conda env create -f pridict2_repo.yml
        
    # activate the created environment
    conda activate pridict2
    
    	### ONLY FOR M1 (or newer) Mac you need to additionally run the following conda install command (tensorflow): 
    	conda install conda-forge::tensorflow
    	# optional (only if encountering error with libiomp5.dylib on MacOS):
    	pip uninstall numpy
    	pip install numpy==1.22.1
    	###
    	
	
    # run desired PRIDICT2.0 command (manual or batch mode, described below)
    python pridict2_pegRNA_design.py manual --sequence-name seq1 --sequence 'GCCTGGAGGTGTCTGGGTCCCTCCCCCACCCGACTACTTCACTCTCTGTCCTCTCTGCCCAGGAGCCCAGGATGTGCGAGTTCAAGTGGCTACGGCCGA(G/C)GTGCGAGGCCAGCTCGGGGGCACCGTGGAGCTGCCGTGCCACCTGCTGCCACCTGTTCCTGGACTGTACATCTCCCTGGTGACCTGGCAGCGCCCAGATGCACCTGCGAACCACCAGAATGTGGCCGC'
    # results are stored in 'predictions' folder
    ```

* `PRIDICT2.0` environment only has to be installed once. When already installed, follow the following commands to use `PRIDICT2.0` again:
    ```shell
    # open Terminal/Command Line
    # navigate into repository
    # activate the created environment
    conda activate pridict2
    # run desired PRIDICT2.0 command (manual or batch mode, described below)
    python pridict2_pegRNA_design.py manual --sequence-name seq1 --sequence 'GCCTGGAGGTGTCTGGGTCCCTCCCCCACCCGACTACTTCACTCTCTGTCCTCTCTGCCCAGGAGCCCAGGATGTGCGAGTTCAAGTGGCTACGGCCGA(G/C)GTGCGAGGCCAGCTCGGGGGCACCGTGGAGCTGCCGTGCCACCTGCTGCCACCTGTTCCTGGACTGTACATCTCCCTGGTGACCTGGCAGCGCCCAGATGCACCTGCGAACCACCAGAATGTGGCCGC'
    # results are stored in 'predictions' folder
    ```

--------------------------

### Running PRIDICT2.0 in 'manual' mode:
  ####  Required:
  -  `--sequence-name`: name of the sequene (i.e. unique id for the sequence)
  -  `--sequence`: target sequence to edit in quotes (format: `"xxxxxxxxx(a/g)xxxxxxxxxx"`; minimum of 100 bases up and downstream of brackets are needed; put unchanged edit-flanking bases *outside* of brackets (e.g. xxxT(a/g)Cxxx instead of xxx(TAC/TGC)xxx)
  ####  Optional:
  -  `--output-dir`: output directory where results are dumped on disk (default: `./predictions`; directory must already exist before running)
  -  `--use-5folds`: Use all 5-folds trained models. Default is to use fold-1 model
  -  `--cores`: Number of cores to use for multiprocessing. Maximum 3 cores due to memory limitations. Default value 0 uses 3 cores if available.
  -  `--nicking`: Additionally, design nicking guides for edit (PE3) with DeepSpCas9 prediction.
  -  `--ngsprimer`: Additionally, design NGS primers for edit based on Primer3 design.
```shell

python pridict2_pegRNA_design.py manual --sequence-name seq1 --sequence 'GCCTGGAGGTGTCTGGGTCCCTCCCCCACCCGACTACTTCACTCTCTGTCCTCTCTGCCCAGGAGCCCAGGATGTGCGAGTTCAAGTGGCTACGGCCGA(G/C)GTGCGAGGCCAGCTCGGGGGCACCGTGGAGCTGCCGTGCCACCTGCTGCCACCTGTTCCTGGACTGTACATCTCCCTGGTGACCTGGCAGCGCCCAGATGCACCTGCGAACCACCAGAATGTGGCCGC'
``` 
--------------------------

### Running in batch mode:
  ####  Required:
  -  `--input-fname`: input file name - name of csv file that has two columns [`editseq`, `sequence_name`]. See `batch_template.csv` in the `./input` folder
  ####  Optional:
  -  `--input-dir` : directory where the input csv file is found on disk
  -  `--output-dir`: directory on disk where to dump results (default: `./predictions`)
  -  `--output-fname`: output filename used for the saved results
  -  `--combine-results`: Compile all results in one dataframe
  -  `--use-5folds`: Use all 5-folds trained models. Default is to use fold-1 model
  -  `--cores`: Number of cores to use for multiprocessing. Maximum 3 cores due to memory limitations. Default value 0 uses 3 cores if available.
  -  `--nicking`: Additionally, design nicking guides for edit (PE3) with DeepSpCas9 prediction.
  -  `--ngsprimer`: Additionally, design NGS primers for edit based on Primer3 design.
```shell

 python pridict2_pegRNA_design.py batch --input-fname batch_template.csv --output-fname batchseqs

``` 
--------------------------
