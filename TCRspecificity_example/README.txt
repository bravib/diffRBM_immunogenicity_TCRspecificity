The scripts are accompanied by the folder Example_diffRBM, where data and outputs are saved. Example_diffRBM contains training_data and test_data folders, 
where training data and test data should be placed. The files should contain the name of the peptide (for training peptide-specific models) and the peptide will be given as an input argument. See the folders for examples of the format of training and test data.  

STEP1: data preprocessing

Run the script:
python tcr_preprocessing.py -pep 'GLCTLVAML'

for the pre-processing (i.e., alignment) of the training and test data for a given peptide (here GLCTLVAML). The routines for the alignment are in the Align_Utils folder - the path to the right folder must be set in these scripts.  Requirements: the alignment routines require matlab and matlab engine API for Python https://fr.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html


STEP2: train and evaluate the models

Download PGM and diffRBM package from github:
https://github.com/jertubiana/PGM
https://github.com/cossio/diffRBM

and place them in the same folder as Example_diffRBM. 

To train and/or evaluate the models on the test data, run: 

python diffrbm_tcr_run.py -pep 'GLCTLVAML' -train 1

Arguments:
For tcr_preprocessing.py and diffrbm_tcr_run.py:
-pep: name of the peptide target for training a peptide-specific model (Required)
-alpha: Choose CDR3 on alpha chain instead of beta chain (default: False)

For diffrbm_tcr_run.py:
-train: Whether to train a model with the peptide-specific training data (default: False)
-LR: Choose the Left-Right format instead of the alignment (default: False)
-hu: Number of Top RBM hidden units (default: 20)
-l12: L^1_2 Top RBM regularization (default: 0.01)
-verbose: Switch off verbose mode

The output is in the folder diffRBM_GLCTLVAML and consists of the test data with the assigned scores (see 'test_data_scored ...' file)
