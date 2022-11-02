# diffRBM_immunogenicity_TCRspecificity
This repo contains data, codes and other results related to the manuscript 'Learning the differences: a machine learning approach to
predicting antigen immunogenicity and T-cell receptor specificity'.

'Immunogenicity_model' contains:
- training/test data and models for HLA-A\*02:01-presented peptides (same for HLA-B\*35:01, HLA-B\*07:02);
- Tables with contact positions on PDB structures and models' predictions;
- Results on the immunogenic vs non-immunogenic discrimination performance.

NB: To re-train the immunogenicity models, the file tcell_full_v3.zip is needed (from here: https://www.iedb.org/database_export_v3.php)

'TCR_specificity_model' contains:
- training/test data and models for NLVPMVATV-specific CDR3beta (same for GILGFVFTL, GLCTLVAML, YLQPRTFLL);
- Tables with contact positions on PDB structures and models' predictions;
- Results on the antigen-specific vs generic discrimination performance (for diffRBM and other methods).

'Notebooks' contains:
- a notebook to produce the results for the diffRBM immunogenicity model and the immunogenicity classifier (notebook_immunogenicity.ipynb);
- a notebook to produce the results for the diffRBM T-Cell Receptor specificity model (notebook_TCR_specificity.ipynb);
- a notebook to produce the figures of the manuscript.

'Align_utils' contains routines for sequence alignment.

'TCR_specificity_example' contains codes that illustrate how to train and evaluate the diffRBM approach to modeling TCR specificity.

'Immunogenicity_example' contains codes that illustrate the diffRBM approach to modeling immunogenicity.
