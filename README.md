# HGSOC_ITH
Software to reproduce analyses carried out for research article: 'Spatial and temporal intra-tumoural heterogeneity in advanced High-Grade Serous Ovarian Cancer: implications on surgical and clinical outcomes'

Prerequisites:
- R version >= 3.4.1
- R packages GenomicRanges, ggplot2, reshape2, NMF, YAPSA, flexmix, ASCAT
- CN signature helper functions in files 'helper_functions.R' and 'main_functions.R' (from https://github.com/markowetzlab/CNsignatures)
- text file containing clinical information for patients from study 'ITHclinicaldata.csv'
- text files containing mutation data for samples 'somatic_HighModImact_variants.csv' and 'germline_highimpact_variants.csv'

Steps:
1. Run ASCAT on raw data files to create ASCAT object (described at https://github.com/VanLoo-lab/ascat)

2. Load required functions contained in files 'HGSOC_ITH_functions.R'

3. Load clinical data by running 'ITH_clinicalprocess.R'

4. Annotated ASCAT-processed CN profiles for all samples using 'AnnotateASCATobject.R'

5. Calculate number-of-event distances between each pair of tumours, for each patient using 'getNEventDists.R'                 

6. Estimate number of clones detected in the tumours of each patient using 'ITH_nClones.R' 

7. Calculate copy-number signatures per sample using 'generateSignaturesITH.R'

8. Visualize copy-number signatures per patient, and perform analyses required to create Fig 3 using 'analyseSignaturesITH.R'

9. Evaluate probability of observing mixed HRD scores for a patient when measuring different deposits, based on measured score in a single deposit, using 'HRDanalysis.R'

10. Map somatic mutations to patient samples using 'mutationAnalysis.R'
