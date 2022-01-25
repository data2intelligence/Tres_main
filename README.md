Example code of run Tres analysis and reproduce main results  

# Stage 1: download single-cell datasets  
Please CD into the src folder and run "./download.py" first to download all data files, including 31 single-cell training datasets and 10 T-cell transcriptomics validation datasets.  
  
  
So, you may run "./run.py inx 31" where inx is a number between 0 and 31 to compute the Tres results for each dataset.    
We used the NIH high performance cluster (HPC) for parallel computation as "./hpc_submit.py 31". You may re-write this file for your local HPC.

The output will be available in data/output, named with the cancer type followed with database accession ID.  

# Stage 2: merge computation scores into signature files  
After computing the correlation results for each dataset, you can run "./run.py" to generate the merged Tres signature files in data/output:  
1, merge.TGFB1, merge.TRAIL, merge.PGE2: Tres scores among all samples for each dataset.  
2, merge.Median: Median signatures among all immunosuppressive signals enumerated above.  
3, merge.Median.AUC: Area under the ROC curve to measure the quality signatures in "merge.Median", using T-cell persistance markers (Krishna et al., Science 2021) as the evaluation standard.  
4, merge.signature: One overall Tres signature, which is the median across all datasets.  

# Stage 3: predict clinical response of immune checkpoint blockade, CAR T, and adoptive cell transfer  
After creating the Tres signature, you can run "./predict.py" to evaluate qualities of therapy response predictions (Figure 2 in the manuscript). All results files are named as *.pdf in data/output.  

**Task 1**: predict responder versus non-responder status using pretreatment samples:  
This task will analyze cohorts:  
Atezolizumab+Paclitaxel\_Pre\_Zhang2021\_TNBC  
CD19CAR\_Infusion\_Fraietta2018\_CLL  
ICB\_Post\_SadeFeldman2018\_Melanoma  
ICB\_Pre\_SadeFeldman2018\_Melanoma  
Nivolumab\_Post\_Caushi2021\_NSCLC  
anti-PD1\_Pre\_Yost2019\_BCC  
The clinical outcomes are responders or non-responders. Thus, we can generate Receiver Operating Characteristic (ROC) curves to compare signatures on predicting therapy response.  
The result for each dataset will contains one boxplot (*.boxplot.pdf) and one ROC curve (*.ROC.pdf).  

**Task 2**: predict cell therapy survival outcome using data from pre-manufacture samples:  
This task will analyze cohorts:  
CD19CAR\_Pre-manufacture_Chen2021\_Bcell  
ACT\_Pre-expansion\_Lauss2017\_Melanoma  
The clinical outcomes are survival durations and the sample expression profiles were taken before manufacturing the therapeutic T cells.  
Thus, we can analyze whether correlation with the overall Tres signature can predict the survival outcome through Kaplan-Meier plots and Cox-PH regression.  
The result for each dataset will contains one Kaplan-Meier plot (*.kmplot.pdf) and one barplot (*.bar.pdf) presenting Wald test risk z-scores.  

**Task 3**: differentiate mild versus severe symptoms of COVID19 using diagnosis samples:  
This task will analyze cohorts:  
COVID19 severity\_Diagnosis\_SchulteSchrepping2020\_PBMC  
COVID19 severity\_Diagnosis\_Su2020\_PBMC  
The clinical outcomes are mild or severe patients. Thus, we can generate Receiver Operating Characteristic (ROC) curves to compare signatures on differentiate infection outcomes.  
The result for each dataset will contains one boxplot (*.boxplot.pdf) and one ROC curve (*.ROC.pdf).  
