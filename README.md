Example code of run Tres analysis and reproduce main results  

# Stage 1: download single-cell datasets  
Please CD into the src folder and run "./download.py" first to download all data files, including 16 single-cell training datasets and 5 validation datasets.  
  
So, you may run "./run.py inx 16" where inx is a number between 0 and 16 to compute the Tres results for each dataset.  
We used the NIH high performance cluster (HPC) for parallel computation as "./hpc_submit.py 16". You may re-write this file for your local HPC.  

The output will be available in data/output, named with the cancer type followed with database accession ID.  

# Stage 2: merge computation scores into signature files  
After computing the correlation results for each dataset, you can run "./run.py" to generate the merged Tres signature files in data/output:  
1, merge.TGFB1, merge.TRAIL, merge.PGE2: Median signatures among all samples for each dataset.  
2, merge.Median: Median signatures among all immunosuppressive signals enumerated above.  
3, merge.Median.AUC: Area under the ROC curve to measure the quality signatures in "merge.Median", using T-cell persistance markers (Krishna et al., Science 2021) as the evaluation standard.  
4, merge.Median.signature: One overall Tres signature, which is the median among datasets with AUC > 0.7 in "merge.Median".  

# Stage 3: predict clinical response of immune checkpoint blockade, CAR T, and adoptive cell transfer  
After creating the Tres signature, you can run "./predict.py" to evaluate qualities of therapy response predictions (Figure 2 in the manuscript). All results files are named as Validation.*.pdf in data/output.  

**Task 1**: predict responder versus non-responder status:  
This task will analyze cohorts: ICB\_Caushi2021, ICB\_SadeFeldman2018.Exhausted.Post, ICB\_SadeFeldman2018.Exhausted.Pre, CD19CAR_Fraietta2018.  
The clinical outcomes are responders or non-responders. Thus, we can generate Receiver Operating Characteristic (ROC) curves to compare signatures on predicting therapy response.  
The result for each dataset will contains one boxplot and one ROC curve.  

**Task 2**: predict cell therapy survival outcome using data from pre-manufacture samples:  
This task will analyze cohorts: ACT\_Lauss2017 and CD19CAR\_Chen2021.  
The clinical outcomes are survival durations and the sample expressoin profiles were taken before manufacturing the therapeutic T cells.  
Thus, we can analyze whether correlation with the overall Tres signature can predict the survival outcome through Kaplan-Meier plots and Cox-PH regression.  
