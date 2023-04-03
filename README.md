# TNBC-ICI

TNBC-ICI is a novel, random forest-based classifier than can efficiently predict response to immunotherapy in early triple-negative breast cancer patients. Here you can find all the employed code to generate this classifier. You can also test your own data to check the efficiency of TNBC-ICI in new cohorts. **Do not use this classifier for medical reasons and to change the current treatment of any patient.**

> Triple-Negative Breast Cancer (TNBC) is the breast cancer subtype with the highest mortality. In the last years, immunotherapy has arisen as a promising treatment in various cancers, including TNBC, but a subset of TNBC patients still does not respond to this treatment. Here, we employed machine learning algorithms to create a classifier that predicts response to immunotherapy plus chemotherapy in primary TNBC patients. This classifier has higher accuracy than other tested informative biomarkers, displaying a lower efficacy in non-TNBC breast cancer patients and other non-breast cancers. This classifier may be used to better select patients for immunotherapy, upfront, decreasing the side effects and the costs of the treatment.

![image](https://user-images.githubusercontent.com/46361666/213405474-2f353c64-4adb-4e97-aaf1-5832f6e5ae2a.png)

## Replicate the classifier

To validate the code employed to generate the classifier, enter the folder "Validate TNBC-ICI code" and use the available script, loading the uploaded environment. Please, be aware that this classifier has been generated using random forest. As its name says, this method has a certain inherent _randomness_. Thus, the obtained classifier from your iterations might be entirely different from the 37-gene one published in the paper. An additional step has been added to the script to directly employ the 37-gene signature to perform the analysis and replicate the results from the publication (with potential minor changes due to the random effect of the algorithm).

NOTE: The replication of the entire analysis can be extremely long due to the number of iterations of the random forest algorithm (>6h). If you want to do a shorter analysis, reduce the number of iterations. Also, if you want to replicate the results of the manuscript using the 37-gene signature, the time will be extremely reduced to under 5-10 minutes.

## Test your samples

If you want to test the efficiency of TNBC-ICI in your immunotherapy-treated cohort, enter the folder "Test new samples". In this folder you can find the necessary script and environment to train the classifier employing the available genes in your cohort. In case that your cohort has missing data or does not have some of the genes in the 37-gene signature, the script is prepared to train the algorithm using the available genes. **However, we do not guarantee that the efficiency of the classifier will be maintained using other combinations of genes.**

## Software repository

All analysis have been performed in Windows 11 using R v.4.0.2 in RStudio. The following R packages have been used:
- pROC_1.16.2
- tidyverse_1.3.0
- varSelRF_0.7-8
- randomForest_4.6-14
- lamisc_0.0.0.9000
- splitTools_0.3.2
- ggplot2_3.3.6
- M3C_1.10.0
- sva_3.38.0
- verification_1.42

Instructions for the installation of all packages are included in the "Validate TNBC-ICI code/Validate TNBC-ICI Simplified.R" script.

## Contact
If there is any issue with the code, contact Miquel Ensenyat: m.ensenyat.mendez@gmail.com

Furthermore, if you use this classifier, please contact us so we can improve the accuracy of the classifier with new data. Be a part of it!


