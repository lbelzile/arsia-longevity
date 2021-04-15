# arsia-longevity

This folder contains the code used to produce the analysis and figures in 

> Belzile, L.R., Davison, A.C., Gampe, J. Rootzén, H and D. Zholud (2022). Is there a cap on longevity?, submitted to ARSIA

The **R** package [`longevity`](https://github.com/lbelzile/longevity) can be installed from Github. It contains two datasets, `dutch` and `japanese`, previously published in other papers. 

The current version of the International Database on Longevity (www.supercentenarian.org) includes a large share of the records used in the paper, but these data are continuously updated with little documentation to track these changes (including the sampling frame provided in metadata files). Other data sources (including data formerly on IDL, and data from IStat and ONS), cannot be shared. As such, we cannot guarantee reproducibility of the analysis.

Recent versions of the IDL include both semi-supercentenarians (i.e., people who died age [105, 110)) who died between
calendar dates `d1` and `d2` and supercentenarians (i.e., who died aged 110+) who died between calendar dates `c1` and `c2` . The exact lifetime (in days) varies among individuals due to leap years. Jdanov et al. (2021) discuss the varying sampling frame, but the IDL data available from supercentenarian.org at the time of writing only includes records of dead individuals. Thus, we treated all IDL observations as interval truncated. 

For the purpose of replicability and assuming that the format of the IDL dataset does not change, `idl_preprocessing.R` shows how to transform the IDL database `idl_complete.csv` into an **R** dataset and calculate the truncation bounds for the (potentially) doubly interval truncated records.

---

Versions of the following scripts were used to obtain the figures and tables in the main text:

- Netherlands: `4.1-Dutch_analysis.R`
- Japanese: `4.2-Japanese_analysis.R`
- England and Wales: `4.3-EnglandWales_analysis.R`
- Profile log likelihoods for the endpoint by country: `6-profile_endpoint.R`

The tables in Section 5 of the Supplementary Material (Appendix F for the arXived version) were created using the Matlab toolbox **LATool** using slighly different data.

The script `SM-5-Comparisons.R` illustrates how similar calculations could be performed using the **R** package `longevity`.

