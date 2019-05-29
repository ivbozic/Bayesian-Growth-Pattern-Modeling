# Bayesian-Growth-Pattern-Modeling

This code performs Bayesian inference of parameters of a logistic growth model using Gibbs sampling.

The main script is log_Bayesian_inference.m. It requires two input files:
* WBC_data.xlsx, with columns PatientID TimeSinceDx WBC
* Dx_to_Tx_data.xlsx, with columns PatientID TimeFromDxToTx

Inferred parameter values for each patient are saved in the file results.txt. Plots of the parameter posterior distributions for each patient are saved into posterior_histograms.pdf.
