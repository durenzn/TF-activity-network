# TF-activity-network
To study systems-level properties of the cell, it is necessary to go beyond individual regulators and target genes to study the regulatory network among transcription factors (TFs). However, it is difficult to directly measure TF activity or DNA-binding fraction of TF by experiment. Here, we proposed a hierarchical graphical model to estimate TF activity from mRNA expression by building TF complexes with protein co-factors and inferring TF’s downstream regulatory network simultaneously. 


This is the MARLAB source code for infering Tf activity network. 
The main program for TF activity network inference is in the file named TFActivyNetwork.m
The directory /demo gives an small example by running demo.m 
The input files in /demo directory can serve as the templates the parameters in demo.m can be the default ones.


# Inputs:

'LA' is the gene expression data, in which rows denote genes and columns denote samples;
'TFid' is the index of the TFs.
'ppi' is the protein-protein interaction network. The first column is the index of TF and the second column is the index of the gene.
'arfa' is the p-value cutoff of the Transcriptional Regulatory Network inferred in Step one.
'lambda' is the paramater balancing the sparsity and likelihood.

# Outputs:

'NW' is the matrix of Transcriptional Regulatory Network with m TFs and N target genes.
'ccmi' is the condotional mutual information.
'Modulator' is the index matrix of modulators.
'AA' is the effect of modulator on TF.
'TFA' is the activity of TF.
'BB' is the effect of TF on target gene.
'TFtN' is the information of TF-target network,it will be printed in the file named TF-target network.txt
'MTFtN' is the information of Modulator-TF-target network,it will be printed in the file named Modulator-TF-target network.txt
'MTFNet and TFtNet' is the network of Modulator-TF and TF-target.
