#!/bin/sh
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(1,1)_p50.RData' r=2 n=512 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2 " Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(0,5)_p50.RData' r=2 n=512 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(1,1)_p50.RData' r=2 n=2^11 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2 " Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(0,5)_p50.RData' r=2 n=2^11 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(1,1)_p50.RData' r=2 n=190 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2 " Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(0,5)_p50.RData' r=2 n=190 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(1,1)_p50.RData' r=2 n=2^12 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2 " Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(0,5)_p50.RData' r=2 n=2^12 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt


R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(1,1)_p50.RData' r=2 n=512 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2 " Simulation_Part_Coherence_Test_FDR_Lambda_Seq_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(0,5)_p50.RData' r=2 n=512 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Lambda_Seq_Repro.r Output.txt



R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(1,1)_p100.RData' r=2 n=190 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(0,5)_p100.RData' r=2 n=190 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(1,1)_p100.RData' r=2 n=2^9 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(0,5)_p100.RData' r=2 n=2^9 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(1,1)_p100.RData' r=2 n=2^11 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(0,5)_p100.RData' r=2 n=2^11 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(1,1)_p100.RData' r=2 n=2^12 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(0,5)_p100.RData' r=2 n=2^12 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt

R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(1,1)_p200.RData' r=2 n=190 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(0,5)_p200.RData' r=2 n=190 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(1,1)_p200.RData' r=2 n=2^9 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(0,5)_p200.RData' r=2 n=2^9 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(1,1)_p200.RData' r=2 n=2^11 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(0,5)_p200.RData' r=2 n=2^11 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(1,1)_p200.RData' r=2 n=2^12 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt
R CMD BATCH --no-save --no-restore "--args file='DGP_VARMA(0,5)_p200.RData' r=2 n=2^12 VAR_order=round(log10(n)) Kernel=Modified_Bartlett B=1000 M1=NULL Delta1=0 Delta2=0.2" Simulation_Part_Coherence_Test_FDR_Repro.r Output.txt



