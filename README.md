


# Instructions for Running R Scripts
- The R script `cross_fitting.R` implements the functions used in the algorithms described in the paper.
- The directory `Script` contains the R scripts used to reproduce the results reported in the paper.
- The directory `data` stores the generated output files, including `.Rdata`, `.pdf`, and `.csv` files.

Before running the R scripts that generate the results reported in the paper, please follow the steps below.  All commands should be executed in the directory `Simulation`.

## 1. Monte Carlo Simulations
For \( n = 5000 \), run the following R scripts in the terminal before you run the Rscripts for Section 6:
```bash
nohup Rscript Script/MonteCarlo1.R > logs/MonteCarlo1.log 2>&1 &
nohup Rscript Script/MonteCarlo2.R > logs/MonteCarlo2.log 2>&1 &
nohup Rscript Script/MonteCarlo3.R > logs/MonteCarlo3.log 2>&1 &
nohup Rscript Script/MonteCarlo4.R > logs/MonteCarlo4.log 2>&1 &
nohup Rscript Script/MonteCarlo5.R > logs/MonteCarlo5.log 2>&1 &
nohup Rscript Script/MonteCarlo6.R > logs/MonteCarlo6.log 2>&1 &
```
For \( n = 5000 \), run the following R scripts in the terminal before you run the Rscripts for Section 6:
```bash
nohup Rscript Script/MonteCarlo1_10000.R > logs/MonteCarlo1_10000.log 2>&1 &
nohup Rscript Script/MonteCarlo2_10000.R > logs/MonteCarlo2_10000.log 2>&1 &
nohup Rscript Script/MonteCarlo3_10000.R > logs/MonteCarlo3_10000.log 2>&1 &
nohup Rscript Script/MonteCarlo4_10000.R > logs/MonteCarlo4_10000.log 2>&1 &
nohup Rscript Script/MonteCarlo5_10000.R > logs/MonteCarlo5_10000.log 2>&1 &
nohup Rscript Script/MonteCarlo6_10000.R > logs/MonteCarlo6_10000.log 2>&1 &
```
Run the following R scripts in the terminal before you run the Rscripts for Section 7:
```bash
nohup Rscript Script/JTPA_9.R > logs/JTPA_9.log 2>&1 &
nohup Rscript Script/JTPA_11.R > logs/JTPA_11.log 2>&1 &
nohup Rscript Script/JTPA_12.R > logs/JTPA_12.log 2>&1 &
nohup Rscript Script/JTPA_9_deg.R > logs/JTPA_9_deg.log 2>&1 &
nohup Rscript Script/JTPA_11_deg.R > logs/JTPA_11_deg.log 2>&1 &
nohup Rscript Script/JTPA_12_deg.R > logs/JTPA_12_deg.log 2>&1 &
```

## 2. Figures
Run the following R scripts in the terminal to generate figures in Section 2:

```bash
Rscript Script/Figure2.1\(a\).R
Rscript Script/Figure2.1\(b\).R
Rscript Script/Figure2.2\(a\).R
Rscript Script/Figure2.2\(b\).R
```

Run the following R scripts in the terminal to generate figures for simulations:
```bash
nohup Rscript Script/Figure6.1\(a\).R
nohup Rscript Script/Figure6.1\(c\).R
nohup Rscript Script/Figure6.1\(e\).R
nohup Rscript Script/Figure6.1\(b\).R
nohup Rscript Script/Figure6.1\(d\).R
nohup Rscript Script/Figure6.1\(f\).R

nohup Rscript Script/FigureC.1\(a\).R
nohup Rscript Script/FigureC.1\(b\).R
nohup Rscript Script/FigureC.1\(c\).R
```

Run the following R scripts in the terminal to generate figures in real data analysis:
```bash
nohup Rscript Script/Figure7.1\(a\)\(b\).R
nohup Rscript Script/Figure7.1\(c\).R
```
## 3. Tables

To generate the tables in Appendix A, run the following R scripts in the terminal:

```bash
nohup Rscript Script/TabelC.1\(n\=5000\).R
nohup Rscript Script/TabelC.2\(n\=5000\).R
nohup Rscript Script/TabelC.3\(n\=5000\).R
nohup Rscript Script/TabelC.1\(n\=10000\).R
nohup Rscript Script/TabelC.2\(n\=10000\).R
nohup Rscript Script/TabelC.3\(n\=10000\).R

nohup Rscript Script/TabelC.4.R
```

