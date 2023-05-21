# Conformalized Matrix Prediction

Each simulation/data example has its own main file.  See details below.

<h2>Reference:</h2>

* Distribution-Free Matrix Prediction Under Arbitrary Missing Pattern <br />
<i>Meijia Shao, Yuan Zhang</i><br>

BibTeX
```bibtex
@article{shao2023distribution,
    AUTHOR = {Shao, Meijia and Zhang, Yuan},
     TITLE = {Distribution-Free Matrix Prediction Under Arbitrary Missing Pattern},
   JOURNAL = {ArXiv preprint arXiv:},
      YEAR = {2023},
}


```

<h2>Contents:</h2>
This README file contains step-by-step instructions on how to reproduce the simulation and data example results in this paper.


<h2>Remarks:</h2>
<ul>
  <li> The original coding files for implementing <i>Sportisse et al. (2020)</i> in 'simu1_PPCA.R' are not provided in the following coding, please check on the following github link to get the code: 
https://github.com/AudeSportisse/PPCA_MNAR 
  <li> Data example: if you would like to reproduce the results for data example 2, please contact the owner of that data set (see "Data_permission.txt" for more details) before running the code.
</ul>


<h2>Subroutine list</h2>

1. 'transforming_A.m'
2. 'NeighborhoodSmoothing.m'
3. 'kernel_Y.m'
4. 'graphon.m'  
5. 'generate_randW.m'
6. 'data_generate.R'
7. 'data_generate_error.R'
8. 'shadedErrorBar.m'  (publically available routine, https://github.com/raacampbell/shadedErrorBar)



<h2>Instructions for result reproduction</h2>

<h3>Preparation:</h3>

Before running our code, please:

1. Place all coding files in main working folder, and create three subfolders named 'result','figure' and 'data';
2. Place the data files in a subfolder named 'data'.

<h3>To reproduce simulation 1 results in Section 4:</h3>

1. Run 'simu1_our_algorithm1.m';
2. Run 'simu1_our_algorithm2.m';
3. Run 'simu1_soft.R';
4. Run 'simu1_missMDA.R';
5. Run 'simu1_mice.R';
6. Run 'simu1_PPCA.R';
7. Run 'plot_simu1_main.m' to reproduce Figure 1;
8. Run 'plot_boxplot_simulation1.m' to reproduce Figure 3.  


<h3>To reproduce simulation 2 results in Section 4:</h3>

1. Run 'simu2_our_algorithm1.m';
2. Run 'simu2_our_algorithm2.m';
3. Run 'simu2_soft.R';
4. Run 'simu2_missMDA.R';
5. Run 'simu2_mice.R';
6. Run 'simu2_PPCA.R';
7. Run 'plot_simu2_main.m' to reproduce Figure 2;
8. Run 'plot_boxplot_simulation2.m' to reproduce Figure 3. 

<h3>To reproduce data results in Section 4:</h3>

1. Run 'data_our_algorithm1.m';
2. Run 'data_softImpute.R';
3. Run 'plot_data_result.R' to reproduce Figure 4. 

<h3>To reproduce Figure 1 \& 2 in Supplementary Materials:</h3>

1. Run 'plot_simu1_supp.m' to reproduce Figure 1. 
2. Run 'plot_simu2_supp.m' to reproduce Figure 2. 

<h3>Computing resources:</h3>

1. Simulations and data example were run on (anonymized university's) computing server, 10 parallel Intel(R) Xeon(R) CPU's, (model specification anonymized), 1GB requested memory for each task, MATLAB R2022a.
2. Simulation and data example time cost: reported as part of results.



