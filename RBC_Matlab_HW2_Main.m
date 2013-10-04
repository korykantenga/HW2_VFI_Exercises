cd('g:\tsprofiles\kantenga\Desktop\Dropbox\ECON714\HW2')

diary('RBC_Model_Main_Ouput')

tic;

%{ Notes: Running only 1, 2, 3 and 5 takes about an hour. %}

%1. Comparison of Initial Guess

RBC_Matlab_HW2_1_Guesses

%2. Comparison of Grids

RBC_Matlab_HW2_2_Grids

%3. Accelerator

RBC_Matlab_HW2_3_Accelerator

%4. Multigrid Algoritm

RBC_Matlab_HW2_4_Multigrid

%5. Stochastic Grid

RBC_Matlab_HW2_5_StochasticGrid

%6. Endogenous Grid



toc;

diary off