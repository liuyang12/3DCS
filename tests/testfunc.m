%TESTFUNC test the functions of the project

%% [0] add *algorithms* and *packages* to operation path
addpath('../algorithms')
addpath(genpath('../packages'))

%% [3.1.1] test the TV operator generator
m = 3; % number of rows
n = 4; % number of columns
H = gentvop(m,n,1);
