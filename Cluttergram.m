clear
clc
close all


% Resolution in number of data points e.g. resolution = 2, we will use half
% of the orbiter positions.
resolution=0.5;
% Resolution of the MOLA topography in steps
mola_res=20;
% Open the file with the position of the orbiter
[filee,path] = uigetfile('*.txt');
m = importdata([path,filee]);
% Open the file with Mars topography
M = csvread( 'MOLA_Jezero_rec.csv' );




[f_x,ti,RR] = Cluttergram_mat(resolution,mola_res,m,M);







