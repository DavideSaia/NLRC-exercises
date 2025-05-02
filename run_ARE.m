%Script for testing ARE function
clear;
clc;
close all;

A = [1 1 0; 0 1 1; 1 0 1];  
B = [1 1 0;3 1 0; 0 0 1];     
Q = eye(3);    
R = 1;  
n = size(A,1);
S = zeros(n,n);

%Calculate P using ARE function
P_ARE = ARE(A,B,Q,R,S);

