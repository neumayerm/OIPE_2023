% Code for the OIPE Paper
% 
% Functions
%
% EMS 2024
% Contact: neumayer@tugraz.at
clear all, close all, clc

addpath .\lib

load SENSOR_FEM_C.mat


% Plot of the material distributions
%func_FEMplot(FEM)

% Solver to simulate the sensor
[FEM] = func_solve_C(FEM);

% Comutation of complex valued Jacobian (reference method)
 Javm = func_getJ_AVM_C(FEM);
 
 % Comutation of complex valued Jacobian using Green's functions
 Jgreen = func_getJ_Green_C(FEM);

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of J transpose times r

% Complex residual vector
ru =  randn(size(FEM.y)) + 1i*randn(size(FEM.y)) ;



 Jc = Jgreen;

 % Real valued transpose of  Jacobian times residual vector
Jtr1= [real(Jgreen), imag(Jc);real(1i*Jgreen), imag(1i*Jc);]'*[real(ru);imag(ru)];

% Transpose of Jacobian operation
[Jtr2] = func_Jtop_C(FEM,ru);


figure
    plot(Jtr1,'b')
    hold on
    plot([real(Jtr2);imag(Jtr2)],'--r')
    




