% Code for the OIPE Paper
% 
% Precomputation steps for an exemplary Complex ECT sensor
%
% EMS 2024
% Contact: neumayer@tugraz.at
clear all, close all, clc

addpath .\lib

load MESH

figure, hold on
for ii = 1:size(MESH.Inzid,1)
    inz  = MESH.Inzid(ii,:);
    XY = MESH.Node(inz,:);
    XY = [XY;XY(1,:)];
    plot(XY(:,1),XY(:,2),'k');    
    
end
axis square


m = unique(MESH.Material);

id = find(MESH.Material == 0);
for ii = 1:length(id)
    inz  = MESH.Inzid(id(ii),:);
    XY = MESH.Node(inz,:);
    XY = mean(XY);
    plot(XY(:,1),XY(:,2),'.r');    
    
end
colormap gray

return
tmp = zeros(size(MESH.Material));
id = find(MESH.Material == 0); tmp(id) = 150;     %Backside 
id = find(MESH.Material == 150); tmp(id) = 50;    %Pipe 

MESH.Material = tmp;



id = find(MESH.Material == 0);
for ii = 1:length(id)
    inz  = MESH.Inzid(id(ii),:);
    XY = MESH.Node(inz,:);
    XY = mean(XY);
    plot(XY(:,1),XY(:,2),'.g');    
    
end
index_in_matval = find(MESH.Material == 0);
%L = fem_make_reg_matrix_1st_order(MESH.Inzid,index_in_matval,MESH.Node);