close all; clear all; clc;

file_name = 'results_fmincon/hex_ellipsoid_coarse.vtk';
file_name = 'results_fmincon/sing1.vtk';
file_name = 'results_fmincon/sing2.vtk';
file_name = 'results_fmincon/sing3.vtk';
    
%% Load mesh
mesh = load_vtk(file_name);
V = mesh.points; H = mesh.cells; nV = size(V,1); nH = size(H,1);
[F, H2F] = hex2face(H); 
[E, H2E] = hex2edge(H); 
    
%% build selector matrices
U = [[0,0,1]; [0,1,1]; [0,1,0];[0,0,0]; [1,0,1]; [1,1,1];[1,1,0]; [1,0,0];];

ii = repmat(1:8,1,nH);
jj = repelem(1:nH,8,1);
kk = H';
ll = kk*0+1;
iijj = sub2ind([8 nH], ii(:), jj(:));
Hselector = sparse(iijj,kk,ll);

Hv = Hselector*V;

%% get dof of system
[Hv U']
