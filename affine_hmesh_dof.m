close all; clear all; clc;

% file_name = 'results_fmincon/hex_ellipsoid_coarse.vtk';
% file_name = 'results_fmincon/sing1.vtk';
% file_name = 'results_fmincon/sing2.vtk';
% file_name = 'results_fmincon/sing3.vtk';
% file_name = 'meshes/blood_vessel.mesh';
file_name = 'meshes/bunny.vtk';
% file_name = 'meshes/metatron.mesh';
% file_name = 'meshes/pinion.mesh';
% file_name = 'meshes/torus.mesh';
% file_name = 'meshes/sphinx.mesh';
    
%% Load mesh
[dname,fname,ext] = fileparts(file_name);
if strcmp(ext,'.vtk')
    mesh = load_vtk(file_name);
elseif strcmp(ext,'.mesh')
    mesh = ImportHexMesh(file_name);
end
V = mesh.points; H = mesh.cells; 
% V = randn(8,3); H=1:8;
nV = size(V,1); nH = size(H,1);
[F, H2F] = hex2face(H); 
[E, H2E] = hex2edge(H); 

figure; hold all; axis equal; rotate3d on;
patch('vertices',V,'faces',F,'facealpha',.1,'facecolor','green');

%% build selector matrices
U = sparse([[0,0,1]; [0,1,1]; [0,1,0];[0,0,0]; [1,0,1]; [1,1,1];[1,1,0]; [1,0,0];]);

ii = repmat(1:8,1,nH);
jj = repelem(1:nH,8,1);
kk = H';
ll = kk*0+1;
iijj = sub2ind([8 nH], ii(:), jj(:));
Hselector = sparse(iijj,kk,ll);
RT = sparse((1:(8*nH))', repelem((1:nH)',8,1), 1,8*nH,nH);
RU = kron(speye(nH),U);

Hv = Hselector*V;


%% use cvx to project
cvx_begin
    cvx_solver mosek
    variable A(3*nH,3)
    variable t(nH,3)
    variable Vo(nV,3)
    AU = reshape(permute(reshape(A*U',3,nH,8),[3,2,1]),8*nH,3);
    minimize norm(V-Vo,'fro')
    subject to
        RU*A + RT*t == Hselector*Vo
cvx_end

%% visualize projection
ft = figure; 
ax1 = subplot(1,2,1); hold all; axis equal; rotate3d on; axis off;
ptc = patch('vertices',V,'faces',F,'facealpha',.1,'facecolor','green');
ax2 = subplot(1,2,2); hold all; axis equal; rotate3d on; axis off;
ptc = patch('vertices',Vo,'faces',F,'facealpha',.1,'facecolor','green');
ft.UserData = linkprop([ax1,ax2],{'xlim','ylim','zlim','CameraPosition','CameraTarget'});

%% get dof of system.
norm(Hselector*Vo - RU*A - RT*t,'fro')
sysA = full([Hselector, RU, RT]);
nsp = null(sysA); % num cols is dofs per coordinate. X3 for full dofs.
dofs = size(nsp,2)*3;
% verified dofs with singular vertices. num interior edges x3 + 3 for translation of center.

%% visualize dofs
figure; hold all; axis equal; rotate3d on; axis off;
T=200; ts = linspace(0,2,T)*2*pi;
for i=1:size(nsp,2)
    for j=1:3
        dV = zeros(size(V));
        dV(:,j) = nsp(1:nV,i);
        
        for ti=1:T
            Vt = Vo + dV*sin(ts(ti));
            try; delete(ptc); delete(qvr); catch; end;
            ptc = patch('vertices',Vt,'faces',F,'facealpha',.1,'facecolor','green');
            qvr = quiver3(Vt(:,1),Vt(:,2),Vt(:,3),dV(:,1),dV(:,2),dV(:,3),'r','linewidth',2)
            drawnow; pause(.0001)
            title(num2str(ts(ti)))
        end
        pause
    end
end


figure; hold all; axis equal; rotate3d on; axis off;
T=200; ts = linspace(0,1,T);
for j=1:10
    % patch('vertices',V,'faces',F,'facealpha',.1,'facecolor','red','edgealpha',.1);
    for i=1:T
        Vt = Vo*ts(i) + V*(1-ts(i));
        try; delete(ptc); catch; end;
        ptc = patch('vertices',Vt,'faces',F,'facealpha',.1,'facecolor','green');
        drawnow; pause(.0001)
    end
end













