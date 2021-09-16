clear all; close all; clc;
% file_name = 'meshes/double-torus.vtk';
% file_name = 'meshes/joint.vtk';
% file_name = 'meshes/rockarm.vtk';
% file_name = 'meshes/hex_ellipsoid_coarse.vtk';
% file_name = 'meshes/hex_tetrahedron.vtk';
% file_name = 'meshes/hex_sphere.vtk';
% file_name = 'meshes/kitten.mesh';
% file_name = 'meshes/bunny.vtk';
file_name = 'meshes/hex_ellipsoid_coarse.vtk';

[dname,fname,ext] = fileparts(file_name);

if strcmp(ext,'.vtk')
    inmesh = load_vtk(file_name);
elseif strcmp(ext,'.mesh')
    inmesh = ImportHexMesh(file_name);
end

out_fname = ['results_fmincon/' fname];
outmesh = load_vtk([out_fname '.vtk']);
assert(all(inmesh.cells(:)==outmesh.cells(:)))

H = inmesh.cells;
F = hex2face(H); [Fu,ia,ic] = unique(sort(F,2),'rows'); F = F(ia,:);

%% visualize linear interpolation
N = 100;
t = linspace(0,1,N);
figure; 
hold all; axis equal off; rotate3d on;
while true
    for i=1:N
        Vt = outmesh.points*t(i) + inmesh.points*(1-t(i));
        try; delete(ptc); catch; end;
        ptc = patch('vertices',Vt,'faces',F,'facealpha',1,'facecolor','green','edgealpha',1)
        title(sprintf('%d out of %d',i,N));
        drawnow;
    end
end





