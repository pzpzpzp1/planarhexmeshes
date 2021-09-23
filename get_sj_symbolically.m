%{ 
THE FINDINGS OF THIS FILE ARE THAT MATLAB SYM DIFF IS TOO SLOW TO USE EVEN WITH VECTORIZED EVALUATION. 
ALSO THAT THE JACDET AT A CORNER IS JUST THE DET OF THE MATRIX WHOSE COLUMNS ARE THE THREE EDGE VECTORS COMING OUT OF THAT CORNER.
WHICH IS PROPORTIONAL OF THAT TET VOLUME. MEANING FACE_PLANARITY ALREADY COMPUTES SOME OF THIS FOR US. JUST GOTTA LEVERAGE IT.
%}

% get scaled jacobian symbolically
% file_name = 'meshes/unit.vtk';
% mesh = load_vtk(file_name);
% V = mesh.points(mesh.cells,:);
clear all; close all;

uvw = [1     0     1;...
    1     1     1;...
    1     1     0;...
    1     0     0;...
    0     0     1;...
    0     1     1;...
    0     1     0;...
    0     0     0];
figure; hold all; rotate3d on; axis equal;
for i=1:8
    text(uvw(i,1),uvw(i,2),uvw(i,3),['<---' num2str(i)])
end

syms v11 v12 v13 ...
     v21 v22 v23 ...
     v31 v32 v33 ...
     v41 v42 v43 ...
     v51 v52 v53 ...
     v61 v62 v63 ...
     v71 v72 v73 ...
     v81 v82 v83 real;
V = [v11 v12 v13 ;...
     v21 v22 v23 ;...
     v31 v32 v33 ;...
     v41 v42 v43 ;...
     v51 v52 v53 ;...
     v61 v62 v63 ;...
     v71 v72 v73 ;...
     v81 v82 v83];
% index mismatch between standard vtk and 'symbolic_jacobian_det'
ia = [6 7 3 2 5 8 4 1];
ib = [8 4 3 7 5 1 2 6];

for i=1:8;
    % J(i,1) = simplify(symbolic_jacobian_det(V([5 6 7 8 1 2 3 4],:), uvw(i,:)));
    J(i,1) = symbolic_jacobian_det(V(ib,:), uvw(i,:));
end
corneredgeinds = reshape([1 2, 1 5, 1 4; 2 1, 2 6, 2 3; 3 4 3 7 3 2; 4 1,4 8,4 3; 5 1 5 6 5 8; 6 7 6 5 6 2; 7 8 7 6 7 3; 8 7 8 5 8 4]',2,3,8);
% ismember(V,symvar(J(1,:)))

%% plot to verify edge connections
figure; hold all; rotate3d on; axis equal off;
for i=1:8
    cla; hold all; rotate3d on; axis equal off;
    for ii=1:8
        text(uvw(ii,1),uvw(ii,2),uvw(ii,3),['<---' num2str(ii)])
    end
    for j=1:3
        vs = uvw(corneredgeinds(:,j,i),:);
        plot3(vs(:,1),vs(:,2),vs(:,3),'k')
    end
    drawnow;
end

diffs = V(corneredgeinds(1,:,:),:)-V(corneredgeinds(2,:,:),:);
elens = reshape(sqrt(sum(diffs.^2,2)),3,8)';
sj_denom = elens(:,1).*elens(:,2).*elens(:,3);

SJ = J./sj_denom;

%% Get derivatives w.r.t. J, SJ
for i=1:24
    gradJ(:,i) = diff(J,V(i));
end
size(gradJ)
for i=1:24
    gradSJ(:,i) = diff(SJ,V(i));
end
size(gradSJ)

out.gradJ = simplify(gradJ);
out.gradSJ = simplify(gradSJ);
out.J=simplify(J); 
out.SJ=simplify(SJ);

%% save symbolic output
figure;spy(gradJ)
for i=1:8
    assert(nnz(gradJ(i,:))==12)
end

%% SCALED TET VOLUME FORMULATION
v1 = V(1,:);
v2 = V(2,:);
v3 = V(4,:);
v4 = V(5,:);
tetvol1 = dot(cross(v2-v1,v3-v1),v4-v1);
tetvol2 = det([v3' v4' v2']-v1');
simplify(tetvol1 - J(1,:))
simplify(tetvol2 - J(1,:))

Vrs = permute(V([1 2 4 5],:),[3 2 1]);
Vvals = randn(1,3,4);
[p1, pgrad] = face_planarity(Vvals);

v11 = Vvals(1,1,1);
v12 = Vvals(1,2,1);
v13 = Vvals(1,3,1);
v21 = Vvals(1,1,2);
v22 = Vvals(1,2,2);
v23 = Vvals(1,3,2);
v41 = Vvals(1,1,3);
v42 = Vvals(1,2,3);
v43 = Vvals(1,3,3);
v51 = Vvals(1,1,4);
v52 = Vvals(1,2,4);
v53 = Vvals(1,3,4);
p1_s = double(subs(J(1)))
pgrad_s = double(subs(gradJ(1,:)))
pgrad_s = reshape(pgrad_s,8,3);
pgrad_s = pgrad_s([1 2 4 5],:);
pgrad_s - permute(pgrad,[2 3 1])'



