clear all; close all; clc;

%% load mesh
file_name = 'meshes/bunny.vtk';
mesh = load_vtk(file_name);
V0 = mesh.points; V = V0; nV = size(V,1);
H = mesh.cells; nH = size(H,1);

%% Get faces and face extractor
F = hex2face(H);
[Fu,ia,ic] = unique(sort(F,2),'rows');
F = F(ia,:);
nF = size(F,1);
Fmat = sparse(repmat(1:nF,1,4),F,F*0+1,nF,nV);
Fmatselector = sparse(1:(4*nF), F', F*0+1, 4*nF, nV);

figure; t=tiledlayout(1,2);t.TileSpacing = 'compact';t.Padding = 'compact';
maxiters = 10000;
for i=1:maxiters
    %% get pre preprojection stats
    Fv = reshape(Fmatselector*V,4,nF,3); % [4 nf 3]
    Fv_pr = permute(Fv,[2 3 1]);
    preprojection_planarity(:,i) = face_planarity(Fv_pr);
    
    %% local projection and stats
    FvProj_pr = planar_face_projections(Fv_pr, false);
    FvProj = permute(FvProj_pr, [3 1 2]);
    postprojection_planarity(:,i) = face_planarity(FvProj_pr);

    %% aggregate by least squares system
    V = Fmatselector\reshape(FvProj,[],3);
    globalproj_error(i) = norm(Fmatselector*V - reshape(Fv,[],3));
    
    %% visualize
    try; delete(ptc); delete(ptc2); catch; end;
    nexttile(1); hold all; axis equal off; rotate3d on; 
    ptc = patch('faces', reshape(1:(nF*4),4,[])', 'vertices', reshape(FvProj,[],3), 'facecolor','green','facealpha',1,'edgealpha',1)
    title(['Planarity: ' num2str(globalproj_error(i)) ])
    nexttile(2); hold all; axis equal off; rotate3d on; 
    ptc2 = patch('vertices',V,'faces',F,'facealpha',1,'facecolor','green','edgealpha',1)
    title(['Planarity: ' num2str(norm(preprojection_planarity(:,i))) ])
    drawnow;
    
end

figure; t=tiledlayout(3,1);t.TileSpacing = 'compact';t.Padding = 'compact';
nexttile; hold all; title('pre projection');
plot(vecnorm(preprojection_planarity,2,1))
nexttile; hold all; title('post projection');
plot(vecnorm(postprojection_planarity,2,1))
nexttile; hold all; title('global projection error');
plot(globalproj_error)









