function planarize_hmesh_fmincon(file_name)

    if nargin==0
        % file_name = 'meshes/bunny.vtk';
%         file_name = 'meshes/double-torus.vtk';
        % file_name = 'meshes/joint.vtk';
        % file_name = 'meshes/rockarm.vtk';
        % file_name = 'meshes/hex_sphere.vtk';
        % file_name = 'meshes/hex_tetrahedron.vtk';
        file_name = 'meshes/hex_ellipsoid_coarse.vtk';
        % file_name = 'meshes/kitten.mesh';
    end

    %% load mesh
    [dname,fname,ext] = fileparts(file_name);
    if strcmp(ext,'.vtk')
        mesh = load_vtk(file_name);
    elseif strcmp(ext,'.mesh')
        mesh = ImportHexMesh(file_name);
    end
    V0 = mesh.points; V = V0; nV = size(V,1);
    H = mesh.cells; nH = size(H,1);

    %% Get faces and face extractor
    F = hex2face(H);
    [Fu,ia,ic] = unique(sort(F,2),'rows');
    F = F(ia,:);
    nF = size(F,1);
    Fmat = sparse(repmat(1:nF,1,4),F,F*0+1,nF,nV);
    Fmatselector = sparse(1:(4*nF), F', F*0+1, 4*nF, nV);
    dat.H = H; 
    dat.F = F;
    dat.Fmat = Fmat;
    dat.Fmatselector = Fmatselector;

    %% setup fmincon
    options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'HessianApproximation','lbfgs',...
        'StepTolerance',1e-10,'ConstraintTolerance',1e-5,'Display','iter','OptimalityTolerance',1e-6 ,...
        'CheckGradients',true,'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-3);
    fun = @(x) objfun(V0, x);
    con = @(x) mycon(dat, x);
    [Vproj,fval,exitflag,output,lambda,~,~] = fmincon(fun, V0, [],[],[],[],[],[],con,options);
    

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
        Vs{i} = V;

        %% visualize
        try; delete(ptc); delete(ptc2); catch; end;
        nexttile(1); hold all; axis equal off; rotate3d on; 
        ptc = patch('faces', reshape(1:(nF*4),4,[])', 'vertices', reshape(FvProj,[],3), 'facecolor','green','facealpha',1,'edgealpha',1)
        title(['Glob Proj Error: ' num2str(globalproj_error(i)) ])
        nexttile(2); hold all; axis equal off; rotate3d on; 
        ptc2 = patch('vertices',V,'faces',F,'facealpha',1,'facecolor','green','edgealpha',1)
        title(['Planarity: ' num2str(norm(preprojection_planarity(:,i))) ])
        drawnow;

        if i > 30
            last30 = vecnorm(preprojection_planarity(:,(i-30):i),2,1);
            [p,S,mu] = polyfit(1:numel(last30), last30, 1);
            if p(1) > -convergencethresh % energy decrease slope is too small
                disp('Converged');
                break;
            end
        end
    end

    %% plot results. seems to be logarithmic convergence.
    figure; t=tiledlayout(3,1);t.TileSpacing = 'compact';t.Padding = 'compact';
    nexttile; 
    semilogx(vecnorm(preprojection_planarity(:,1:i),2,1))
    hold all; title('pre projection');
    nexttile; 
    plot(vecnorm(postprojection_planarity(:,1:i),2,1))
    hold all; title('post projection');
    nexttile; 
    cla; semilogx(globalproj_error(1:i))
    hold all; title('global projection error');

    %% save output
    [dname,fname,ext] = fileparts(file_name);
    out_fname = ['results/' fname];
    outmesh.cells = mesh.cells; outmesh.points = V;
    save_vtk(outmesh, [out_fname '.vtk'])
    save([out_fname '.mat'])

end

function [c, ceq, gradc, gradceq] = mycon(dat, V)
    %% get pre preprojection stats
    F = dat.F;
    nF = size(F,1);
    nV = size(V,1);
    Fv = reshape(dat.Fmatselector*V,4,nF,3); % [4 nf 3]
    Fv_pr = permute(Fv,[2 3 1]);
    [ceq, gradceq_parts] = face_planarity(Fv_pr);
    
    ii = repelem((1:nF),12,1);
%     jj = repmat(F,1,3)';
    jj = [F,F+nV,F+2*nV]';
    kk = permute(gradceq_parts,[3 2 1]);
    gradceq = sparse(ii(:), jj(:), kk(:), nF, nV*3)';
    
    c=0;
    gradc = zeros(nV*3,1);
end

function [val, grad] = objfun(V,x)
    grad = -(V(:)-x(:));
    val = .5 * grad' * grad;
end




