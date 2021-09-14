function planarize_hmesh_fmincon(file_name)

    if nargin==0
%         file_name = 'meshes/bunny.vtk';
%         file_name = 'meshes/double-torus.vtk';
%         file_name = 'meshes/joint.vtk';
        % file_name = 'meshes/rockarm.vtk';
        % file_name = 'meshes/hex_sphere.vtk';
        % file_name = 'meshes/hex_tetrahedron.vtk';
        file_name = 'meshes/hex_ellipsoid_coarse.vtk';
%         file_name = 'meshes/kitten.mesh';
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
        'CheckGradients',false,'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-3);
    fun = @(x) objfun(V0, x);
    con = @(x) mycon(dat, x);
    [Vproj,fval,exitflag,output,lambda,~,~] = fmincon(fun, V0, [],[],[],[],[],[],con,options);
    [~,final_planarity] = con(Vproj);
    
    figure; hist(final_planarity); hold all; title('final planarities');
    
    f1=figure; t=tiledlayout(1,2);t.TileSpacing = 'compact';t.Padding = 'compact';
    ax1 = nexttile(1); hold all; axis equal off; rotate3d on; 
    ptc1 = patch('vertices',V0,'faces',F,'facealpha',1,'facecolor','green','edgealpha',1)
    title('Original')
    ax2 = nexttile(2); hold all; axis equal off; rotate3d on; 
    ptc2 = patch('vertices',Vproj,'faces',F,'facealpha',1,'facecolor','green','edgealpha',1)
    title(['After Projection Distance: ' num2str(fval) ' Planarity: ' num2str(norm(final_planarity)) ])
    f1.UserData=linkprop([ax1 ax2],{'XLim','Ylim','Zlim','CameraTarget','CameraPosition'})
    
    %% save output
    [dname,fname,ext] = fileparts(file_name);
    out_fname = ['results_fmincon/' fname];
    outmesh.cells = mesh.cells; outmesh.points = Vproj;
    save_vtk(outmesh, [out_fname '.vtk'])
    
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




