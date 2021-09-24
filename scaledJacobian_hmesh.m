% input hex mesh as V,H.
% input p which is a norm parameter. 1 means uniform sum over all scaled jacobians. higher p means prioritize the worst performing element.
% output is the aggregated scaled jacobian score and gradient perturbation on V.
function [energy, grad,out] = scaledJacobian_hmesh(V,H,p)
    if nargin==0
        file_name = 'meshes/hex_ellipsoid_coarse.vtk';
        file_name = 'test.vtk';
        [dname,fname,ext] = fileparts(file_name);
        if strcmp(ext,'.vtk')
            mesh = load_vtk(file_name);
        elseif strcmp(ext,'.mesh')
            mesh = ImportHexMesh(file_name);
        end
        V = mesh.points;
        H = mesh.cells; 
        p=1;
    end
    nH = size(H,1);
    nV = size(V,1);
    
    HV = permute(reshape(V(H(:),:),nH,8,3),[2 3 1]); % 8 3 nH;
    
    % hardcoded ordering in reference to canonical hex vert ordering. each quadruplet is the central vertex and then right hand rule points into the cube cycling around its 3 adjacent vertices.
    corneredgeinds = [1 2 4 5; 2 6 3 1; 3 2 7 4; 4 1 3 8; 5 1 8 6; 6 5 7 2; 7 6 8 3; 8 5 4 7]; % 8 x 4;
    THV = HV(corneredgeinds,:,:); % 8*4 3 nH
    rV = reshape(permute(reshape(THV,8,4,3,[]),[4 1 3 2]),[],3,4); %nH*8 3 4

    if nargout > 1
        [SJ,SJgrad, J,Jgrad] = scaledJacobian_vec(rV);
    else
        SJ = scaledJacobian_vec(rV);
    end
    
    % minimizing energy will increase all scaled jacobians. p=1 uniformly weights all hexes. higher p gives more weight to worse hexes.
    energy = sum((1- SJ).^p);
    
    if nargout > 1
        out.SJ = SJ;
        out.SJgrad = SJgrad;
        out.J = J;
        out.Jgrad = Jgrad;
        % min over corners per hex first.
        hexalabSJ = min(reshape(SJ,[],8)')'; % [min(hexalabSJ) max(hexalabSJ)];
        out.hexalabSJ = hexalabSJ;

        % reweight and aggregate gradient
        fracturedGrad = reshape(permute((p*(1-SJ).^(p-1)).*-SJgrad, [1 3 2]), nH,8,4,3); % nH,8,4,3
        accuminds = reshape(H(:,corneredgeinds),nH,8,4);
        grad = zeros(size(V));
        grad(:,1) = accumarray(accuminds(:),reshape(fracturedGrad(:,:,:,1),[],1),[nV,1]);
        grad(:,2) = accumarray(accuminds(:),reshape(fracturedGrad(:,:,:,2),[],1),[nV,1]);
        grad(:,3) = accumarray(accuminds(:),reshape(fracturedGrad(:,:,:,3),[],1),[nV,1]);
    end

end

function verify
    file_name = 'meshes/hex_ellipsoid_coarse.vtk';
    mesh = load_vtk(file_name);
    V = mesh.points; V0=V;
    H = mesh.cells; 
    
    % finite diff verify
    p=4; eps=1e-6;
    pert = randn(size(V));
    [~, grad] = scaledJacobian_hmesh(V,H,p);
    [ep] = scaledJacobian_hmesh(V+eps*pert,H,p);
    [em] = scaledJacobian_hmesh(V-eps*pert,H,p);
    fdiff = (ep-em)/(2*eps);
    adiff = dot(grad(:),pert(:));
    [adiff-fdiff]
    
    %% grad desc verify
    V=V0; p=1; clear energy; Vprev = V;
    figure; hold all; axis equal; rotate3d on; axis off;
    for i=1:100
        data = processhmesh(V,H,0);
        try; delete(ptc); catch; end;
        ptc = patch('vertices',V,'faces',data.F,'facealpha',.1,'facecolor','green');
        drawnow;
        
        if i==1; dt = min(data.edgelengths)/300; end;
        [energy(i), grad] = scaledJacobian_hmesh(V,H,p);
        V = V - dt*grad;
        
        % adaptive timestep
        if i>1
            if energy(i) > energy(i-1)
                dt = dt/2;
                V = Vprev;
            else
                dt = dt * 2^(1/4);
                Vprev = V;
            end
        end
    end
end



