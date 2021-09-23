%% computes scaled tet volume.
%% V is n x 3 x 4
%SJ, and SJgrad are scaled tet vol and its grad
%J, and Jgrad are unscaled tet vol and its grad. same as face_planarity.
% SJgrad_simple is SJgrad but we pretend edgelengths are independent of vertex pos. Makes the denom essentially a constant scaling factor. 
% my intuition says that wont work as a good proxy because it incentivizes increase jacobian by uniform scaling up, which ignores the relative angles entirely. 
function [SJ,SJgrad, J,Jgrad, SJgrad_simple] = scaledJacobian_vec(V)
    if nargin==0
        file_name = 'meshes/hex_ellipsoid_coarse.vtk';
        [dname,fname,ext] = fileparts(file_name);
        if strcmp(ext,'.vtk')
            mesh = load_vtk(file_name);
        elseif strcmp(ext,'.mesh')
            mesh = ImportHexMesh(file_name);
        end
        V0 = mesh.points; V = V0; nV = size(V,1);
        H = mesh.cells; nH = size(H,1);
        
        HV = permute(reshape(V(H(:),:),nH,8,3),[2 3 1]); % 8 3 nH;
        
        corneredgeinds = [1 3 4 5; 2 6 3 1; 3 2 7 4; 4 1 3 8; 5 1 8 6; 6 5 7 2; 7 6 8 3; 8 5 4 7]; % 8 x 4;
        THV = HV(corneredgeinds,:,:) % 8*4 3 nH
        V = reshape(permute(reshape(THV,8,4,3,[]),[4 1 3 2]),[],3,4); %nH*8 3 4
    end
    
    if nargout > 2
        [J, Jgrad] = face_planarity(V);
        [denom, denomgrad] = denomfunc(V);
    else
        [J] = face_planarity(V);
        [denom] = denomfunc(V);
    end
    SJ = J./denom;
    
    if nargout > 2
        SJgrad = (denom .* Jgrad - J .* denomgrad)./denom.^2;
        SJgrad_simple = Jgrad./denom;    
    end
end


function verifySJ
    V0 = randn(10,3,4);
    pert = randn(size(V0));
    eps = 1e-5;
    [~, dgrad] = scaledJacobian_vec(V0);
    [dp] = scaledJacobian_vec(V0 + eps*pert);
    [dm] = scaledJacobian_vec(V0 - eps*pert);
    fdiff = (dp-dm)/(2*eps);
    adiff = sum(pert.*dgrad,[2 3]);
    assert(norm(fdiff - adiff)<.001);
end

function [denom, denomgrad] = denomfunc(V)
    v1 = V(:,:,1);
    v2 = V(:,:,2);
    v3 = V(:,:,3);
    v4 = V(:,:,4);
    e12 = v1-v2;
    e13 = v1-v3;
    e14 = v1-v4;
    el12 = vecnorm(e12,2,2);
    el13 = vecnorm(e13,2,2);
    el14 = vecnorm(e14,2,2);
    denom = el12.*el13.*el14;
    
    if nargout > 1
        grad_el12 = zeros(size(V));
        grad_el12(:,:,1) = e12./(el12);
        grad_el12(:,:,2) = -e12./(el12);

        grad_el13 = zeros(size(V));
        grad_el13(:,:,1) = e13./(el13);
        grad_el13(:,:,3) = -e13./(el13);

        grad_el14 = zeros(size(V));
        grad_el14(:,:,1) = e14./(el14);
        grad_el14(:,:,4) = -e14./(el14);

        denomgrad = el12.*el13.*grad_el14 + el12.*el14.*grad_el13 + el13.*el14.*grad_el12;
    end
end

function verifydenom
    V0 = randn(10,3,4);
    pert = randn(size(V0));
    eps = 1e-5;
    [~, dgrad] = denomfunc(V0);
    [dp] = denomfunc(V0 + eps*pert);
    [dm] = denomfunc(V0 - eps*pert);
    fdiff = (dp-dm)/(2*eps);
    adiff = sum(pert.*dgrad,[2 3]);
    assert(norm(fdiff - adiff)<.001);
end



