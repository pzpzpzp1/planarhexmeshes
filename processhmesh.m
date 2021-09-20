function data = processhmesh(V,H,visualize)
    if nargin ==0
        file_name = 'results_fmincon/hex_ellipsoid_coarse.vtk';
        mesh = load_vtk(file_name);
        % mesh = ImportHexMesh(file_name);
        V = mesh.points;
        H = mesh.cells;
        visualize = 1;
    end
    nV = size(V,1); 
    nH = size(H,1);
    
    %% hex-face
    [F, H2F] = hex2face(H); 
    nF = size(F,1);
    data.F = F; data.H2F = H2F;
    %% hex-edge
    [E, H2E] = hex2edge(H); 
    data.E = E; data.H2E = H2E;
    %% edge-face
    nE = size(E,1);
    faceedges = reshape(permute(reshape(F(:,[1 2, 2 3, 3 4, 4 1]),nF,2,4),[1 3 2]),nF*4,2); % nF 4 2
    [alltrue, faceInEdgelist] = ismember(sort(faceedges,2), sort(E,2),'rows');
    assert(all(alltrue));
    E2F = sparse(repmat(1:nF,1,4), faceInEdgelist, faceInEdgelist*0+1, nF, nE)';
    data.E2F = E2F;
    %% edge-vert   hex-vert
    E2V = sparse([1:nE, 1:nE],E,E*0+1,nE,nV);
    H2V = sparse(repmat(1:nH,1,8),H,H*0+1,nH,nV);
    data.E2V = E2V;
    data.H2V = H2V;
    
    %% get boundary
    isBoundaryFace = full((sum(H2F)==1)');
    isBoundaryEdge = full(sum(E2F(:,isBoundaryFace)')'~=0);
    isBoundaryVertex = false(nV,1); isBoundaryVertex(unique(F(isBoundaryFace,:)))=true;
    isBoundaryHex = any(isBoundaryVertex(H),2);
    data.isBoundaryFace = isBoundaryFace;
    data.isBoundaryEdge = isBoundaryEdge;
    data.isBoundaryVertex = isBoundaryVertex;
    data.isBoundaryHex = isBoundaryHex;
    
    %% get singularities
    efdeg = full(sum(E2F,2)); % edge face degree
    isSingularEdge = full((efdeg~=4 & ~isBoundaryEdge) | (efdeg~=3 & isBoundaryEdge));
    isSingularVertex = false(nV,1); isSingularVertex(unique(E(isSingularEdge,:)))=true;
    
    intSingEdgesPerVert = sum(E2V(isSingularEdge & ~isBoundaryEdge,:))';
    boundarySingEdgesPerVert = sum(E2V(isSingularEdge & isBoundaryEdge,:))';
    isSingularNode = isSingularVertex & ...
                    ((~isBoundaryVertex & intSingEdgesPerVert~=2) | ...
                    (isBoundaryVertex & ~((intSingEdgesPerVert==1 & boundarySingEdgesPerVert==0 | ...
                                           intSingEdgesPerVert==0 & boundarySingEdgesPerVert==2))));
    data.isSingularEdge = isSingularEdge;
    data.isSingularVertex = isSingularVertex;
    data.isSingularNode = isSingularNode;
    
    if visualize
        figure; hold all; axis equal off; rotate3d on;
        patch('vertices',V,'faces',F(isBoundaryFace,:),'facecolor','green','facealpha',.1,'edgealpha',0)
        patch('vertices',V,'faces',E(isBoundaryEdge,[1 2 1]))
        patch('vertices',V,'faces',E(isSingularEdge,[1 2 1]),'linewidth',2,'edgecolor','blue')
        scatter3(V(isSingularNode,1),V(isSingularNode,2),V(isSingularNode,3),50,'r','filled')
    end
    
end