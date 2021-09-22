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
    data.V=V;
    data.H=H;
    data.nV = nV;
    data.nH = nH;
    
    %% hex-face
    [F, H2F, H2Farray, F2Harray, H2F6, H2F_flip, Fall] = hex2face(H);
    nF = size(F,1);
    data.nF = nF;
    data.F = F; data.H2F = H2F;
    data.H2Farray = H2Farray;
    data.F2Harray = F2Harray;
    data.H2F6 = H2F6;
    data.H2F_flip = H2F_flip;
    data.Fall = Fall; % (nH*6, 4)
    
    %% hex-edge
    [E, H2E] = hex2edge(H); 
    data.E = E; data.H2E = H2E;
    %% edge-face
    nE = size(E,1);
    faceedges = reshape(permute(reshape(F(:,[1 2, 2 3, 3 4, 4 1]),nF,2,4),[1 3 2]),nF*4,2); % nF 4 2
    [alltrue, faceInEdgelist] = ismember(sort(faceedges,2), sort(E,2),'rows');
    assert(all(alltrue));
    E2F = sparse(repmat(1:nF,1,4), faceInEdgelist, faceInEdgelist*0+1, nF, nE)';
    F2Earray = reshape(faceInEdgelist,nF,4);
    data.nE = nE;
    data.E2F = E2F;
    data.F2Earray = F2Earray;
    %% edge-vert   hex-vert
    E2V = sparse([1:nE, 1:nE],E,E*0+1,nE,nV);
    H2V = sparse(repmat(1:nH,1,8),H,H*0+1,nH,nV);
    F2V = sparse(repmat(1:nF,1,4),F,F*0+1,nF,nV);
    data.E2V = E2V;
    data.H2V = H2V;
    data.F2V = F2V;
    
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
    data.efdeg = efdeg;
    
    intSingEdgesPerVert = sum(E2V(isSingularEdge & ~isBoundaryEdge,:))';
    boundarySingEdgesPerVert = sum(E2V(isSingularEdge & isBoundaryEdge,:))';
    isSingularNode = isSingularVertex & ...
                    ((~isBoundaryVertex & intSingEdgesPerVert~=2) | ...
                    (isBoundaryVertex & ~((intSingEdgesPerVert==1 & boundarySingEdgesPerVert==0 | ...
                                           intSingEdgesPerVert==0 & boundarySingEdgesPerVert==2))));
    data.isSingularEdge = isSingularEdge;
    data.isSingularVertex = isSingularVertex;
    data.isSingularNode = isSingularNode;
    
    %% geometric stuff
    cellBarycenters = H2V*V/8;
    faceBarycenters = F2V*V/4;
    edgelengths = vecnorm(V(E(:,1),:)-V(E(:,2),:),2,2);
    data.cellBarycenters = cellBarycenters;
    data.faceBarycenters = faceBarycenters;
    data.edgelengths = edgelengths;
    
    fv1 = V(F(:,1),:);
    fv2 = V(F(:,2),:);
    fv3 = V(F(:,3),:);
    fv4 = V(F(:,4),:);
    v123n = @(v1,v2,v3) cross(v2-v1,v3-v1);
    n1 = v123n(fv1,fv2,fv3);
    n2 = v123n(fv2,fv3,fv4);
    n3 = v123n(fv3,fv4,fv1);
    n4 = v123n(fv4,fv1,fv2);
    n = (n1+n2+n3+n4)/4;
    n = n./vecnorm(n,2,2);
    faceNormals = n;
    data.faceNormals = faceNormals;
    
    bvertNormals = F2V(isBoundaryFace,isBoundaryVertex)' * data.faceNormals(isBoundaryFace,:);
    bvertNormals = bvertNormals./vecnorm(bvertNormals,2,2);
    data.bvertNormals = bvertNormals;
    
    %% visualize
    if visualize
        figure; hold all; axis equal off; rotate3d on;
        patch('vertices',V,'faces',F(isBoundaryFace,:),'facecolor','green','facealpha',.1,'edgealpha',0)
        patch('vertices',V,'faces',E(isBoundaryEdge,[1 2 1]))
        patch('vertices',V,'faces',E(isSingularEdge,[1 2 1]),'linewidth',3,'edgecolor','blue')
        scatter3(V(isSingularNode,1),V(isSingularNode,2),V(isSingularNode,3),100,'r','filled')
        
%         quiver3(faceBarycenters(isBoundaryFace,1),faceBarycenters(isBoundaryFace,2),faceBarycenters(isBoundaryFace,3),...
%             faceNormals(isBoundaryFace,1),faceNormals(isBoundaryFace,2),faceNormals(isBoundaryFace,3))
%         quiver3(V(isBoundaryVertex,1),V(isBoundaryVertex,2),V(isBoundaryVertex,3),...
%             bvertNormals(:,1),bvertNormals(:,2),bvertNormals(:,3))
    end
end