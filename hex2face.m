function [F, H2F, H2Farray, F2Harray, H2F6, H2F_flip, Fall] = hex2face(H)
%     Fall = [H(:,1),H(:,2),H(:,3),H(:,4); H(:,5),H(:,1),H(:,2),H(:,6);...
%           H(:,8),H(:,5),H(:,6),H(:,7); H(:,4),H(:,3),H(:,7),H(:,8);...
%           H(:,6),H(:,2),H(:,3),H(:,7); H(:,5),H(:,1),H(:,4),H(:,8)];
    nH = size(H,1);
    Fall = reshape(permute(reshape(H(:,[4 3 2 1, 5 1 2 6, 8 7 3 4, 5 6 7 8, 7 6 2 3, 5 8 4 1]),[],4,6),[1 3 2]),[],4); % (nH*6, 4)
    Fall_rs = reshape(Fall,[],6,4);
    H2F_flip = [4 3 2 1 6 5];
    
    hexfaceinds = repmat([1 2 3 4 5 6]', 1, nH)';
    hexfaceindsud = flipud(hexfaceinds(:));
    fallfinds = repmat((1:nH)',1,6);
    
    [~,ia,ic] = unique(sort(Fall,2),'rows');
    F = Fall(ia,:);
    nF = size(F,1);
    
    finds = repmat(1:nF,1,4);
    H2Farray = reshape(finds(ic),[],6);
    H2F = sparse(repmat(1:nH,1,6), H2Farray, H2Farray*0+1, nH, nF);
    
    H2F6 = sparse(repmat(1:nH,1,6), H2Farray, hexfaceinds, nH, nF);
    [ii,jj,kk] = find(H2F6);
    for i=1:numel(ii)
        A = squeeze(Fall_rs(ii(i),kk(i),:))';
        B = F(jj(i),:);
        assert(numel(intersect(A,B))==4)
    end
    
    %% build F2Harray. for boundary, first column indicates which hex the face is part of. second column indicates which face of that hex it is as 1-6. same for interior but now theres two hexes.
    isBoundaryFace = full(sum(H2F)'==1);
    F2Harray = zeros(nF,4);
    [ii,jj] = find(H2F(:,isBoundaryFace)'); [ii,perm] = sort(ii); jj=jj(perm);
    F2Harray(isBoundaryFace,1) = jj;
    [ii,jj] = find(H2F(:,~isBoundaryFace)'); [ii,perm] = sort(ii); jj=jj(perm);
    F2Harray(~isBoundaryFace,[1 3]) = reshape(jj,2,[])';
    F2Harray(isBoundaryFace,2) = H2F6(sub2ind([nH,nF], F2Harray(isBoundaryFace,1), find(isBoundaryFace)));
    F2Harray(~isBoundaryFace,2) = H2F6(sub2ind([nH,nF], F2Harray(~isBoundaryFace,1), find(~isBoundaryFace)));
    F2Harray(~isBoundaryFace,4) = H2F6(sub2ind([nH,nF], F2Harray(~isBoundaryFace,3), find(~isBoundaryFace)));
    
end