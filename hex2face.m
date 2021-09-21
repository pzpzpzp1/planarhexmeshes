function [F, H2F] = hex2face(H)
%     Fall = [H(:,1),H(:,2),H(:,3),H(:,4); H(:,5),H(:,1),H(:,2),H(:,6);...
%           H(:,8),H(:,5),H(:,6),H(:,7); H(:,4),H(:,3),H(:,7),H(:,8);...
%           H(:,6),H(:,2),H(:,3),H(:,7); H(:,5),H(:,1),H(:,4),H(:,8)];
      
    Fall = reshape(permute(reshape(H(:,[4 3 2 1, 5 1 2 6 8 7 3 4 5 6 7 8 7 6 2 3 5 8 4 1]),[],4,6),[1 3 2]),[],4);
      
    [~,ia,ic] = unique(sort(Fall,2),'rows');
    F = Fall(ia,:);
    nH = size(H,1);
    nF = size(F,1);
    
    finds = repmat(1:nF,1,4);
    hfinds = reshape(finds(ic),[],6);
    H2F = sparse(repmat(1:nH,1,6), hfinds, hfinds*0+1, nH, nF);
end