function [F, H2F] = hex2face(H)
    Fall = [H(:,1),H(:,2),H(:,3),H(:,4); H(:,5),H(:,1),H(:,2),H(:,6);...
          H(:,8),H(:,5),H(:,6),H(:,7); H(:,4),H(:,3),H(:,7),H(:,8);...
          H(:,6),H(:,2),H(:,3),H(:,7); H(:,5),H(:,1),H(:,4),H(:,8)];
      
    [~,ia,ic] = unique(sort(Fall,2),'rows');
    F = Fall(ia,:);
    nH = size(H,1);
    nF = size(F,1);
    
    finds = repmat(1:nF,1,4);
    hfinds = reshape(finds(ic),[],6);
    H2F = sparse(repmat(1:nH,1,6), hfinds, hfinds*0+1, nH, nF);
end