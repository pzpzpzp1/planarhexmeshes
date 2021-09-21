function [Vo,Ho] = minimizeMesh(Vo,Ho)
    fullVerts = 1:size(Vo,1);
    referencedVerts = unique(Ho(:));
    unreferencedVerts = setdiff(fullVerts, referencedVerts);
    
    inds = 0*fullVerts;
    inds(referencedVerts) = 1:numel(referencedVerts);
    Ho = inds(Ho);
    Vo(unreferencedVerts,:) = [];
end