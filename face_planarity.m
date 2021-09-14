% quantifies planarity of a face with vertices V: [nF 3 4]
% essentially computes unsigned volume of the tet formed by 4 verts.
function [p1, pgrad] = face_planarity(V)
    if nargin ==0
        nF = 100;
        V = randn(nF,3,4);
    end
    nF = size(V,1);
    
    v1 = V(:,:,1);
    v2 = V(:,:,2);
    v3 = V(:,:,3);
    v4 = V(:,:,4);
    
    c1 = cross(v2-v1, v3-v1);
    c2 = cross(v3-v2, v4-v2);
    c3 = cross(v4-v3, v1-v3);
    c4 = cross(v1-v4, v2-v4);
    svol1 = dot(c1, v4-v1, 2);
%     svol2 = dot(c2, v1-v2, 2);
%     svol3 = dot(c3, v2-v3, 2);
%     svol4 = dot(c4, v3-v4, 2);
%     p1 = abs(svol1);
    p1 = svol1;
%     p2 = abs(svol2);
%     p3 = abs(svol3);
%     p4 = abs(svol4);
    i1 = sign(svol1);
%     i2 = sign(svol2);
%     i3 = sign(svol3);
%     i4 = sign(svol4);
    
    pgrad = zeros(nF,3,4);
%     pgrad(:,:,4) = i1 .* c1;
%     pgrad(:,:,1) = -i1 .* c2;
%     pgrad(:,:,2) = i1 .* c3;
%     pgrad(:,:,3) = -i1 .* c4;
    
    pgrad(:,:,4) = c1;
    pgrad(:,:,1) = -c2;
    pgrad(:,:,2) = c3;
    pgrad(:,:,3) = -c4;
end

%% finite diff check
function verify()
    nF = 100;
    V = randn(nF,3,4);
    [p1, pgrad] = face_planarity(V);
    dV = randn(nF,3,4);
    eps = 1e-3;
    p1p = face_planarity(V+eps*dV);
    p1m = face_planarity(V-eps*dV);
    fdiff = (p1p-p1m)/(2*eps);
    adiff = dot(reshape(pgrad,nF,[]), reshape(dV,nF,[]), 2);
    norm([fdiff - adiff])

end