nF = 20;
V = randn(nF,3,4);


f1 = figure; axis equal; hold all; rotate3d on;
dt = .0001;
for i=1:1000
    cla; axis equal off; hold all; rotate3d on;
    for j=1:nF
        Vj = permute(V(j,:,:),[3 2 1]);
        patch('vertices',Vj,'faces',[1 2 3 4],'facecolor','green')
    end
    drawnow; pause(.01);
    
    [p1, pgrad] = face_planarity(V);
    V = V - dt * pgrad;
    
end