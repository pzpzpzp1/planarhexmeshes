function stencils = convert_hmesh_to_stencils(file_name, uniformkeys)
    if nargin==0
        close all;
%         file_name = 'results_fmincon/hex_ellipsoid_coarse.vtk';
%         file_name = 'results_fmincon/sing1.vtk';
%         file_name = 'results_fmincon/sing2.vtk';
        file_name = 'results_fmincon/sing3.vtk';
%         file_name = 'results_fmincon/bunny.vtk';
        uniformkeys = true;
    end
    
    %% magic parameters
    pfactor = .5;
    keybufferperc = .2; % reserved space on each corner

    %% Load mesh
    mesh = load_vtk(file_name);
    V = mesh.points; H = mesh.cells; nV = size(V,1); nH = size(H,1);
    [F, H2F] = hex2face(H); 
    nF = size(F,1);
    Fmatselector = sparse(1:(4*nF), F', F*0+1, 4*nF, nV);
    
    %% Test planarity
    Fv = reshape(Fmatselector*V,4,nF,3); % [4 nf 3]
    Fv_pr = permute(Fv,[2 3 1]);
    [ceq] = face_planarity(Fv_pr);
    if(norm(ceq) > 1e-6)
        error(['Faces of this hex mesh are not planar enough. planarity: ' num2str(norm(ceq))]);
    end

    %% Get edge structures and build keys
    [E, H2E] = hex2edge(H); nE = size(E,1);
    faceedges = reshape(permute(reshape(F(:,[1 2, 2 3, 3 4, 4 1]),nF,2,4),[1 3 2]),nF*4,2); % nF 4 2
    [alltrue, faceInEdgelist] = ismember(sort(faceedges,2), sort(E,2),'rows');
    assert(all(alltrue));
    E2F = sparse(repmat(1:nF,1,4), faceInEdgelist, faceInEdgelist*0+1, nF, nE)';
    efdegs = sum(E2F,2); % edge face degree
    elens = vecnorm(V(E(:,1),:)-V(E(:,2),:),2,2);
    
    E2Fc = cell(nE,1);
    for i=1:nE
        deg = efdegs(i);
        keyorder = efdeg_to_keyorder(deg);
        nk = numel(keyorder);
        
        % generate key sizes.
        if uniformkeys
            keylens = ones(nk,1)';
        else
            keylens = rand(nk,1)'+.5;
        end
        keylens = keylens/sum(keylens);
        
        estr.keylens = keylens;
        estr.keyorder = keyorder;
        estr.deg = deg;
        estr.faces = find(E2F(i,:));
        E2Fc{i} = estr;
    end

    %% Build planar face outlines using edge keys
    keybuffer = min(elens) * keybufferperc;
    for i=1:nF;              %         i=randi(nF)
        % build keys
        fvinds = F(i,:);
        for j=1:4
            jjp1 = [j j+1];
            if j==4; jjp1 = [4 1]; end;
            [isforward, eindf] = ismember(fvinds(jjp1), E, 'rows');
            [~,         eindb] = ismember(fvinds(flip(jjp1)), E, 'rows');
            eind = eindf+eindb;
            
            estr = E2Fc{eind};
            keyorder = estr.keyorder;
            keylens = estr.keylens;
            if ~isforward
                keylens = flip(keylens);
                keyorder = flip(keyorder);
            end
            
            iamfacen = find(estr.faces == i);
            [key{j}, keylensA{j}] = buildkey(keyorder,keylens,iamfacen,keybuffer,elens(eind),pfactor);
        end
        
        % fit key to embedded edge
        fv = V(fvinds,:);
        fvout = rotate_face_flat(fv);
        for j=1:4
            jjp1 = [j j+1];
            if j==4; jjp1 = [4 1]; end;
            
            x = key{j}(:,1);
            y = key{j}(:,2);
            % yscaled = y * -keybuffer;
            yscaled = y * -keylensA{j}(2);
            startp = fvout(jjp1(1),:);
            endp = fvout(jjp1(2),:);
            dir = endp-startp;
            angle = atan2(dir(2), dir(1));
            R2 = axang2rotm([0 0 1 angle]);
            line{j} = R2(1:2,1:2)*[x'; yscaled'] + startp';
            % plot(line{j}(1,:),line{j}(2,:),'k.-','linewidth',2)
        end
        
        % combine into polyshape
        xys = zeros(2,0);
        for j=1:4
            xys = [xys line{j}];
        end
        ps = polyshape(xys(1,:),xys(2,:));
        
        pieces{i}.keypiece = ps;
        pieces{i}.vinds = fvinds;
        pieces{i}.quad = fvout;
        
    end
    
    %% plot all pieces together.
    % plot individual pieces flat on one figure.
    colors = uint8(rand(nF,3)*255);
    maxelen = max(elens);
    sqind = ceil(sqrt(nF));
    xx = (1:sqind)*maxelen*sqrt(2);
    [XX,YY] = meshgrid(xx,xx);
    
    figure; hold all; axis equal off;  drawnow;
    for i=1:nF
        piece = pieces{i};
        [ii,jj] = ind2sub([sqind,sqind],i);
        newcenter = [XX(ii,jj) YY(ii,jj)];
        
        fvout = pieces{i}.quad + newcenter;
        pstranslated = pieces{i}.keypiece; pstranslated.Vertices = pstranslated.Vertices + newcenter;
        patch('faces',[1 2 3 4],'vertices',fvout,'facecolor','green','facealpha',.2,'edgealpha',.3);
        plot(pstranslated,'facecolor',double(colors(i,:))/255)
        drawnow;
    end
    
    % plot the full mesh.
    figure; hold all; axis equal off; rotate3d on;
    ptc = patch('faces',F,'vertices',V,'facecolor','flat','facealpha',.3,'edgealpha',1, 'FaceVertexCData', colors);
    drawnow; pause(.1)
    
end

% build key.
function [key, keylensA] = buildkey(keyorder,keylens,iamfacen,keybuffer,elen,pfactor)
    keyorderA = [0 keyorder 0];
    keytotallength = elen-2*keybuffer;
    keylensA = keylens/sum(keylens);
    keylensA = keylensA*keytotallength;
    keylensA = [keybuffer keylensA keybuffer];


    raised = keyorderA == iamfacen;
    nk = numel(keyorderA);
    x=0;
    y=0;
    curr = 0;
    for i=1:nk
        if raised(i)==curr
            x = [x, x(end)+keylensA(i)];
            y = [y, curr];
        else
            x = [x, x(end)];
            y = [y, ~curr];
            
            x = [x, x(end)+keylensA(i)];
            y = [y, y(end)];
        end
        curr = raised(i);
    end
    y=y-.5;
    y(y>0)=y(y>0)*pfactor; % make protrusion less than indentation
    key = [x' y'];
end

function fvout = rotate_face_flat(fv)
    v1=fv(1,:);
    v2=fv(2,:);
    v3=fv(3,:);
    v4=fv(4,:);
    
    normal = cross(v2-v1, v3-v1); normal=normal/norm(normal);
    rdir = cross(normal,[0 0 1]);
    rdir=rdir/norm(rdir);
    if any(isnan(rdir))
        R = eye(3);
    else
        ra = acos(dot(normal,[0 0 1]));
        R = axang2rotm([rdir ra]);
    end
    
    fvout = fv*R';
    fvout = fvout - mean(fvout);
    fvout = fvout(:,[1 2]);
end

function keyorder = efdeg_to_keyorder(efdeg)
    assert(efdeg > 1);
    
    keyorder = [];
    for i=1:efdeg
        if i~=efdeg
            keyorder = [keyorder, i, i+1, i];
        else
            keyorder = [keyorder, i, 1, i];
        end
    end
    
end
