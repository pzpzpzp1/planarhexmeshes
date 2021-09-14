% projects faces with vertices V: [nF 3 4]
% to be planar. projections are decoupled.
function [Vproj, Vin] = planar_face_projections(Vin, visualize)
    % clc; close all;
    if nargin ==0
        visualize = 1;
        nF = 20;
        Vin = randn(nF,3,4);
        Vin(:,1,:) = Vin(:,1,:) + [1:nF]'*2; % space quads apart by 2.
    end
    nF = size(Vin,1);

    options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'HessianApproximation','lbfgs',...
        'StepTolerance',1e-10,'ConstraintTolerance',1e-5,'Display','iter','OptimalityTolerance',1e-6 ,...
        'CheckGradients',false,'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-3);

    fun = @(x) objfun(Vin,x);
    [Vproj,fval,exitflag,output,lambda,~,~] = fmincon(fun,Vin(:),[],[],[],[],[],[],@mycon,options);
    Vproj = reshape(Vproj,[],3,4);
    
    %% visualize
    if visualize
        figure; axis equal off; hold all; rotate3d on;
        for i=1:nF
            Vi  = permute(Vin(i,:,:),[3 2 1]);
            patch('vertices',Vi,'faces',[1 2 3 4],'facecolor','red')
            Viproj  = permute(Vproj(i,:,:),[3 2 1]);
            patch('vertices',Viproj,'faces',[1 2 3 4],'facecolor','green')
        end
    end
    
    %% check optimality
    if nargin == 0
        [planarity, pgrad] = face_planarity(Vproj);
        diff = Vin - Vproj;
        dualvars = (reshape(diff(:)./pgrad(:),nF, []));
        fmincondualvars = lambda.eqnonlin;
        [mean(dualvars,2) fmincondualvars std(dualvars')']
        display('dual variables and stdev of component dual vars. 3rd col should be 0. first 2 should be identical.')
    end
    
end

function [val, grad] = objfun(V,x)
    grad = -(V(:)-x(:));
    val = .5 * grad' * grad;
end

function [c, ceq, gradc, gradceq] = mycon(Vin)
    c=0;
    Vin = reshape(Vin,[],3,4);
    nF = size(Vin,1);
    
    [planarity, pgrad] = face_planarity(Vin);
    
    ceq = planarity';
    gradceq = sparse(1:numel(Vin), repmat([1:nF]',12,1), pgrad(:), numel(Vin), nF);
    
    gradc = zeros(numel(Vin),1);
end





