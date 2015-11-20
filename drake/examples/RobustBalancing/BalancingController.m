classdef BalancingController < DrakeSystem
% File remains a work in progress
  properties
    plant
    coeffs
    pows
    z_nom
    % what is pitch-nominal for the feedback?
  end
  
  methods
    function obj = BalancingController(p,coeffs,pows)
      nx = p.getNumStates;
      nu = p.getNumInputs;
      obj = obj@DrakeSystem(0,0,nx,nu,true,false);
      obj.plant = p;
      obj = obj.setOutputFrame(p.getInputFrame);
      obj = obj.setInputFrame(p.getStateFrame);
      obj.coeffs = coeffs;
      obj.pows = pows;
    end
    
    function u = output(obj,t,~,x)
      zddot = 0;
      angddot = 0;
      swingddot = [0;0];
      
      nq = obj.plant.getNumPositions;
      q = x(1:nq);
      qd = x(nq+1:end);
      [com,Jcom] = obj.plant.getCOM(q);
      comdot = Jcom*qd;      
      kinsol = obj.plant.doKinematics(q,qd);
      
      % crude guess at stance leg
      [phi,normal,d,xA,xB,idxA,idxB,mu,n,D] = obj.plant.contactConstraints(kinsol);
      phidot = n*qd;
      
      if phi(1)^2 + phidot(1)^2 < phi(2)^2 + phidot(2)^2
        [r_swing,J_swing] = obj.plant.forwardKin(kinsol,idxB(2),xB(:,2));
        r_stance = obj.plant.forwardKin(kinsol,idxB(2),xB(:,2));
      else
        [r_swing,J_swing] = obj.plant.forwardKin(kinsol,idxB(1),xB(:,1));
        r_stance = obj.plant.forwardKin(kinsol,idxB(1),xB(:,1));
      end
      rdot_swing = J_swing*qd;
      
      z = [com(1) - r_stance(1);comdot(1);r_swing(1);rdot_swing(1)];
      dVdswingddot = obj.coeffs*prod(repmat(z',length(obj.coeffs),1).^obj.pows,2);
      swingddot(1) = sign(dVdswingddot);
      
      
      % todo--choose a reasonable swingddot(2)
      %     -add feedback to zddot and angddot
      %     -code this up better 
      %     -get real initial conditions
      
      u = obj.controlFunction(t,[],x,zddot,angddot,swingddot);
    end
    
    function u = controlFunction(obj,t,~,x,zddot,angddot,swingddot)
      nq = obj.plant.getNumPositions;
      q = x(1:nq);
      qd = x(nq+1:end);
      nu = obj.plant.getNumInputs;
      nf = 2;
      
      kinsol = obj.plant.doKinematics(q,qd);
      
      % get momentum matrix and convert to planar
      % [angmom;xmom;zmom]
      A = obj.plant.centroidalMomentumMatrix(kinsol);
      A = A([2;4;6],:);
      Adotv = obj.plant.centroidalMomentumMatrixDotTimesV(kinsol);
      Adotv = Adotv([2;4;6]);
      
      % crude guess at stance leg
      [phi,normal,d,xA,xB,idxA,idxB,mu,n,D] = obj.plant.contactConstraints(kinsol);
      phidot = n*qd;
      
      if phi(1)^2 + phidot(1)^2 < phi(2)^2 + phidot(2)^2
        [~,J_swing] = obj.plant.forwardKin(kinsol,idxB(2),xB(:,2));
        Jdotv_swing = obj.plant.forwardJacDotTimesV(kinsol,idxB(2),xB(:,2));
        
        [~,J_stance] = obj.plant.forwardKin(kinsol,idxB(1),xB(:,1));
        Jdotv_stance = obj.plant.forwardJacDotTimesV(kinsol,idxB(1),xB(:,1));
        
        J_f = [n(1,:)+mu(1)*D{1}(1,:);n(1,:)+mu(1)*D{2}(1,:)];
      else
        [~,J_swing] = obj.plant.forwardKin(kinsol,idxB(1),xB(:,1));
        Jdotv_swing = obj.plant.forwardJacDotTimesV(kinsol,idxB(1),xB(:,1));
        
        [~,J_stance] = obj.plant.forwardKin(kinsol,idxB(2),xB(:,2));
        Jdotv_stance = obj.plant.forwardJacDotTimesV(kinsol,idxB(2),xB(:,2));
        
        J_f = [n(2,:)+mu(2)*D{1}(2,:);n(1,:)+mu(2)*D{2}(2,:)];
      end
      
      swingddot = swingddot - 100*J_swing([1;3],:)*qd;
      
      M = [A([1;3],:); J_swing([1;3],:); J_stance([1;3],:)];
      w = [angddot;zddot;swingddot;0;0] - [Adotv([1;3]);Jdotv_swing([1;3],:);Jdotv_stance([1;3],:)];
      
      [H,C,B] = obj.plant.manipulatorDynamics(q,qd);
      Hinv = inv(H);
      
      
      Q = diag([ones(nu,1);1e-3*ones(nf,1)]);
      h = zeros(nu+nf,1);
      
      A = [M*Hinv*B M*Hinv*J_f'];
      b = w + M*Hinv*C;
      [z,fval,flag]=quadprog(Q,h,[],[],A,b,[obj.plant.umin;zeros(nf,1)],[obj.plant.umax;inf(nf,1)]);
      u = z(1:nu);
      f = z(nu+1:nu+nf);
      
      if flag ~= 1
        keyboard
      end
    end
  end
  
end

