classdef AtlasCMPController < MIMODrakeSystem
  % A simple jacobian-based standing controller that regulates the center of
  % mass while achieving end effector goals in the nullspace
  methods
    function obj = AtlasCMPController(r, dt)
      % @param sys the PD-controlled closed-loop system
      % @param r the time-stepping manipulator
      % @param k_com (optional) the COM objective gain
      % @param k_nom (optional) the nominal configuration gain
      
      typecheck(r,'TimeSteppingRigidBodyManipulator');
      
      input_frame = AtlasState(r);
      output_frame = AtlasInput(r);
      
      obj = obj@MIMODrakeSystem(0,output_frame.dim,input_frame,output_frame,false,true);
      obj = setSampleTime(obj,[dt; 0]); % sets controller update rate
      obj = setInputFrame(obj,input_frame);
      obj = setOutputFrame(obj,output_frame);
      obj.manip = r;
      obj.footIndices = find(~cellfun(@isempty,strfind(r.getLinkNames(),'foot')));
      
      obj.nq = obj.manip.getNumDOF();
      obj.nv = obj.manip.getNumDOF();
      
      obj.tau0 = zeros(output_frame.dim, 1);
      x0 = r.getInitialState();
      obj.q0 = x0(1 : obj.nq);
      obj.ankleJointIndices = find(~cellfun(@isempty,strfind(output_frame.coordinates,'leg_ak')));
    end
    
    function tau0 = getInitialState(obj)
      tau0 = obj.tau0;
    end
    
    function tau = mimoUpdate(obj,t,x,varargin)
      state = varargin{1};
      q = state(1 : obj.nq);
      v = state(obj.nq + 1 : end);
      kinsol = doKinematics(obj.manip, q, false, false);
      com = getCOM(obj.manip,kinsol);
      [A, Adot] = getCMM(obj.manip, kinsol, v);
      m = obj.manip.getMass();
      g = obj.manip.getGravity();
      Wg = [zeros(3, 1); m * g];
      
      % PD control
      J = eye(obj.nv);
      p = obj.kP * (obj.q0 - q) - obj.kD * v;

      % Contact point wrench matrix, Jacobians
      nFeet = length(obj.footIndices);
      nRho = 0;
      tauContactWrenches = zeros(obj.nv, 1);
      QFeet = cell(nFeet, 1);
      JFeetQFeet = cell(nFeet, 1);
      for i = 1 : nFeet
        footTransform = kinsol.T{obj.footIndices(i)};
        footTransform(1:3, 4) = footTransform(1:3, 4) - com;
        
        contactPoints = obj.manip.getBodyContacts(obj.footIndices(i));
        QFoot = computeContactWrenchMatrix(obj.mu, footTransform, contactPoints, obj.nBasisVectors);
        QFeet{i} = QFoot;
        nRho = nRho + size(QFeet{i}, 2);
        
        [JFoot, tauIndices] = geometricJacobian(...
          obj.manip.getManipulator(), kinsol, 1, obj.footIndices(i), obj.footIndices(i));
        JFoot = transformTwists(footTransform, JFoot); % express in CoM frame
        JFootFull = zeros(6, obj.nv);
        JFootFull(:, tauIndices) = JFoot;
        JFeetQFeet{i} = JFootFull' * QFoot;
      end
      rhoFootSize = nRho / nFeet;

      % Equations of motion
      [H,C,B] = manipulatorDynamics(obj.manip, q, v, false);
      
      % Quadratic program
      Wrho = obj.wRho * eye(nRho);
      W = blkdiag((J' * J), Wrho);
      f = [-2 * J * p; zeros(nRho, 1)];
      Aeq = [A, -horzcat(QFeet{:})];
      beq = Wg - Adot * v;
      lb = [-inf(obj.nv, 1); zeros(nRho, 1)];
      
      % append ankle torque stuff
%       S = zeros(length(obj.ankleJointIndices), size(B, 2));
%       for i = 1 : length(obj.ankleJointIndices)
%         S(i, obj.ankleJointIndices(i)) = 1;
%       end
%       BankleT = S * B';
%       Aeq = [Aeq; BankleT * H, -BankleT * horzcat(JFeetQFeet{:})];
%       beq = [beq; -BankleT * C];
      
      qp = QuadraticProgram(W,f,[],[],Aeq,beq,lb,[]);
      [xi, ~, exitflag] = solve(qp, zeros(obj.nv + nRho, 1), false);
      vd = xi(1 : obj.nv);
      rho = xi(obj.nv + 1 : end);

      columnStart = 1;
      for i = 1 : nFeet
        rhoFoot = rho(columnStart : columnStart + rhoFootSize - 1);
        tauContactWrenches = tauContactWrenches + JFeetQFeet{i} * rhoFoot;
        columnStart = columnStart + rhoFootSize;
      end
      
      % Inverse dynamics
      tau = B' * (H * vd + C - tauContactWrenches);
    end
    
    function y = mimoOutput(obj,t,x,varargin)
      y = x;
    end
    
  end
  
  properties
    manip;
    tau0;
    q0;
    nq;
    nv;
    kP = 100;
    kD = 20;
    mu = 0.5;
    epsilon = 1;
    nBasisVectors = 4;
    footIndices;
    wRho = 0.001;
    ankleJointIndices;
  end
end

function Q = computeContactWrenchMatrix(mu, transform, points, nVectorsPerPoint)
wrenchSize = 6;
nPoints = size(points, 2);
Q = zeros(wrenchSize, nVectorsPerPoint * nPoints);
startColumn = 1;
for j = 1 : size(points, 2)
  Qf = transform(1:3,1:3) * computeSupportVectors(mu, nVectorsPerPoint);
  point = transform(1:3,1:3) * points(:, j) + transform(1:3, 4);
  Qtau = vectorToSkewSymmetric(point) * Qf;
  QPoint = [Qtau; Qf];
  Q(:, startColumn : startColumn + nVectorsPerPoint - 1) = QPoint;
  startColumn = startColumn + nVectorsPerPoint;
end
end

function Beta = computeSupportVectors(mu, nVectors)
angleIncrement = 2 * pi / nVectors;
Beta = zeros(3, nVectors);
for i = 1 : nVectors
  angle = (i - 1) * angleIncrement;
  beta = [mu * cos(angle);
          mu * sin(angle);
          1];
  Beta(:, i) = beta / norm(beta);
end
end