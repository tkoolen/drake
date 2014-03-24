function A = getCMM2(robot, kinsol)

com = robot.getCOM(kinsol.q);
nv = robot.getNumStates() / 2;
nBodies = length(robot.body);
momentumSize = 6;

transformedMotionSubspaces = zeros(momentumSize, nv);
A = zeros(momentumSize, nv);
ancestorMatrix = false(nBodies, nBodies);
for i = 2 : nBodies % don't process 'world' body
  body = robot.body(i);
  I = body.I;
  HBodyToCoM = kinsol.T{i};
  HBodyToCoM(1:3, 4) = HBodyToCoM(1:3, 4) - com;
  HCoMToBody = transformInverse(HBodyToCoM);
  dofnum = body.dofnum;
  X = motionSubspace(body, kinsol.q(dofnum));
  transformedMotionSubspaces(:, dofnum) = transformTwists(HBodyToCoM, X);
  
  ancestorMatrix(i, :) = ancestorMatrix(body.parent, :);
  ancestorMatrix(i, i) = true;
  ancestors = ancestorMatrix(i, :); % including self, not including world
  
  dofs = vertcat(robot.body(ancestors).dofnum);
  J = transformedMotionSubspaces(:, dofs);
  AdH = transformAdjoint(HCoMToBody);
  A(:, dofs) = A(:, dofs) + AdH' * I * AdH * J;
end

end

function HInv = transformInverse(H)
  HInv = zeros(4, 4);
  R = H(1:3, 1:3);
  p = H(1:3, 4);
  HInv(1:3, 1:3) = R';
  HInv(1:3, 4) = -R' * p;
  HInv(4, 4) = 1;
end
