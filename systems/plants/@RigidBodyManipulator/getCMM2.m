function A = getCMM2(robot, kinsol)

com = robot.getCOM(kinsol.q);
nv = robot.getNumStates() / 2;
nBodies = length(robot.body);
momentumSize = 6;

transformedMotionSubspaces = zeros(momentumSize, nv);
A = zeros(momentumSize, nv);
for i = 2 : nBodies % don't process 'world' body
  body = robot.body(i);
  I = body.I;
  HBodyToCoM = kinsol.T{i};
  HBodyToCoM(1:3, 4) = HBodyToCoM(1:3, 4) - com;
  HCoMToBody = inv(HBodyToCoM);
  dofnum = body.dofnum;
  X = motionSubspace(body, kinsol.q(dofnum));
  transformedMotionSubspaces(:, dofnum) = transformTwists(HBodyToCoM, X);
  
  ancestors = findAncestorBodies(robot, i);
  % prepend current body, get rid of world body:
  chainToWorld = [i; ancestors(1 : end - 1)];
  
  J = zeros(6, nv);
  dofs = vertcat(robot.body(chainToWorld).dofnum);
  J(:, dofs) = transformedMotionSubspaces(:, dofs);
  AdH = transformAdjoint(HCoMToBody);
  A = A + AdH' * I * AdH * J;
end

end
