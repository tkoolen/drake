function testCentroidalMomentumMatrix2()
testAtlas('rpy');
end

function testAtlas(floatingJointType)
robot = createAtlas(floatingJointType);

nq = robot.getNumStates() / 2; % TODO
nv = robot.getNumStates() / 2; % TODO

t1 = 0; t2 = 0;

nTests = 100;
for i = 1 : nTests
  q = randn(nq, 1);
  v = randn(nv, 1);
  kinsol = robot.doKinematics(q,false,false, v);

  tic;
  A1 = robot.getCMM(kinsol);
  t1 = t1 + toc;

  tic;
  A2 = robot.getCMM2(kinsol);
  t2 = t2 + toc;

  valuecheck(A1, A2, 1e-10);
end

displayComputationTime = false;
if displayComputationTime
  fprintf('Old method time: %0.5f\n', t1/nTests);
  fprintf('New method time: %0.5f\n', t2/nTests);
end

end