function testDoKinematicsComputationTime()
robot = createAtlas('rpy');

nq = robot.getNumStates() / 2; % TODO
nv = robot.getNumStates() / 2; % TODO

nTests = 1000;

time = 0;
for i = 1 : nTests
  if mod(i, 100) == 0
    fprintf('test %d\n', i);
  end
  
  q = randn(nq, 1);
  v = randn(nv, 1);
  
  tic
  kinsol = robot.doKinematics(q,false,false, v);
  time = time + toc * 1e3;
end
time = time / nTests;

fprintf('took %0.3f ms per call\n', time);

end


function robot = createAtlas(floatingJointType)
options.floating = floatingJointType;
w = warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits');
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
robot = RigidBodyManipulator(fullfile('../../../examples/Atlas/urdf/atlas_minimal_contact.urdf'),options);
warning(w);
end
