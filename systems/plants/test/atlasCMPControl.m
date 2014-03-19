function atlasCMPControl()
robot = createAtlas('rpy');
n = robot.getNumStates();
% x0 = zeros(n, 1);

% TODO: cascade a Trajectory object describing u(t) with the system, then simulate.
traj = simulate(robot, [0 10]);
fnplt(traj(1));
end