function runAtlasCMPControl()
options.dt = 0.001;
r = Atlas('urdf/atlas_minimal_contact.urdf',options);

load('data/atlas_fp.mat');
r = r.setInitialState(xstar);
c = AtlasCMPController(r, options.dt);

sys = FeedbackSystem(r, c);

T = 5; % sec
S=warning('off','Drake:DrakeSystem:UnsupportedSampleTime');
% output_select(1).system=1;
% output_select(1).output=1;

% sys = mimoCascade(sys,v,[],[],output_select);
warning(S);
traj = simulate(sys,[0 T]);

visualize = true;
if visualize
  v = r.constructVisualizer();
  v.display_dt = .05;
  sys = cascade(sys,v);
  simulate(sys, [0 T]);
  playback(v,traj,struct('slider',true));
end
end