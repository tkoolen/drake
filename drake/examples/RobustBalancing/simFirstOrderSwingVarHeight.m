clear all
data_0 = load('V0_FirstOrderSwingVariableHeight');
data_1 = load('V1_FirstOrderSwingVariableHeight');
data = {data_0;data_1};
% data = {data_1};
model = data_0.model;
p = HybridCapturabilityPlant(model,data,true);
t = msspoly('t',1);
x = msspoly('x',model.num_states);
%%
% x0 = [2;.2;0;0;.9;0];
% x0 = [1;-.15;0;0;.6;0];
x0 = [2;0;0;0;.6;0];
x0 = [2;.3;0;0;.4;0];
x0 = [2;.1;0;-.1;0;0];
x0 = [2;0;0;0;.8;0];
x0 = [2;0;0;0;.9;0];
% x0 = [2;.3;0;0;.4;0];
traj = p.simulate([0 1.3],[x0;0;0;0]);
plant = NStepCapturabilityPlant(model);
v = NStepCapturabilityVisualizer(plant);
v = v.setInputFrame(p.getOutputFrame);
v.playback_speed = .5;

v.playback(traj);

%%
t_sim = traj.pp.breaks;
x_sim = traj.eval(t_sim);
V_sim = msubs(data_1.Vsol,[t;x],[t_sim;x_sim(1:end-3,:)]);
% plot(t_sim,V_sim);
%%
% sample sim
