data_0 = load('V0_VariableHeightPointMass2D');
data_1 = load('V1_VariableHeightPointMass2D');
data_2 = load('V2_VariableHeightPointMass2D');
data = {data_0;data_1;data_2};
% data = {data_1};
model = data_0.model;
p = HybridCapturabilityPlant(model,data);
%%
x0 = [2;-.26;0;1.8;0]; 
% x0 = [2;0;0;0;0];
% x0 = [3;-.9;0;1.4;0];
% x0 = [3;0;0;1;0];
traj = p.simulate([0 .3*(x0(1)-1)+1],[x0;0;0;0]);
plant = NStepCapturabilityPlant(model);
v = NStepCapturabilityVisualizer(plant);
v = v.setInputFrame(p.getOutputFrame);
v.playback_speed = .5;

v.playback(traj);
