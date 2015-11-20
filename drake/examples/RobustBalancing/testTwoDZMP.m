clear all
mode_1 = load('controller_V1');
mode_2 = load('controller_V2');

mode_data{1} = mode_1;
mode_data{2} = mode_2;
%%

z = 1;  
p = TwoDZMPPlant(z,mode_data);

%%

% x0 = [2;-.55;0;0;.7;0;0];
x0 = [2;-.1;0;0;.4;0;0];
% x0 = [3;0;0;0;.5;0;0];
traj = p.simulate([0 2],x0);

%%
v = TwoDZMPVisualizer(p);
v.playback(traj);