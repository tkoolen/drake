% clear all
% load V1_LIPM3D
V0 = load('V0_LIPM2D');
data = {V0};
model = V0.model;
% data{1}.u_sol = usol;
% p = HybridCapturabilityPlant(model,data);

%%
t = msspoly('t',1);
x = msspoly('x',model.num_states);

plant = NStepCapturabilityPlant(model);
controller = PolyController(plant,t,x,num,den);
p = feedback(plant,controller);

% x0 = [1;-.1;.5];
% x0 = [1;-.29;.31];
x0 = [0;.1];
% x0 = [1;-0;.1];
% x0 = [1;0;0];
% x0 = [1;-.2;0;.75;0];
% x0(2) = -x0(4)/sqrt(model.g);
T = 1;

traj = p.simulate([0 T],[x0;0;0;0]);
plant = NStepCapturabilityPlant(model);
v = NStepCapturabilityVisualizer(plant);
v = v.setInputFrame(p.getOutputFrame);
v.playback_speed = T;
figure(25)
v.playback(traj);

%% 6

t_sim = traj.pp.breaks;
x_sim = traj.eval(t_sim);
V0 = msubs(V0.Vsol,[t;x],[t_sim;x_sim(1:model.num_states,:)]);
figure(5)
subplot(3,1,1)
plot(t_sim,x_sim(1:model.num_states,:))
subplot(3,1,2)
plot(t_sim,V0)
grid on
subplot(3,1,3)
plot(t_sim,sqrt(sum(x_sim(1:2,:).*x_sim(1:2,:),1)))
grid on

min(sqrt(sum(x_sim(1:2,:).*x_sim(1:2,:),1)))

norm(min(sqrt(sum(x_sim(1:2,:).*x_sim(1:2,:),1))))