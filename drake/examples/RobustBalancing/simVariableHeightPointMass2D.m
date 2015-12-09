function simVariableHeightPointMass2D
load V1_VariableHeightPointMass2D
V0 = load('V0_VariableHeightPointMass2D');

g = 10;
z_nom = 1;
step_max = .7;
step_time = 0.3;
cop_max = 0.0; % set to 0 to get point foot model with no continuous inputs

f_max = 1.1;
f_min = .9;

model = VariableHeightPointMass2D(g, z_nom, step_max, step_time, f_max, f_min);

t = msspoly('t',1);
x = msspoly('x',model.num_states);
u = msspoly('u',model.num_inputs);
s = msspoly('s',model.num_reset_inputs);
xdot = model.dynamics(t,x,u);

Vdot = diff(Vsol,x)*xdot + diff(Vsol,t);
dVdotdu = diff(Vdot,u);

[pows,coeffs] = decomp_ordered(dVdotdu,[t;x]);

plant = NStepCapturabilityPlant(model);
controller = NStepCapturabilityController(plant,coeffs,pows);


sys_cl = plant.feedback(controller);


xp = model.reset(t, x, s);

V0p = subs(V0.Vsol,[x;t],[xp;0]);

%%
x0 = [0;0;.9;0];
sys_cl = sys_cl.setSimulinkParam('Solver','ode4');
sys_cl = sys_cl.setSimulinkParam('FixedStep','.001');
traj = sys_cl.simulate([0 step_time],x0);


t_sim = traj.pp.breaks;
x_sim = traj.eval(t_sim);
V_sim = msubs(Vsol,[t;x],[t_sim;x_sim]);
figure(1)
subplot(2,1,1)
plot(t_sim,x_sim)
legend('x_1','x_2','x_3','x_4')

%%
V0p_sim = subs(V0p,x,x_sim(:,end));
S = linspace(-step_max,step_max,1000);
V0p_vals = msubs(V0p_sim,s,S);
[V0_opt,i_opt] = max(V0p_vals);
V0_opt
s_opt = S(i_opt)

%% 
x0_0step = double(subs(xp,[x;s],[x_sim(:,end);s_opt]));
% x0_0step = [0-.5;0;1.623;0];
Vdot_0step = diff(V0.Vsol,x)*xdot + diff(V0.Vsol,t);
dVdotdu_0step = diff(Vdot_0step,u);

[pows_0step,coeffs_0step] = decomp_ordered(dVdotdu_0step,[t;x]);
controller_0step =  NStepCapturabilityController(plant,coeffs_0step,pows_0step);

sys_cl_0step = plant.feedback(controller_0step);
sys_cl_0step = sys_cl_0step.setSimulinkParam('Solver','ode4');
sys_cl_0step = sys_cl_0step.setSimulinkParam('FixedStep','.001');
traj_0step = sys_cl_0step.simulate([0 1],x0_0step);

t_sim_0step = traj_0step.pp.breaks;
x_sim_0step = traj_0step.eval(t_sim_0step);

V0_sim_0step = msubs(V0.Vsol,[t;x],[t_sim_0step;x_sim_0step]);

figure(1)
subplot(2,1,2)
plot(t_sim_0step,x_sim_0step)
legend('x_1','x_2','x_3','x_4')
end