load V1_LIPM3D
V0 = load('V0_LIPM3D');

g = 10;
z_nom = 1;
step_max = .7;
step_time = 0.3;
cop_max = 0.1; % set to 0 to get point foot model with no continuous inputs

model = LIPM3D(g, z_nom, step_max, step_time, cop_max);

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
x0 = [-.15;0;1.61;0];
% x0 = [-.2;0;1.67;0];
sys_cl = sys_cl.setSimulinkParam('Solver','ode4');
sys_cl = sys_cl.setSimulinkParam('FixedStep','.001');
traj = sys_cl.simulate([0 step_time],x0);


t_sim = traj.pp.breaks;
x_sim = traj.eval(t_sim);
V_sim = msubs(Vsol,[t;x],[t_sim;x_sim]);
figure(1)
subplot(2,1,1)
plot(t_sim,x_sim)
subplot(2,1,2)
plot(t_sim,V_sim)

%%
V0p_sim = subs(V0p,x,x_sim(:,end));
[R,THETA] = meshgrid(linspace(0,step_max,301),linspace(-pi,pi,301));
S1 = R(:).*cos(THETA(:));
S2 = R(:).*sin(THETA(:));
V0p_vals = msubs(V0p_sim,s,[S1';S2']);
[V0_opt,i_opt] = max(V0p_vals);
V0_opt
s_opt = [S1(i_opt);S2(i_opt)];