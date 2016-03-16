g = 10;
z_nom = 1;
step_max = .7;
step_time = 0.3;
cop_max = .1; % set to 0 to get point foot model with no continuous inputs
fz_range = .3;
fx_range = .3;
inertia_ratio = .6^2/2;

model = TransformedFull2DModel(g, inertia_ratio, z_nom, step_max, step_time, fz_range, fx_range, cop_max);


%% Get an initial quadratic Lyapunov candidate
x = msspoly('x',model.num_states);
t = msspoly('t',1);
u = msspoly('u',model.num_inputs);
f = model.dynamics(t,x,u);

A = double(subs(diff(f,x),[t;x;u],zeros(1+model.num_states+model.num_inputs,1)));
B = double(subs(diff(f,u),[t;x;u],zeros(1+model.num_states+model.num_inputs,1)));

Q = 100*eye(model.num_states);
R = eye(model.num_inputs);
[K,S] = lqr(A,B,Q,R);

%%
V0 = x'*S*x;
[V,u_fn] = quadraticControlLyapunovAlternations(x,u,f,V0*100)
%%
for i=1:15,
  [V,u_fn] = quadraticControlLyapunovAlternations(x,u,f,V);
  contourSpotless(V,x(1),x(2),[-3 3],[-1 1],[t;x(3:end)],zeros(model.num_states-1,1),1);
end;

%%
s = msspoly('s',model.num_reset_inputs);
r = model.reset(t,x,s);

% get inner approximation of pre-reset
% s = (A'SA)^-1 * A'Sx
S=double(diff(diff(V,x)',x))/2;
T=eye(model.num_states) + A*inv(A'*S*A)*A'*S;
T=T'*Q*T;

% S < T
% s = (A'SA)^-1 * A'Sx > 1 ==> (x+A)'S(x+A)

% Vm(x) <= 1 ==> x'Tx <= 1

%  quadraticControlAlternationsWithReset(x,u,s,f,r,V,step_time,V)
