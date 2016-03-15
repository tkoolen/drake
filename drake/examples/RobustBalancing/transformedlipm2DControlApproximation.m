clear all
data = load('V0_TransformedLIPM2D');
model = data.model;
R_diag = data.R_diag;
x=msspoly('x', model.num_states); % state
u=msspoly('u', model.num_inputs);
t=msspoly('t', 1); 

xdot = model.dynamics(t,x,u);
g = diff(xdot,u);

options.T = data.T;
options.time_varying = true;
options.degree = 6;
usol = msspoly();
for i=1:model.num_inputs,
  num = g(:,i)'*diff(data.Vsol,x)';
%   [p,q] = normApproximation(vec'*vec,R_diag,options);
%   usol = [usol;g(:,i)'*diff(data.Vsol,x)'*p];
  den = normUpperBound(num,R_diag,options);
end

%%
goal_radius = .2;
target = @(x) goal_radius^2 - x'*x;
options.target = target;
options.degree = 6;
options.u_den = den;
% [Vsol,Wsol] = finiteTimeInnerApproximation(model,options.T,num,R_diag,target,options);
[Vsol,qsol] = lyapunovInnerApproximation(model,data.Vsol,options.T,num,R_diag,options);

A = diag(1./(R_diag.^2));
h_X = 1 - x'*A*x;
figure(1)
hold off
contourSpotless(data.Vsol,x(1),x(2),[-1 1], [-1 1],t,0,0)
hold on
contourSpotless(Vsol,x(1),x(2),[-1 1], [-1 1],t,0,1)
contourSpotless(h_X,x(1),x(2),[-1 1], [-1 1],t,0,0)
hold off

%%
num2= -g(:,i)'*diff(Vsol,x)';
den2 = normUpperBound(num2,R_diag,options);

%%
options2 = options;
options2.u_den = den2;
[Vsol2,qsol2] = lyapunovInnerApproximation(model,data.Vsol,options.T,num2,R_diag,options2);

%%
A = diag(1./(R_diag.^2));
h_X = 1 - x'*A*x;
figure(2)
hold off
contourSpotless(data.Vsol,x(1),x(2),[-1 1], [-1 1],t,0,0)
hold on
contourSpotless(Vsol2,x(1),x(2),[-1 1], [-1 1],t,0,1)
contourSpotless(h_X,x(1),x(2),[-1 1], [-1 1],t,0,0)
hold off