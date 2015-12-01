function lipm3DPlotFun(n, Vsol, Wsol, h_X, R_diag, t, x, model_params)
T = model_params.T;
g = model_params.g;
z_nom = model_params.z_nom;
step_max = model_params.step_max;

q = x(1 : 2);
v = x(3 : 4);

sub_vars = [q(2);v(2);t];
sub_val = [0;0;0];
plot_vars = [q(1);v(1)];

figure(1)
contourSpotless([Wsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[1 0],{'b','r'});
xlabel('q_1')
ylabel('v_1')
title('W(x)')

% from Koolen et. al IJRR
% regions should depend on the instantaneous capture point
r_ic = q + v*sqrt(z_nom/g);
dN = captureLimit(T, 0, step_max, z_nom, g, n); % theoretical max ICP distance

figure(n*10+2)
contourSpotless([Vsol;h_X;r_ic'*r_ic],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[0 0 dN^2],{'b','r','g'});
xlabel('q_1')
ylabel('v_1')
title('V(0,x)')
end
