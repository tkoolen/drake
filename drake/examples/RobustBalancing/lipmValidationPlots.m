N = 6;
for n=1:N,
  data{n} = load(sprintf('V%d_LIPM3D',n-1));
end

x = msspoly('x',4);
t = msspoly('t',1);

q = x(1 : 2);
v = x(3 : 4);

sub_vars = [q(2);v(2);t];
sub_val = [0;0;0];
plot_vars = [q(1);v(1)];
figure(1)
clf
hold on
for n=N:-1:1,
  Vsol = data{n}.Vsol;
  obj = data{n}.model;
  R_diag = data{n}.R_diag;
  A = diag(1./(R_diag.^2))
  h_X = 1 - x'*A*x;
  
  % from Koolen et. al IJRR
  % regions should depend on the instantaneous capture point
  r_ic = q + v*sqrt(obj.z_nom / obj.g);
  dN = lipmCaptureLimit(obj.T, obj.cop_max, obj.step_max, obj.z_nom, obj.g, n); % theoretical max ICP distance
%   contourSpotless([Vsol;h_X;r_ic'*r_ic],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val,[0 0 dN^2],{'b','r','g'});
  h=contourSpotless([Vsol],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val,[0],{'b','r','g'},100,100,n/2);
  set(h,'Fill','On')
  h=contourSpotless([h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val,[0],{'k'});
  set(h,'LineWidth',5)
end
%%

set(gca,'LooseInset',get(gca,'TightInset'))

  xlabel('q_1','FontSize',24)
  ylabel('v_1','FontSize',24)
  title('Capture Regions','FontSize',24)
hold off