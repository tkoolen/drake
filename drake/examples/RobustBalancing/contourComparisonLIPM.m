variableHeight = load('V1_VariableHeightPointMass2D');
variableHeightAndPitch = load('V1_VariableHeightandPitch2D');
lipm = load('V1_LIPM3D');

x = msspoly('x',6);
t = msspoly('t',1);

figure(1)
clf
hold on

h=contourSpotless([lipm.Vsol],x(1),x(3), [-2 2],[-2 2],[t;x([2;4])],zeros(3,1),0,{'g'});
h=contourSpotless([variableHeight.Vsol],x(1),x(3), [-2 2],[-2 2],[t;x([2;4;5;6])],zeros(5,1),0,{'b'});
h=contourSpotless([variableHeightAndPitch.Vsol],x(1),x(4), [-2 2],[-2 2],[t;x([2;3;5;6])],zeros(5,1),0,{'k'});
% plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val,[0],{'b','r','g'},100,100,n/2);
