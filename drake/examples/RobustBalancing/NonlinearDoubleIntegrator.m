classdef NonlinearDoubleIntegrator < NStepCapturabilitySOSSystem
  % From Korda et al 2014
  
  properties
  end
  
  methods
    function obj = NonlinearDoubleIntegrator(umax)
      obj@NStepCapturabilitySOSSystem(2, 1, 0);
    end
    
    function xdot = dynamics(obj, t, x, u)
      xdot = [x(2) + .1*x(1)^3;.3*u]; 
    end
    
    function xp = reset(obj, t, xm, s)
      % control input changes q only
      % qp = qm - u
      xp = xm;
    end
    
    function ret = inputLimits(obj, u, x)
      ret = 1-u^2;
    end    
        
    function[umin,umax,A] = simpleInputLimits(obj,x)
      umin = -1;
      umax = 1;
%       umax = [];
%       A = eye(2)/obj.cop_max^2;
A = [];
    end
    
    % Transform u
    % y = A*u + b s.t. limits on y are |y_i| <= 1
    % u = C*y + d inverse transform
    function [A,b,C,d] = unitBoxInputTransform(obj)
      A = eye(1);
      b = zeros(1,1);
      C = eye(1);
      d = zeros(1,1);
    end
    
    
    function ret = resetInputLimits(obj, s)
      ret = [];
    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x, u_sol)      
      sub_vars = [t];
      sub_val = [0];
      plot_vars = [x];
      
      
      figure(1)
      contourSpotless([Wsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[1 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('W(x)')
      
      figure(n*10+2)
      contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[0 0],{'b','r','g'});
      xlabel('q_1')
      ylabel('v_1')
      title('V(0,x)')
      
      
      figure(n*10+4)
      hold off
      h=contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[0 0],{'k','r'});
      hold on
      h=contourSpotless(u_sol,plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val);
      xlabel('q_1')
      ylabel('v_1')
      title('u(0,x)')
    end
    
    function draw(obj,t,x)
      x_stance = x(end-2:end-1);
      % draw line from origin to COM
      h=line([x_stance(1);x(1)],[x_stance(2);x(2)]);
      set(h,'LineWidth',3,'Color','red')
      radius = .1;
      rectangle('Position',[x(1:2)+x_stance-radius/2;radius;radius],'Curvature',[1,1], 'FaceColor','k')
      xlim([-2 2])
      ylim([-2 2])
      axis off
    end
    
  end
end
