classdef UnivariateCubic < NStepCapturabilitySOSSystem
  methods
    function obj = UnivariateCubic()
      obj@NStepCapturabilitySOSSystem(1, 0, 0);
    end
  end
  
  methods
    function xdot = dynamics(obj, t, x, u)
      xdot = x * (x - 0.5) * (x + 0.5);
    end
    
    function xp = reset(obj, t, xm, s)
      xp = xm;
    end
    
    function ret = inputLimits(obj, u, x)
      ret = zeros(1, 1, 'like', u);
    end
    
    function [umin,umax,A] = simpleInputLimits(obj,x)
      umin = zeros(0, 1);
      umax = zeros(0, 1);
      A = [];
    end
    
    function ret = resetInputLimits(obj, s)
      ret = zeros(1, 1, 'like', s);
    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x, ~)
      % figure 1 in http://arxiv.org/pdf/1208.1751.pdf
      figure(1)
      xs = linspace(-1, 1, 100);
      Ws = msubs(Wsol,[x;t],[xs;zeros(size(xs))]);
%       Ws = zeros(length(xs),1);
%       for i = 1 : length(xs)
%         Ws(i) = subs(Wsol, [x; t], [xs(i); 0]);
%       end
      plot(xs, Ws);
      figure(2)
      Vs = msubs(Vsol,[x;t],[xs;zeros(size(xs))]);
      plot(xs, Vs);      
    end
    
    function[umin,umax,A] = simpleInputLimits(obj,x)
%       q = x(1:2);
%       z = q(2) + obj.z_nom;
      umin = [];
      umax = [];
      A = [];
    end
  end
end