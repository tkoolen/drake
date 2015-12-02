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
    
    function ret = resetInputLimits(obj, s)
      ret = zeros(1, 1, 'like', s);
    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x)
      % figure 1 in http://arxiv.org/pdf/1208.1751.pdf
      xs = linspace(-1, 1, 100);
      Ws = zeros(length(xs));
      for i = 1 : length(xs)
        Ws(i) = subs(Wsol, [x; t], [xs(i); 0]);
      end
      plot(xs, Ws);
    end
  end
end