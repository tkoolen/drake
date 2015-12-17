classdef PolyController < DrakeSystem
  properties
    u_coeffs
    u_pows
    sos_plant
  end
  methods
    function obj = PolyController(plant,coeffs,pows)
      sos_plant = plant.sos_plant;
      obj = obj@DrakeSystem(0,0,sos_plant.num_states+3,sos_plant.num_inputs);
      obj.u_coeffs = coeffs;
      obj.u_pows = pows;
      obj = obj.setOutputFrame(plant.getInputFrame);
      obj = obj.setInputFrame(plant.getOutputFrame);
      obj.sos_plant = sos_plant;
    end
    
    function u = output(obj,t,~,x)
      t = t - x(end-2);
      x = x(1:end-3);
      u = obj.u_coeffs*prod(repmat([t;x]',length(obj.u_coeffs),1).^obj.u_pows,2)
      
      [umin, umax, A] = obj.sos_plant.simpleInputLimits(x);
%       u = max(u,umin);
%       u = min(u,umax);
%       u=1;
    end
  end
end