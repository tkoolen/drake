classdef PolyController < DrakeSystem
  properties
    u_coeffs_num
    u_pows_num
    u_coeffs_den
    u_pows_den
    sos_plant
    nBaseStates
  end
  methods
    function obj = PolyController(plant,t_mss,x_mss,u_num,base_states,u_den)
      if nargin < 5
        base_states = true;
      end
      if nargin < 6
        u_den = msspoly(1);
      end
      
      if base_states
        nBaseStates = 3;
      else
        nBaseStates = 0;
      end
      
      sos_plant = plant.sos_plant;
      obj = obj@DrakeSystem(0,0,sos_plant.num_states+nBaseStates,sos_plant.num_inputs);
      
      obj.nBaseStates = nBaseStates;
      
      [pows_num,coeffs_num] = decomp_ordered(u_num,[t_mss;x_mss]);
      [pows_den,coeffs_den] = decomp_ordered(u_den,[t_mss;x_mss]);

      obj.u_coeffs_num = coeffs_num;
      obj.u_pows_num = pows_num;
      obj.u_coeffs_den = coeffs_den;
      obj.u_pows_den = pows_den;
      obj = obj.setOutputFrame(plant.getInputFrame);
      obj = obj.setInputFrame(plant.getOutputFrame);
      obj.sos_plant = sos_plant;
    end
    
    function u = output(obj,t,~,x)
      if obj.nBaseStates > 0
        t = t - x(end-2);
        x = x(1:end-3);
      end
      num = obj.u_coeffs_num*prod(repmat([t;x]',length(obj.u_coeffs_num),1).^obj.u_pows_num,2);      
      den = obj.u_coeffs_den*prod(repmat([t;x]',length(obj.u_coeffs_den),1).^obj.u_pows_den,2);
      u = num/den;
      [umin, umax, A] = obj.sos_plant.simpleInputLimits(x);
      
      if u > umax
        display('u>umax');
        u
      elseif u < umin
        display('u<umin')
        u
      end
      
      u = max(u,umin);
      u = min(u,umax);
    end
  end
end