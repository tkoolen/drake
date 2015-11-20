classdef TwoDZMPStancePlant < SecondOrderSystem
  % states: [cm_x;swing_x;stance_x]
  properties
    z_nom
    coeffs
    pows
    swingddot_lim = 1;
  end
  
  methods
    function obj = TwoDZMPStancePlant(z,coeffs,pows)
      obj = obj@SecondOrderSystem(3,0,true);
      obj.z_nom = z;
      if nargin > 1
        obj.coeffs = coeffs;
        obj.pows = pows;
      end
    end
    
    function qdd = sodynamics(obj,t,q,qd,u)
      q = q([1;2]);
      qd = qd([1;2]);
      if ~isempty(obj.coeffs)
        dVdswingddot = obj.coeffs*prod(repmat([q;qd]',length(obj.coeffs),1).^obj.pows,2);
        swingddot = obj.swingddot_lim*sign(dVdswingddot);
      else
        swingddot = 0;
      end
      qdd = [q(1)/obj.z_nom;swingddot;0];      
    end
  end
  
end

