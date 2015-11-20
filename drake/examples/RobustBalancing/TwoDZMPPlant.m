classdef TwoDZMPPlant < HybridDrakeSystem
  
  properties
    z_nom
    swingddot_lim = 1;
  end
  
  methods
    function obj = TwoDZMPPlant(z,mode_data)
      obj = obj@HybridDrakeSystem(0,6);
      obj.z_nom = z;
     
      final_plant = TwoDZMPStancePlant(z);
      final_plant = final_plant.setOutputFrame(obj.getOutputFrame);
      obj = obj.addMode(final_plant);
      
      for i=1:length(mode_data),
        stance_plant_i = TwoDZMPStancePlant(z,mode_data{i}.dVducoeffs,mode_data{i}.dVdupows);
        stance_plant_i = stance_plant_i.setOutputFrame(obj.getOutputFrame);
        obj = obj.addMode(stance_plant_i);
        
        guard_i = @(obj,t,x,u) guard(obj,t,x,u,mode_data{i}.Vpcoeffs,mode_data{i}.Vppows,mode_data{i}.dVpdxcoeffs,mode_data{i}.dVpdxpows);
        
        obj = obj.addTransition(1+i,guard_i,@transition,false,true);
      end                  
      
    end
    
    function phi = guard(obj,t,x,u,Vpcoeffs,Vppows,dVpdxcoeffs,dVpdxpows)
      xdot = obj.modes{1}.dynamics(t,x,u);  
      xdot = xdot([1;2;4;5]);
      x = x([1;2;4;5]);
      Vp = Vpcoeffs*prod(repmat(x',length(Vpcoeffs),1).^Vppows,2);
      
      % calculate Vpdot
      
      Vpdot = xdot'*dVpdxcoeffs*prod(repmat(x',length(dVpdxcoeffs),1).^dVpdxpows,2);
      
      phi = max(-Vp,Vpdot);
      display(sprintf('t: %f Vp: %f \t Vpdot: %f',t,Vp,Vpdot))
    end
    
    function [to_mode_xn,to_mode_num,status] = transition(obj,from_mode_num,t,x,u)
      to_mode_num = from_mode_num-1
      status = 0;
      to_mode_xn = [-x(2);-x(1);x(1)+x(2)+x(3);x(4);-x(4);0];
    end
  end
  
end

