classdef HybridCapturabilityPlant < HybridDrakeSystem
  % augmented state with [t_offset;x_stance;y_stance]
  properties
    sos_plant
    data
    t
    x
    u
    s
    xp
  end
  
  methods
    function obj = HybridCapturabilityPlant(sos_plant,data)
      obj = obj@HybridDrakeSystem(sos_plant.num_inputs,sos_plant.num_states+3);
      obj.sos_plant = sos_plant;
      
      N = length(data);
      
      ol_plant = NStepCapturabilityPlant(sos_plant);
      obj = obj.setOutputFrame(ol_plant.getOutputFrame);
      obj = obj.setInputFrame(ol_plant.getInputFrame);
      
      t = msspoly('t',1);
      x = msspoly('x',sos_plant.num_states);
      u = msspoly('u',sos_plant.num_inputs);
      s = msspoly('s',sos_plant.num_reset_inputs);
      
      xdot = sos_plant.dynamics(t,x,u);
      xp = sos_plant.reset(t, x, s);
      
      for i=1:N,
%         ol_plant = NStepCapturabilityPlant(sos_plant);
%         ol_plant = ol_plant.setOutputFrame(ol_plant.getOutputFrame);
%         ol_plant = ol_plant.setInputFrame(ol_plant.getInputFrame);
        
        Vdot = diff(data{i}.Vsol,x)*xdot + diff(data{i}.Vsol,t);
        dVdotdu = diff(Vdot,u);
        [pows,coeffs] = decomp_ordered(dVdotdu,[t;x]);
        controller = NStepCapturabilityController(ol_plant,coeffs,pows);
        plant = ol_plant.feedback(controller);
        obj = obj.addMode(plant);
        
        if i > 1
          data{i}.V0p = subs(data{i-1}.Vsol,[x;t],[xp;0]);
          
          time_guard_i = @(obj,t,x,u) time_guard(obj,t,x,u,data{i}.T);
          obj = obj.addTransition(i,time_guard_i,@transition,false,false);
        end
      end
      
      obj.data = data;
      obj.t = t;
      obj.x = x;
      obj.u = u;
      obj.s = s;
      obj.xp = xp;
      
      obj = obj.setSimulinkParam('Solver','ode4');
      obj = obj.setSimulinkParam('FixedStep','.001');
    end
    
    function phi = time_guard(obj,t,x,u,T_max)
      t = t - x(end-2);
      phi = T_max - t;
    end
    
    function [to_mode_xn,to_mode_num,status] = transition(obj,from_mode_num,t,x,u)
      r_orig = x(end-1:end);
      t_orig = t;
      t = t - x(end-2);
      x = x(1:end-3);
      to_mode_num = from_mode_num-1;
      status = 0;
      % extract s via exhaustive search
      [smin,smax] = simpleResetInputLimits(obj.sos_plant,x);
      if length(smin) > 1
        error('Not yet implemented')
      end
      
      if length(smin) == 1
        S = linspace(smin,smax,100);
        
        V0p = subs(obj.data{from_mode_num}.V0p,obj.x,x);
        V0p_vals = msubs(V0p,obj.s,S);
        [V0_opt,i_opt] = max(V0p_vals);
        t_orig
        V0_opt
        s_opt = S(i_opt)
      else
        s_opt = [];
      end
      
      to_mode_xn = [double(subs(obj.xp,[obj.t;obj.x;obj.s],[t;x;s_opt]));t_orig;obj.sos_plant.stanceUpdate(x,r_orig,s_opt)];
    end
  end
  
  % guard: start with fixed-time, implement state based later
  % transition: brute force it
  % visualizer: modify state first
end