classdef NStepCapturabilityVisualizer < Visualizer
  properties
    plant
  end
  
  methods
    function obj = NStepCapturabilityVisualizer(plant)
      obj = obj@Visualizer(plant.getOutputFrame);
      obj.plant = plant;
    end
    
    function draw(obj,t,y)
      obj.plant.sos_plant.draw(t,y);
    end
  end
end

