classdef TwoDZMPVisualizer < Visualizer
  
  properties
  end
  
  methods
    function obj = TwoDZMPVisualizer(plant)
      obj = obj@Visualizer(plant.getOutputFrame);
    end
    
    function draw(obj,t,y)
      h=line(y(3)+[0;y(1)],[0;1]);
      set(h,'LineWidth',3,'Color','red')
      h=line(y(3)+[y(1);y(1)+y(2)],[1;0]);
      set(h,'LineWidth',3,'Color','black')
      
      xlim([-1 1])
      ylim([-1 1])
    end
  end
  
end

