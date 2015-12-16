function createRotatingVideo(hFig, filename)
writer = VideoWriter(filename);
set(writer, 'Quality', 90)

open(writer);
dtheta = 1;
for i = 1 : 360 / dtheta
  camorbit(dtheta, 0); drawnow;
  writeVideo(writer, getframe(hFig));
end
close(writer);
end
