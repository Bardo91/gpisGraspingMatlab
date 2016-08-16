function showAxis( origin, frameMatrix)
%SHOWAXIS This method receives a rotation matrix and plot it in the current
%figure.

v1 = frameMatrix(:,1);
v2 = frameMatrix(:,2);
v3 = frameMatrix(:,3);

hold on;
quiver3(origin(1), origin(2), origin(3), v1(1), v1(2), v1(3), 'r');
quiver3(origin(1), origin(2), origin(3), v2(1), v2(2), v2(3), 'g');
quiver3(origin(1), origin(2), origin(3), v3(1), v3(2), v3(3), 'b');

end

