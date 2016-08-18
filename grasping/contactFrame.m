function [ contactFrame ] = contactFrame( contactNormal )
%CONTACTFRAME given the normal of a contact, compute the contact frame
% based on the global frame (it returns a rotation matrix)
%   This function tries to put the y axis of the contact frame in the
%   the y axis of the global frame if the given normal is not "very"  
%   parallel to the x axis. If it is then it tries to align the x with the
%   x.

if abs(contactNormal'*[1;0;0]) < 0.8
    z = contactNormal;
    y = cross(contactNormal, [1;0;0]);
    y = y/norm(y);
    x = cross(y, z);
    x = x/norm(x);
else
    z = contactNormal;
    x = cross(contactNormal, [0;1;0]);
    x = x/norm(x);
    y = cross(z, x);
    y = y/norm(y);
end

contactFrame = [x,y,z];

end

