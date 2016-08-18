function [ G ] = graspMatrix( contactPoints, objCenter)
%GRASPMATRIX this function computes the grasp matrix given the
%contactPoints.

G = [];
for i = 1:length(contactPoints)
    cp = contactPoints{i};
    r = cp.pos - objCenter;
    S = -[   0, -r(3), r(2) ;...
            r(3), 0, -r(1) ;...
            -r(2), r(1), 0];
    R = cp.contactFrame;
    
    Ai = [  R, zeros(3,3);
            S*R, R];
    Gi = Ai*cp.selectionMatrix;
        
    G = [G, Gi];
end

end

