function [ G ] = graspMatrix( contactPoints )
%GRASPMATRIX this function computes the grasp matrix given the
%contactPoints.

G = [];
for i = 1:length(contactPoints)
    cp = contactPoints{i};
    r = cp.pos;
    S = [   0, -r(3), r(2) ;...
            r(3), 0, -r(1) ;...
            -r(2), r(1), 0];
    R = cp.contactFrame;
    
    Gi = [  R, zeros(3,3);
            S*R, R];

    G = [G; Gi];
end

end

