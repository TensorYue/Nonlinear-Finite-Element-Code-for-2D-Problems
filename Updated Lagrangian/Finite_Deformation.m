function [F, J] = Finite_Deformation(Na_X,U)

% Evaluate H
H = U'*Na_X;

% Evaluate F
I = eye(2); % nod number of degree = 2
F = I + H;

% Evaluate J
J = det(F);

end
