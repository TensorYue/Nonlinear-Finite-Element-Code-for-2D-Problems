function [Na_X, J_cauchy] = Finite_Derivative(nodel_position,cauchy_position)

% Evaluate integral point location
cau = cauchy_position(1);
yi = cauchy_position(2);

% Evaluate Na,cauchy
cauchy_table = [-1 -1;1 -1;1 1;-1 1];

Na_cauchy = zeros(4,2); % nen number of element nodal = 4
for i=1:4
    cauchy = cauchy_table(i,1);
    yita = cauchy_table(i,2);
    Na_cauchy(i,1) = 0.25*cauchy*(1+yi*yita);
    Na_cauchy(i,2) = 0.25*yita*(1+cau*cauchy);
end

% Evaluate X,cauchy
X_cauchy = nodel_position'*Na_cauchy;

% Evaluate J_cauchy
J_cauchy = det(X_cauchy);

% Evaluate Na,X
Na_X = Na_cauchy*[X_cauchy(2,2) -X_cauchy(1,2); -X_cauchy(2,1) X_cauchy(1,1)]/J_cauchy;

end