% Evaluate_Integral

% Three point quadrature
cauchy=[-sqrt(1/3) sqrt(1/3)];
w=[1 1];

% Elementary Fint and Kint
Fint_element=zeros(8,1);
Tangent_element=zeros(8,8);

for i=1:2
    for j=1:2

        cauchy_position=[cauchy(i) cauchy(j)];

        [Na_X, J_cauchy] = Finite_Derivative(nodel_position,cauchy_position);

        [F, J] = Finite_Deformation(Na_X,U);

        [P, S, Cmat] = Finite_Stress_Tangent(F, J);

        [Fint, Kint] = Finite_Internal_Force_Tangent(Na_X,F,P,S,Cmat);
         
        Fint_element=Fint_element+Fint*w(i)*w(j)*J_cauchy;

        Tangent_element=Tangent_element+Kint*w(i)*w(j)*J_cauchy;
        
    end 
end

% Mass Matrix
M_cell=cell(4,4);
dense=1;
Me = (J_cauchy*4)*dense*diag([0.25 0.25 0.25 0.25]); % Lumped mass
%Me = (J_cauchy*4)*dense*[0.11111 0.05556 0.02778 0.05556
%                        0.05556 0.11111 0.05556 0.02778
%                        0.02778 0.05556 0.11111 0.05556
%                        0.05556 0.02778 0.05556 0.11111]; % Consistent mass
for I=1:4
    for J=1:4
        M_cell{I,J}=[Me(I,J) 0;0 Me(I,J)];
    end
end
M=cell2mat(M_cell);