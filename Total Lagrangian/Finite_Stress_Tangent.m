function [P, S, Cmat] = Finite_Stress_Tangent(F, J)

% Material Property
E = 100; % (MPa)
v = 0.25;
lambda = E*v/(1+v)/(1-2*v);
miu = E/2/(1+v);

% Evaluate Cauchy Green Deformation tensor
C = F'*F;

% Evaluate inv(F)
invF = [F(2,2) -F(1,2); -F(2,1) F(1,1)]/J;

% Evaluate inv(C)
invC = invF*invF';

% Evaluate PK1 (varies with material)
P = miu*F + (lambda*log(J)-miu)*invF';

% Evaluate PK2 (varies with material)
S = invF*P;

% Evaluate Cmat
Cmat = [getCmat(invC,J,1,1,1,1) getCmat(invC,J,1,1,2,2) getCmat(invC,J,1,1,1,2)
        getCmat(invC,J,2,2,1,1) getCmat(invC,J,2,2,2,2) getCmat(invC,J,2,2,1,2)
        getCmat(invC,J,1,2,1,1) getCmat(invC,J,1,2,2,2) getCmat(invC,J,1,2,1,2)];

end

function [Cmat_ijkl] = getCmat(invC,J,i,j,k,l)

% Material Property
E = 100; % (MPa)
v = 0.25;
lambda = E*v/(1+v)/(1-2*v);
miu = E/2/(1+v);

Cmat_ijkl = lambda*invC(i,j)*invC(k,l) + (miu-lambda*log(J))*(invC(i,k)*invC(j,l)+invC(i,l)*invC(j,k));

end

