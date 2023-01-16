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
stress = 1/J*P*F';

% Evaluate Cmat
Cmat = [getcmat(F,invC,J,1,1,1,1)+stress(1,1) getcmat(F,invC,J,1,1,2,2) getcmat(F,invC,J,1,1,1,2)+stress(1,2)
        getcmat(F,invC,J,2,2,1,1) getcmat(F,invC,J,2,2,2,2)+stress(2,2) getcmat(F,invC,J,2,2,1,2)
        getcmat(F,invC,J,1,2,1,1)+stress(2,1) getcmat(F,invC,J,1,2,2,2) getcmat(F,invC,J,1,2,1,2)+stress(2,2)];

end

function [Cmat_ijkl] = getCmat(invC,J,i,j,k,l)

% Material Property
E = 100; % (MPa)
v = 0.25;
lambda = E*v/(1+v)/(1-2*v);
miu = E/2/(1+v);

Cmat_ijkl = lambda*invC(i,j)*invC(k,l) + (miu-lambda*log(J))*(invC(i,k)*invC(j,l)+invC(i,l)*invC(j,k));

end

function [cmat_ijkl] = getcmat(F,invC,detJ,i,j,k,l)
cmat_ijkl = 0;
for I=1:2
    for J=1:2
        for K=1:2
            for L=1:2
                cmat_ijkl = cmat_ijkl + 1/detJ*F(i,I)*F(j,J)*F(k,K)*F(l,L)*getCmat(invC,detJ,I,J,K,L);
            end
        end
    end
end
end

