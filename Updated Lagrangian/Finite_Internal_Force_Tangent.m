function [Fint, Kint] = Finite_Internal_Force_Tangent(Na_X,F,P,S,cmat)

% Evaluate Fint
% Fint_mat = Na_X*P;
% Fint = reshape(Fint_mat',[],1); % Take care of reshape operation, it goes with row(first) - column

% Evaluate B matrix (Boolean matrices)
B = zeros(3,2*4); % nen = 4
% x_X=F(1,1);
% x_Y=F(1,2);
% y_X=F(2,1);
% y_Y=F(2,2);
% for i=1:4
%     Na_X_loop = Na_X(i,1);
%     Na_Y_loop = Na_X(i,2);
%     B(:,2*i-1) = [Na_X_loop*x_X;Na_Y_loop*x_Y;Na_X_loop*x_Y+Na_Y_loop*x_X];
%     B(:,2*i) = [Na_X_loop*y_X;Na_Y_loop*y_Y;Na_X_loop*y_Y+Na_Y_loop*y_X];
% end

% Test for current
invF = [F(2,2) -F(1,2); -F(2,1) F(1,1)]/det(F);
Na_x = Na_X*invF;
for i=1:4
    Na_x_loop = Na_x(i,1);
    Na_y_loop = Na_x(i,2);
    B(:,2*i-1) = [Na_x_loop;0;Na_y_loop];
    B(:,2*i) = [0;Na_y_loop;Na_x_loop];
end

sigma = 1/det(F)*P*F';
sigma_voigt = [sigma(1,1);sigma(2,2);sigma(1,2)];
Fint = B'*sigma_voigt;
Kint = B'*cmat*B; % current c

% Evaluate Kmat
% Kmat = B'*Cmat*B;

% Evaluate Kgeo
% Kgeo_cell = cell(4,4);

% for i=1:4
%     for j=1:4
%         Kgeo_cell{i,j} = eye(2)*(Na_X(i,:)*S*Na_X(j,:)');
%     end
% end
% Kgeo = cell2mat(Kgeo_cell);

% Evaluate Kint(only use Kmat)
% Kint = Kmat + Kgeo;

end