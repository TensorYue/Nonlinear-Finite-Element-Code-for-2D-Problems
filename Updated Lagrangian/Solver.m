% Solver

% Residual Storer
StoreRe_10 = zeros(1,100);
Step_10 = 0;
StoreRe_50 = zeros(1,100);
Step_50 = 0;

% Create Store 
StoreU=zeros(4,2,controlstep+1);

% Qbound
Qbound=zeros(noequation,1);
for i=1:size(Boundary_Q,1)
    j=2*Boundary_Q(i,1)-2+Boundary_Q(i,2);
    Qbound(j)=1;
end

% Hbound
Hbound=zeros(noequation,1);
for i=1:size(Boundary_H,1)
    j=2*Boundary_H(i,1)-2+Boundary_H(i,2);
    Hbound(j)=Boundary_H(i,3);
end

for control=1:controlstep+1

    for n=1:100000
        Evaluate_Integral;
        Fin=Fint_element;
        Tan=Tangent_element;
    
        BoundaryMapping=1:1:noequation;
    
        CU=reshape(U',8,1);
        TU=CU;
        CFin=Fin;
        CFex=(control-1)/controlstep*Hbound;
        CTan=Tan;
    
        for i=1:noequation
            if Qbound(noequation+1-i)==1
               CU(noequation+1-i)=[];
               CFin(noequation+1-i)=[];
               CFex(noequation+1-i)=[];
               CTan(noequation+1-i,:)=[];
               CTan(:,noequation+1-i)=[];
               BoundaryMapping(noequation+1-i)=[];
            end
        end

        Re=CFex-CFin;

        % Store Residual

        if control==10+1
            StoreRe_10(n) = norm(Re);
            Step_10 = Step_10 + 1;
        end
        if control==50+1
            StoreRe_50(n) = norm(Re);
            Step_50 = Step_50 + 1;
        end

        if norm(Re)<tol
            outer_count=outer_count+1
            n

           StoreU(:,:,control)=U;
           break
        end

        dCU=CTan\Re;
        CU=CU+dCU;

        for i=1:length(BoundaryMapping)
            TU(BoundaryMapping(i))=CU(i);
        end
    
        U=(reshape(TU,2,4))';

    end
end