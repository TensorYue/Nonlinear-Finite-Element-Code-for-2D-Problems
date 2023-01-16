% Solver

% Create Store 
StoreD=zeros(4,2,timestep+1);
StoreV=zeros(4,2,timestep+1);
StoreA=zeros(4,2,timestep+1);
StoreD(:,:,1)=Dinitial;
StoreV(:,:,1)=Vinitial;
StoreA(:,:,1)=Ainitial;

% Qbound % Process Q bound to be its own size
Qbound=zeros(noequation,1);
for i=1:size(Boundary_Q,1)
    j=2*Boundary_Q(i,1)-2+Boundary_Q(i,2);
    Qbound(j)=1;
end

% Hbound % Process H bound to be its own size
Hbound=zeros(noequation,1);
for i=1:size(Boundary_H,1)
    j=2*Boundary_H(i,1)-2+Boundary_H(i,2);
    Hbound(j)=Boundary_H(i,3);
end

for t = 1:timestep % Loop over time step
    dold=StoreD(:,:,t);
    vold=StoreV(:,:,t);
    aold=StoreA(:,:,t);
    d_predict=dold+dt*vold+dt^2/2*(1-2*beta)*aold; % Predictor of d
    v_predict=vold+(1-gama)*dt*aold; % Predictor of v
    Cd_predict=reshape(d_predict',8,1);
    dnew=d_predict;
    vnew=v_predict;
    anew=aold;

    for n=1:4
        Cdnew=reshape(dnew',8,1);
        Canew=reshape(anew',8,1);
        U=dnew;
        Evaluate_Integral; % Integral over element (Include Mass Matrix)
        Fin=Fint_element+M*Canew; % Evaluate Finternal
        Tan=Tangent_element*beta*dt^2+M; % Evaluate Consistent Tangent Matrix
    
        BoundaryMapping=1:1:noequation;
    
        CU=reshape(U',8,1);
        
        TU=CU;
        Tanew=Canew;
        CFin=Fin;
        CFex=Hbound;
        CTan=Tan;
    
        for i=1:noequation % Process Matrix to Eleminate Q Boundary
            if Qbound(noequation+1-i)==1
               CU(noequation+1-i)=[];
               Canew(noequation+1-i)=[];
               CFin(noequation+1-i)=[];
               CFex(noequation+1-i)=[];
               CTan(noequation+1-i,:)=[];
               CTan(:,noequation+1-i)=[];
               BoundaryMapping(noequation+1-i)=[];
            end
        end

        Re=CFex-CFin; % Calculate Residual
        if n==1
            normRe0=norm(Re);
        end

        dCanew=CTan\Re; % Calculate increment
        Canew=Canew+dCanew; % Update acceleration

        for i=1:length(BoundaryMapping) % Recover acceleration to its own size
            Tanew(BoundaryMapping(i))=Canew(i);
        end
    
        anew=reshape(Tanew,2,4)'; 
        vnew=v_predict+dt*gama*anew; % Update velocity
        dnew=d_predict+dt^2*beta*anew; % Update displacement

         if norm(Re)<tol*normRe0 % Convergence criteria
            outer_count=outer_count+1;
            n

           StoreD(:,:,t+1)=dnew;
           StoreV(:,:,t+1)=vnew;
           StoreA(:,:,t+1)=anew;
           break
        end
 
    end
end