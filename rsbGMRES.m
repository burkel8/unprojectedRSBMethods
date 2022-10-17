function [resid] = rsbGMRES(A,B,X,shifts,m,k,s,n,tol,U,C,shift_recycle_method)

eyebar = zeros(m*s+s,m*s);
for i=1:m*s
eyebar(i,i)=1;
end

R=zeros(n,s);

for i=1:s
R(:,i) = B(:,i) - (A + shifts(i)*speye(n))*X(:,i);
end
iter = 1;
resid(iter) = norm(R);

shift_monitor = 1;


while(norm(R) > tol)

 CTU = C'*U;
 CTR = C'*R;
        for i = 1:length(shifts)
            rProj = ( eye(size(U,2)) + shifts(i)*CTU ) \ CTR(:,i);
            X(:,i) = X(:,i) + U*rProj;
            R(:,i) = R(:,i) - C*rProj - U*(shifts(i)*rProj);
        end


[W,H,S0,B] = rsbBlock_Arnoldi(A,R,m,C);


 WTU = W(:,1:(m+1)*s)'*U;
        Uhat = U - C*CTU - W(:,1:(m+1)*s)*WTU;
        Vhat = orth(Uhat);
        VTU = Vhat'*U;
        kr = size(Vhat,2);
        RES = zeros((m+1)*s+k+kr,s);
        RHS = eye((m+1)*s+k+kr,s);
        RHS(1:s,:) = S0;
        RHS(end-kr+1:end,:) = Vhat'*R;
        Y = zeros(m*s+k,s);
        
        HH = @(sigma) H(1:(m+1)*s,1:m*s) + sigma*eye((m+1)*s,m*s);
        G = @(sigma) [HH(sigma) sigma*WTU; B(:,1:m*s) eye(k)+sigma*CTU; zeros(kr,m*s) sigma*VTU];
        
        for p=1:s
            if shifts(p) ~= 0
                Gsig = G(shifts(p));
                Y(:,p) = Gsig\RHS(:,p);
                RES(:,p) = RHS(:,p) - Gsig*Y(:,p);
            else
                rhstemp = zeros((m+1)*s,1);
                rhstemp(1:s) = S0(:,p);
                ytemp = H(1:(m+1)*s,1:m*s)\rhstemp;
                Y(1:m*s,p) = ytemp;
                Y(m*s+1:end,p) = - B(:,1:m*s)*ytemp;
                restemp = rhstemp - H(1:(m+1)*s,1:m*s)*ytemp;
                RES(1:(m+1)*s,p) = restemp;
            end
        end

     

 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute approximation at end of cycle %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = X + W(:,1:m*s)*Y(1:m*s,:);
    X = X + U*Y(m*s+1:end,:);
    R = W(:,1:(m+1)*s)*RES(1:(m+1)*s,:);
    R = R + C*RES((m+1)*s+1:(m+1)*s+k,:);
    R = R + Vhat*RES((m+1)*s+k+1:end,:);

iter = iter + 1;
resid(iter) = norm(R);

shift = shifts(shift_monitor);
[U] = block_harm_ritz(H,W,m,s,k,U,C,shift);
C = A*U;
[C,RR] = qr(C,0);
U = U/RR;

 if shift_recycle_method == 1
    shift_monitor = shift_monitor + 1;
    if shift_monitor == size(shifts,2) 
        shift_monitor = 1;
    end
 end


end
end





