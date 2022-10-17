function [W,H,RR] = block_Arnoldi(A,R,m,s,n)

    W = zeros(n,(m+1)*s);
    H = zeros((m+1)*s,m*s);
   
    [W(:,1:s),RR] = qr(R,0);
    
     for k = 1:s:m*s
        V = A*W(:,k:k+s-1);
        W(:,k+s:k+2*s-1) = V;
       
        for j=1:s:k
            H(j:j+s-1,k:k+s-1) = W(:,j:j+s-1)'*W(:,k+s:k+2*s-1);
            W(:,k+s:k+2*s-1) = W(:,k+s:k+2*s-1) - W(:,j:j+s-1)*H(j:j+s-1,k:k+s-1);
        end
        [W(:,k+s:k+2*s-1),H(k+s:k+2*s-1,k:k+s-1)] = qr(W(:,k+s:k+2*s-1),0);

    end

end