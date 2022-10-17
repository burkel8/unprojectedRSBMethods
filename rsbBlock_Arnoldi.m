function [W,H,S0,B] = rsbBlock_Arnoldi(A,R,m,C)


    [n,s] = size(R);
    k = size(C,2);
    H = zeros((m+1)*s,m*s);
    W = zeros(n,(m+1)*s);
    B = zeros(k,m*s);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get first block Arnoldi vector %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [W(:,1:s),S0] = qr(R,0);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize RHS for    %
    % least squares problem %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    RHS = zeros((m+1)*s,s);
    %keyboard

    iter = 0;
    %%%%%%%%%%%%%%%%%%%%
    % Iterate by block %
    %%%%%%%%%%%%%%%%%%%%
    for pp = 1:s:m*s
        iter = iter + 1;
        
        V = A*W(:,pp:pp+s-1);
        W(:,pp+s:pp+2*s-1) = V;
        
        %%%%%%%%%%%%%%%%%%%%%
        % Orthog against C  %
        % Store coeffs in B %
        %%%%%%%%%%%%%%%%%%%%%
        B(:,pp:pp+s-1) = C(:,1:k)'*W(:,pp+s:pp+2*s-1);
        W(:,pp+s:pp+2*s-1) = W(:,pp+s:pp+2*s-1) - C(:,1:k)*B(:,pp:pp+s-1);
        
        %%%%%%%%%%%%%%%%%
        % Block Arnoldi %
        % Process       %
        %%%%%%%%%%%%%%%%%
        for j=1:s:pp
            H(j:j+s-1,pp:pp+s-1) = W(:,j:j+s-1)'*W(:,pp+s:pp+2*s-1);
            W(:,pp+s:pp+2*s-1) = W(:,pp+s:pp+2*s-1) - W(:,j:j+s-1)*H(j:j+s-1,pp:pp+s-1);
        end
        [W(:,pp+s:pp+2*s-1),H(pp+s:pp+2*s-1,pp:pp+s-1)] = qr(W(:,pp+s:pp+2*s-1),0);
    end



end
