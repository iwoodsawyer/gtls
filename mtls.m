function X = mtls(A1,A2,B,varargin)
%MTLS Mixed Total Least Squares
%   X=MTLS(A1,A2,B) solves the mixed total least squares problem, also
%   known as errors in variables, formulated in the overdetermined set of
%   linear equations [A1 A2]X = B, where A1 are the error-free variables, 
%   and A2 = A0 + dA2 and B = B0 + dB are the variables with disturbances.
%   The covariance matrix of the disturbances is expected to be diagonal
%   and denoted by E([dA2 dB]^T[dA2 dB]) = sigma_d.*I. The computation is
%   based on SVD([A1 A2 B]) and any singular values less than a tolerance
%   are treated as zero. The default tolerance is MAX(SIZE(A))*EPS(NORM(A)).
%
%   X=MTLS(A1,A2,B,TOL) uses the tolerance TOL instead of the default.
%
%   X=MTLS(A1,A2,B,TOL,'SV') plot the singular values wrt to tolerance value.
%
%   See also TLS, WTLS, GMTLS. 
%
%   References:
%   [1] S. van Huffel, and J. Vandenwalle, "The Total Least Squares
%   Problem: Computational Aspects and Analysis", SIAM, 1991.
%   [2] S. van Huffel, and J. Vandenwalle, "Analysis and properties of the
%   generalized total least squares problem AX=B when some or all columns
%   in A are subject to error", SIAM J. Matrix Anal. Appl., Vol.10, No.3,
%   pp. 294-315, 1989.

% check input arguments
if nargin < 3
    error('MTLS requires at least three input arguments')
end
[m1,n1] = size(A1);
[m2,n2] = size(A2);
[mb,d] = size(B);
if ~isequal(m1,m2,mb)
    error('The number of rows of matrix A1, A2, and B must be the same size.')
end
m = m1 + m2;
n = n1 + n2;
if n > m
   warning('Problem is underdetermined (m < n). Solution is not unique and will return minimum norm solution.')
end 

A = [A1 A2];
if n1 > 0
    R = triu(qr([A B],0));
    R11 = R(1:n1,1:n1);
    R12 = R(1:n1,n1+1:n1+n2+d);
    R22 = R(n1+1:end,n1+1:n1+n2+d);
    
    [~,S,V] = svd(R22,0);
    if m > 1, s = diag(S);
    elseif m == 1, s = S(1);
    else , s = 0;
    end  
    
    if ((nargin >= 4) && ~isempty(varargin{1}))
        tol = varargin{1};
    else
        tol = max(m,n)*eps(max(s));
    end
    if nargin == 5  % plot singular values wrt to tolerance
        figure,
        semilogy(1:length(s),s,'x',[1; length(s)],[tol; tol],'-r');
        grid;
    end
    
    p = min(sum(s > tol),n2); % rank determination
    if p == 0
        X = zeros(n,d,class(A));
    else      
        if p < n2 % rank deficient
            V22 = V(:,p+1:end);
            V12 = R11\(-R12*V22);
            V2 = [V12; V22];
            V2 = rq(V2,0);
            r = abs(diag(V2(p+1:n+d,:)));
            rtol = max(m,n)*eps(max(r));
            v = sum(r < rtol);
            if v > 0 % solution is non generic -> lower the rank
                p = p - v;
                V22 = V(:,p+1:end);
                V12 = R11\(-R12*V22);
                V2 = [V12; V22];
                V2 = rq(V2,0);
            end
        else % full rank
            V22 = V(:,p+1:end);
            V12 = R11\(-R12*V22);
            V2 = [V12; V22];
        end
        
        V12 = V2(1:n,end-d+1:end);
        V22 = V2(n+1:end,end-d+1:end);
        X = -V12/V22;
    end
else
    % Do full total least squares
    if nargin >= 5
        X = tls(A,B,varargin{1},varargin{2});
    elseif nargin == 4
        X = tls(A,B,varargin{1});
    else
        X = tls(A,B);
    end 
end



