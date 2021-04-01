function [X,A0,B0,dA,dB] = tls(A,B,varargin)
%TLS Total Least Squares
%   X=TLS(A,B) solves the total least squares problem, also known as errors
%   in variables, formulated in the overdetermined set of linear equations
%   (A0+dA)X = (B0+dB), where covariance matrix of the disturbances is
%   expected to be diagonal and denoted by E([dA dB]^T[dA dB]) = sigma_d*I.
%   The computation is based on SVD([A B]) and any singular values less
%   than a tolerance are treated as zero. The default tolerance is
%   MAX(SIZE(A))*EPS(NORM(A)).
%
%   X=TLS(A,B,TOL) uses the tolerance TOL instead of the default.
%
%   X=TLS(A,B,TOL,'SV') plot the singular values wrt to tolerance value.
%
%   [X,A0,B0]=TLS(A,B) returns the estimated A0 and B0 matrices.
%
%   [X,A0,B0,dA,dB]=TLS(A,B) returns the estimated dA and dB disturbances.
%
%   See also MTLS, GTLS, GMTLS. 
%
%   References:
%   [1] S. van Huffel, and J. Vandenwalle, "The Total Least Squares
%   Problem: Computational Aspects and Analysis", SIAM, 1991.
%   [2] S. van Huffel, and J. Vandenwalle, "Analysis and properties of the
%   generalized total least squares problem AX=B when some or all columns
%   in A are subject to error", SIAM J. Matrix Anal. Appl., Vol.10, No.3,
%   pp. 294-315, 1989.

% check input arguments
if nargin < 2 
    error('TLS requires at least two input arguments')
end
[m,n] = size(A);
[mb,d] = size(B);
if ~isequal(m,mb)
    error('The number of rows of matrix A and B must be the same size.')
end
if n > m
   warning('Problem is underdetermined (m < n). Solution is not unique and will return minimum norm solution.')
end 

R22 = [A B];
if ((nargout == 1) && (m > (5/3)*(n+d)))
    % If only solution X is returned, it is more efficient to transform   
    % [A B] to upper triangular form R first using QR factorization.
    R22 = triu(qr(R22,0));
end

[U,S,V] = svd(R22,0);
if m > 1, s = diag(S);
elseif m == 1, s = S(1);
else , s = 0;
end

if ((nargin >= 3) && ~isempty(varargin{1}))
    tol = varargin{1};
else
    tol = max(m,n)*eps(max(s));
end
if nargin == 4  % plot singular values wrt to tolerance
    figure, 
    semilogy(1:length(s),s,'x',[1; length(s)],[tol; tol],'-r');
    grid;
end

p = min(sum(s > tol),n); % rank determination
if p == 0
    X = zeros(n,d,class(A));
else
    if p < n % rank deficient
        V2 = rq(V(:,p+1:end),0);
        r = abs(diag(V2(p+1:n+d,:)));
        rtol = max(m,n)*eps(max(r));
        v = sum(r < rtol);
        if v > 0 % solution is non generic -> lower the rank
            p = p - v;
            V2 = rq(V(:,p+1:end),0);
        end
    else % full rank
        V2 = V(:,p+1:end);
    end
    V12 = V2(1:n,end-d+1:end);
    V22 = V2(n+1:end,end-d+1:end);
    X = -V12/V22;
    if nargout > 1
        AB = U(:,1:p)*S(1:p,1:p)*V(:,1:p)';
        A0 = AB(:,1:n);
        B0 = AB(:,n+1:end);
    end
    if nargout > 3
        AB = U(:,p+1:end)*S(p+1:end,p+1:end)*V(:,p+1:end)';
        dA = AB(:,1:n);
        dB = AB(:,n+1:end);
    end
end











