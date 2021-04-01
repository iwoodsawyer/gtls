% Example with mixed weighted disturbances

n = 10000; % number of data samples

R = rand(16,16); % random noise correlation
E = 0.01.*randn(n,16)*R;
dA = E(:,1:11);
dB = E(:,12:16);

A0 = randn(n,21);
A1 = A0(:,1:10); % error-free
A2 = A0(:,11:21) + dA;

X0 = randn(21,5); 

B0 = A0*X0;
B = B0 + dB;

W = cov([dA dB]); % covariance of disturbances

% perform TLS
tic
X = tls([A1 A2],B);
toc
norm(X-X0,'fro')

% perform MTLS
tic
X = mtls(A1,A2,B);
toc
norm(X-X0,'fro')

% perform GLS
tic
X = gtls([A1 A2],B,blkdiag(0.001.*eye(10),W));
toc
norm(X-X0,'fro')

% perform GMTLS
tic
X = gmtls(A1,A2,B,W);
toc
norm(X-X0,'fro')



