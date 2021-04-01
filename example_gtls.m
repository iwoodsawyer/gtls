% Example with weighted disturbances

n = 10000; % number of data samples

R = rand(26,26); % random noise correlation
E = 0.01.*randn(n,26)*R;
dA = E(:,1:21);
dB = E(:,22:26);

A0 = randn(n,21);
A = A0 + dA;
A1null = zeros(n,0);

X0 = randn(21,5); 

B0 = A0*X0;
B = B0 + dB;

W = cov([dA dB]); % covariance of disturbances

% perform TLS
tic
X = tls(A,B);
toc
norm(X-X0,'fro')

% perform MTLS
tic
X = mtls(A1null,A,B);
toc
norm(X-X0,'fro')

% perform GTLS
tic
X = gtls(A,B,W);
toc
norm(X-X0,'fro')

% perform GMTLS
tic
X = gmtls(A1null,A,B,W);
toc
norm(X-X0,'fro')