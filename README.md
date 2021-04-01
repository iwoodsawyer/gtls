These functions calculates the solution for the generalized and/or mixed total least squares problem.

The total least squares problem, also known as errors in variables, solves the over-determined set of linear equations (A0+dA)X = (B0+dB), where covariance matrix of the unknown disturbances dA and dB is considered to be diagonal and denoted by E([dA dB]^T[dA dB]) = sigma_d.*I.

The mixed total least squares problem solves the over-determined set of linear equations [A1 A2]X = B, where A1 are the error-free variables, and A2 = A0 + dA2 and B = B0 + dB are the variables with disturbances.

The generalized total least squares problem solves the over-determined set of linear equations (A0 + dA)X = (B0 + dB), where the covariance matrix of the disturbances dA and dB is positive definite and given by sigma_d.*W = E([dA dB]^T[dA dB]).

[![View Total Least Squares with mixed and/or weighted disturbances on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/50332-total-least-squares-with-mixed-and-or-weighted-disturbances)
