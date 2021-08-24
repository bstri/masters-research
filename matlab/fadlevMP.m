function [p,Ainv] = fadlevMP(A)
% 
% Faddeev-Leverrier approach to generate coefficients of the
% characteristic polynomial and inverse of a given square matrix
%
% Input:     A - a square matrix
% Output:    p - the coefficient array of the characteristic polynomial 
%               p(1) = 1
%               p(2) = trace(B{1}) = trace(A)
%               p(3) = trace(B{2})/2
%               ...
%            Ainv - the inverse of A
%               Ainv = (B{n-1} - p(n)*I)/p(n+1)
% Example: 
%   mp.Digits(100)
%   A = randi(10,30) 
%   [p, Ainv] = fadlevMP(A)
%
    [n,~] = size(A);
    A = mp(A);
    I = eye(n);
    B = A;
    p = mp(ones(1,n+1));
    for k = 2:n
        % disp(k);
        p(k) = -trace(B)/(k-1);
        M = B + p(k)*I; 
        B = A*M;
    end
    p(n+1) = -trace(B)/n;
    Ainv = -M/p(n+1); 
    % detA = p(n+1)*(-1)^n; % uncomment if needed
end
