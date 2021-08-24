% This file compares takes a matrix A and converts it to lower triangular
% Then, calculating determinant is easy, and stored in d
% Then a fadlev algorithm is used to compare to
% And finally, matlab's det function is used just to make it feel bad
function [d,fadDet,matDet] = testFadlevAccuracy(A)
    A = tril(A);
    d = vpa(1);
    [n,~] = size(A);
    for k = 1:n
        d = d*A(k,k);
    end
    p = fadlevMP(A); % to test an algorithm, just change this
    fadDet = p(n+1)*(-1)^n;
    matDet = det(A);
end
