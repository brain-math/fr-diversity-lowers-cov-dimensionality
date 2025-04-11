function V = triu_ele(A)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

n = size(A, 1);
I = repmat(1:n, [n, 1]);
J = I';
V = A(J<I);

end

