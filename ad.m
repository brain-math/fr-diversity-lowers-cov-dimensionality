function t = ad(data1, data2)
%AD Calculate the Anderson-Darling 2-sample statistic
m = length(data1); n = length(data2); k = m + n;
data = [data1; data2]; data = sort(data);
t = 0;
for i = 1 : k-1
    c = nnz(data1 <= data(i));
    t = t + (k*c-m*i)^2/i/(k-i);
end
t = t/m/n;

end
