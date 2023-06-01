function ML = ML_init(z,awzp,N,LogC)
% ML = ML_init(z,awzp,N,LogC)
% Returns: ML - Logarithmic operator as NxN matrix
    ML = awzp.'.*log(abs(z.'-z));
    ML(1:N+1:N^2) = 0;
    ML = ML+LogC.*awzp.';
end