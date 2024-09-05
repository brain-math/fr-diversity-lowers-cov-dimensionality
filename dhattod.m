function D = dhattod(Dhat, alpha)
%DHATTOD Eliminate the effect of temporal sampling
D = Dhat ./ (1-alpha*Dhat);
    
end
