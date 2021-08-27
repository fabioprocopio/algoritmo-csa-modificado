function vetor_mod = GOBL(vetor, lower, upper)
    vetor_mod = [];
    d = size(vetor, 2);
    
    a_j = min(vetor);
    b_j = max(vetor);
    for i=1:d
        k = rand();    
        vetor_mod(i) = k * (a_j + b_j) - vetor(i);
        
        if vetor_mod(i) > upper || vetor_mod(i) < lower
           vetor_mod(i) = a_j + (b_j - a_j) * rand();
        end
    end