 function oposto_elite = EOBL(elite, vetor)
    oposto_elite = [];  
    d            = size(elite, 2);
    d_aj         = min(vetor);
    d_bj         = max(vetor);
    
    for i=1:d
        oposto_elite(i) = rand * (d_aj + d_bj) - elite(i);
        dimensao_atual  = oposto_elite(i);
        % Verificando se está dentro do limite
        if dimensao_atual < d_aj || dimensao_atual > d_bj
            % Escolhe um valor aleatório no intervalo [d_aj, d_bj]
            oposto_elite(i) = (d_bj - d_aj) * rand + d_aj;
        end
    end 