function dim_return = verifica_funcao_dim_fixa(funcao, dim)
    if strcmp(funcao, 'fun_hartmann6')
        dim_return = 6;
        
    elseif strcmp(funcao, 'fun_hartmann3')
        dim_return = 3;
        
    elseif strcmp(funcao, 'fun_goldstein_price')
        dim_return = 2;
        
    elseif strcmp(funcao, 'fun_kowalik')
        dim_return = 4;
    
    elseif strcmp(funcao, 'fun_six_hump')
        dim_return = 2;
        
    elseif strcmp(funcao, 'fun_shekel5')	
 		dim_return = 4;
        
    elseif strcmp(funcao, 'fun_shekel7')	
        dim_return = 4;
        
    elseif strcmp(funcao, 'fun_shekel10')	
 		dim_return = 4;
    else
        dim_return = dim;
    end;
