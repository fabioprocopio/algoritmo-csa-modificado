function [qntdColum, local] = verifica_funcao_grafico(funcao)

    if strcmp(funcao, 'fun_esfera')
        qntdColum = 4;
        local = 'best';


    elseif strcmp(funcao, 'fun_schwefel222')
        qntdColum = 4;
        local = 'best';


    elseif strcmp(funcao, 'fun_schwefel12')
        qntdColum = 4;
        local = 'best';


    elseif strcmp(funcao, 'fun_rosenbrock')
        qntdColum = 5;
        local = 'best';

    elseif strcmp(funcao, 'fun_step')
        qntdColum = 1;
        local = 'northeast';


    elseif strcmp(funcao, 'fun_rastrigin')
        qntdColum = 3;
        local = 'southeast';


    elseif strcmp(funcao, 'fun_ackley')
        qntdColum = 4;
        local = 'best';

    elseif strcmp(funcao, 'fun_griewank')
        qntdColum = 3;
        local = 'southeast';
        
    elseif strcmp(funcao, 'fun_goldstein_price')
        qntdColum = 5;
        local = 'best';
        
    elseif strcmp(funcao, 'fun_hartmann3')
        qntdColum = 5;
        local = 'best';
        
    elseif strcmp(funcao, 'fun_hartmann6')
        qntdColum = 5;
        local = 'best';
        
    elseif strcmp(funcao, 'fun_hartmann6')
        qntdColum = 5;
        local = 'best';
    
    elseif strcmp(funcao, 'fun_shekel5')
        qntdColum = 5;
        local = 'best'; 
    elseif strcmp(funcao, 'fun_shekel7')
        qntdColum = 5;
        local = 'best';
    elseif strcmp(funcao, 'fun_shekel10')
        qntdColum = 5;
        local = 'best';
    else
        qntdColum = 4;
        local = 'northoutside';
    end;