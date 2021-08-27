function [low, up] = get_espaco_busca(funcao)	
	
    if strcmp(funcao, 'fun_ackley')
        low       = -32;
        up        = 32;
    elseif strcmp(funcao, 'fun_alpine')	
 		low       = -10;        
        up        = 10;
    elseif strcmp(funcao, 'fun_esfera')
        low       = -100;        
        up        = 100;           
    elseif strcmp(funcao, 'fun_griewank')
        low       = -600;        
        up        = 600;    
    elseif strcmp(funcao, 'fun_powel_sum')	
 		low       = -1;        
        up        = 1;     
	elseif strcmp(funcao, 'fun_rastrigin')
        low       = -5.12;        
        up  	  = 5.12; 
    elseif strcmp(funcao, 'fun_rosenbrock')
        low       = -30;        
        up        = 30;                             
    elseif strcmp(funcao, 'fun_salomon')	
 		low       = -100;
        up        = 100; 
	elseif strcmp(funcao, 'fun_schumer_steiglitz')	
 		low       = -100;        
        up        = 100;      
    elseif strcmp(funcao, 'fun_zakharov')
        low       = -5;        
        up        = 10;                                 	
    elseif strcmp(funcao, 'fun_quartic_noisy')	
		low       = -1.28;        
        up        = 1.28;         
    elseif strcmp(funcao, 'fun_rotated_hyper_ellipsoid')	
 		low       = -65.536;        
        up        = 65.536; 
	elseif strcmp(funcao, 'fun_weierstrass')	
 		low       = -1;        
        up        = 1;
    elseif strcmp(funcao, 'fun_step')	
 		low       = -100;        
        up        = 100;
 	
    elseif strcmp(funcao, 'fun_schwefel221')	
 		low       = -100;        
        up        = 100;
    elseif strcmp(funcao, 'fun_schwefel222')	
 		low       = -10;        
        up        = 10;
    elseif strcmp(funcao, 'fun_schwefel12')	
 		low       = -100;        
        up        = 100;
    elseif strcmp(funcao, 'fun_schwefel')	
 		low       = -500;        
        up        = 500;
    elseif strcmp(funcao, 'fun_goldstein_price')	
 		low       = -2;        
        up        = 2;
    elseif strcmp(funcao, 'fun_hartmann3')	
 		low       = 0;        
        up        = 1;
    elseif strcmp(funcao, 'fun_hartmann6')	
 		low       = 0;        
        up        = 1;
    elseif strcmp(funcao, 'fun_kowalik')	
 		low       = -5;        
        up        = 5;
        
    elseif strcmp(funcao, 'fun_six_hump')	
 		low       = -5;        
        up        = 5;
        
    elseif strcmp(funcao, 'fun_shekel5')	
 		low       = 0;        
        up        = 10;
        
    elseif strcmp(funcao, 'fun_shekel7')	
 		low       = 0;        
        up        = 10;
        
    elseif strcmp(funcao, 'fun_shekel10')	
 		low       = 0;        
        up        = 10;
	end;