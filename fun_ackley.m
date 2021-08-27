function fitness = fun_ackley(xx)
    %Valores recomendados para as constantes
    a = 20;
    b = 0.2;
    c = 2 * pi;
    
    d     = length(xx);
    soma1 = sum(xx.^2);
    soma2 = sum(cos(c*xx));

    termo1 = -a * exp(-b*sqrt(soma1/d));
    termo2 = -exp(soma2/d);
	fitness = termo1 + termo2 + a + exp(1);
	
	
    