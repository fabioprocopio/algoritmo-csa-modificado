function fitness = fun_griewank(xx)
    d = length(xx);
    soma = 0;
    produto = 1;

    for ii = 1:d
        xi = xx(ii);
        soma = soma + xi^2/4000;
        produto = produto * cos(xi/sqrt(ii));
    end

    fitness = soma - produto + 1;
	
	