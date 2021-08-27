function fitness = fun_rosenbrock(xx)
    d = length(xx);
    soma = 0;
    for ii = 1:(d-1)
        xi = xx(ii);
        xnext = xx(ii+1);
        new = 100*(xnext-xi^2)^2 + (xi-1)^2;
        soma = soma + new;
    end
    fitness = soma;    

