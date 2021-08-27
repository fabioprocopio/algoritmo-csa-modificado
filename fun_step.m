function fitness = fun_step(xx)
	d = size(xx, 2);
	
	soma = 0;
	for ii=1:d
		soma = soma + floor(xx(ii) + 0.5)^2;
	end;
	fitness = soma;