function fitness = fun_rastrigin(xx)
    %fitness = length(xx) + sum(xx.*xx-cos(2*pi*xx));
	d = length(xx);
	sum = 0;
	for ii = 1:d
		xi = xx(ii);
		sum = sum + (xi^2 - 10*cos(2*pi*xi));
	end

	fitness = 10*d + sum;
