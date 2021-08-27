function fitness = fun_quartic_noisy(xx)   	
	
	d = length(xx);
	f = zeros(1, 1);
	for ii = 1:d
		f = ii*xx(1:1, ii).^4 + f;
	end	
	fitness = f + rand(1, 1);