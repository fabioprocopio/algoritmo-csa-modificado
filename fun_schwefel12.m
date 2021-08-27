function fitness = fun_schwefel12(x)
	fitness = 0;
	n = size(x, 2);
    for i=1:n
    	fitness = fitness + sum(x(1:i))^2;
    end
end