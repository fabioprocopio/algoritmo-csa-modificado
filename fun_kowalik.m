function fitness = fun_kowalik(xx)
    fitness=0;
    [SwarmSize, Dim] = size(xx);
    a=[0.1957 0.1947 0.1735 0.16 0.0844 0.0627 0.0456 0.0342 0.0323 0.0235 0.0246];
    b=1./[0.25 0.5 1 2 4 6 8 10 12 14 16];
    for i=1:11
        fitness = fitness+(a(1,i)-(xx(:,1).*(b(1,i).^2+b(1,i)*xx(:,2)))./(b(1,i).^2+b(1,i)*xx(:,3)+xx(:,4))).^2;
end