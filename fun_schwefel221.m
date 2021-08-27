function fitness = fun_schwefel221(x)
    [SwarmSize, Dim] = size(x);
    for i=1:SwarmSize
        fitness(i,1) = max(abs(x(i,:)));
    end