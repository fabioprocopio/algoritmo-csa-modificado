function fitness = fun_schwefel222(x) %Schwefel 2.22
    n=size(x,2);
    sum1 = 0;
    sum2=1;
    for i=1:n
        sum1=sum1+(x(i)^2)^0.5;
        sum2=sum2*(x(i)^2)^0.5;
    end
    fitness = sum1+sum2;