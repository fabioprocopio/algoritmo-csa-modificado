%Fonte: http://www.ie.itcr.ac.cr/rpereira/mat_ant/Genetic%20Algorithms/AII.pdf
function convergencia = GA(npar, popsize, maxit, funcao)

[varlo, varhi] = get_espaco_busca(funcao);
funcao = str2func(funcao);

%_______________________________________________________
mutrate   = 0.1;       % set mutation rate
selection = 0.5;       % fraction of population kept
Nt        = npar;      % continuous parameter GA Nt=#variables
keep      = floor(selection*popsize); % #population
% members that survive
nmut      = ceil((popsize-1)*Nt*mutrate); % total number of mutations
mating    = ceil((popsize-keep)/2);       % number of matings

%_______________________________________________________
% Create the initial population
par=(varhi-varlo)*rand(popsize,npar)+varlo; % random
for i=1:popsize
    cost(i) = funcao(par(i,:));          
end

[cost,ind] = sort(cost); % min cost in element 1
par        = par(ind,:); % sort continuous
minc(1)    = min(cost);  % minc contains min of
meanc(1)   = mean(cost); % meanc contains mean of population

% Localiza o melhor valor de fitness e o melhor indivíduo
best_fitness  = minc(1);
best          = par(1,:);
convergencia  = [best_fitness];

iga = 1;
while (iga < maxit)    
    iga = iga + 1; % increments generation counter
    
    %_______________________________________________________
    % Pair and mate
    mating    = ceil((popsize-keep)/2); % number of matings
    prob = flipud([1:keep]'/sum([1:keep])); % weights
    % chromosomes
    odds  = [0 cumsum(prob(1:keep))']; % probability distribution function
    pick1 = rand(1, mating); % mate #1
    pick2 = rand(1, mating); % mate #2
    % ma and pa contain the indicies of the chromosomes that will mate
    ic=1;
    while ic<=mating
        for id=2:keep+1
            if pick1(ic)<=odds(id) & pick1(ic)>odds(id-1)
                ma(ic)=id-1;
            end
            if pick2(ic)<=odds(id) & pick2(ic)>odds(id-1)
                pa(ic)=id-1;
            end
        end
        ic=ic+1;
    end
    
    %_______________________________________________________
    % Performs mating using single point crossover
    ix=1:2:keep; % index of mate #1
    xp=ceil(rand(1,mating)*Nt); % crossover point
    r=rand(1,mating); % mixing parameter
    for ic=1:mating
        xy=par(ma(ic),xp(ic))-par(pa(ic),xp(ic)); % ma and pa
        % mate
        par(keep+ix(ic),:)=par(ma(ic),:);   % 1st offspring
        par(keep+ix(ic)+1,:)=par(pa(ic),:); % 2nd offspring
        par(keep+ix(ic),xp(ic))= par(ma(ic),xp(ic))-r(ic).*xy;
        % 1st
        par(keep+ix(ic)+1,xp(ic))=par(pa(ic),xp(ic))+r(ic).*xy;
        % 2nd
        if xp(ic)<npar % crossover when last variable not selected
            par(keep+ix(ic),:)=[par(keep+ix(ic),1:xp(ic)) par(keep+ix(ic)+1,xp(ic)+1:npar)];
            par(keep+ix(ic)+1,:)=[par(keep+ix(ic)+1,1:xp(ic)) par(keep+ix(ic),xp(ic)+1:npar)];
        end
    end
    
    %_______________________________________________________
    % Mutate the population
    mrow=sort(ceil(rand(1,nmut)*(popsize-1))+1);
    mcol=ceil(rand(1,nmut)*Nt);
    for ii=1:nmut
        par(mrow(ii),mcol(ii))=(varhi-varlo)*rand+varlo; % mutation
    end
    
    %_______________________________________________________
    % The new offspring and mutated chromosomes are evaluated
    %cost = calcula_fitness(par, funcao);
    for i=1:popsize,
        cost(i) = funcao(par(i,:));   
    end;    
    
    %_______________________________________________________
    % Sort the costs and associated parameters
    [cost,ind] = sort(cost);
    par        = par(ind,:);
    
    best_fitness = min(cost);
    best         = par(1,:);
    convergencia = [convergencia best_fitness];
end %fim do laço