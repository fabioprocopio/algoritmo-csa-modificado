              %
function convergencia = BA(dimensoes, num_morcegos, max_iteracoes, funcao)

amplitude_max        = 0.5;
amplitude_min        = 0.5;
taxa_max             = 0.5;
taxa_min             = 0.5;
freq_min             = 0;
freq_max             = 1;
gamma                = 0.05;
alpha                = 0.95;

n = num_morcegos;
r0Max = taxa_max;
r0Min = taxa_min;
AMax = amplitude_max;
AMin = amplitude_min;
r = rand(n, 1 ) .* 0.5 + 0;
r0 = rand(n, 1 ) .* (r0Max - r0Min ) + r0Min;
A = rand(n, 1 ) .* (AMax - AMin ) + AMin;

% Definindo os limites 
[Lb, Ub] = get_espaco_busca(funcao);
%Limite Inferior
Lb=Lb*ones(1,dimensoes);
%Limite Superior
Ub=Ub*ones(1,dimensoes);

% Convertendo a função
str_funcao = funcao;
funcao = str2func(funcao);

% Initial arrays
Q = zeros(n, 1);         % Frequency
v = zeros(n, dimensoes); % Velocities

% Initialize the population/solutions
for i=1:n,
  Sol(i,:)   = Lb+(Ub-Lb).*rand(1,dimensoes);
  Fitness(i) = funcao(Sol(i,:));
end

% Find the current best
[best_fitness,I]=min(Fitness);
best=Sol(I,:);

%alpha = 0.9;
%gamma = 0.9;
iter               = 1;
convergencia = [];
while (iter <= max_iteracoes)
    for i=1:n,  
        %Q(i)   = freq_min + (freq_min-freq_max) * rand;
        Q(i)   = freq_min + (freq_max-freq_min) * rand;
        v(i,:) = v(i,:)   + (Sol(i,:) - best) * Q(i);
        S(i,:) = Sol(i,:) + v(i,:);
        
        S(i,:) = simplebounds(S(i,:), Lb, Ub);
        
        % Taxa de pulso
        if rand > r(i)
            % Gera uma solução em torno da melhor solução encontrada até então
            S(i,:) = best + 0.001 * randn(1, dimensoes);
        end

        % Evaluate new solutions
        Fnew = funcao(S(i,:));
        % If the solution improves or not too loudness
        if (rand < A(i)  && Fnew <= Fitness(i))           
            Sol(i,:)   = S(i,:);
            Fitness(i) = Fnew;        
             
            A(i) = A(i) * alpha;
            r(i) = r0(i) * ( 1 - exp( -gamma * iter ) );
        end

        % Update the current best
        if Fnew <= best_fitness,
            best         = S(i,:);
            best_fitness = Fnew;
        end
    end
       
    convergencia = [convergencia best_fitness]; 
    iter = iter + 1;
end      
    
% Application of simple limits/bounds
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound vector
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);

  % Apply the upper bound vector 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;