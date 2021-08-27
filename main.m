funcao = 'fun_ackley';
format shortE;

gera_grafico      = 1;
conj_experimentos = 2; %1: 10 dimensões(configurações do paper CPO); 
                       %2: 30 dimensões(configurações do paper GWOCSA);
                       
%%%% OBSERVAR ESTE PARÂMETRO %%%
teste             = 0; %0: roda 30 vezes
                       %1: roda 3 vezes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if teste == 1 
    num_rodadas = 3;
else
    num_rodadas = 30;
end

if conj_experimentos == 1
   dim           = 10;   
   populacao     = 20;   
   max_iter      = 2000;
   gera_grafico  = 0;
   
else 
   dim       = 30;   
   populacao = 30;   
   max_iter  = 300;
   gera_grafico = 1;
end

dim = verifica_funcao_dim_fixa(funcao, dim);

%%%%%% Definição dos vetores %%%%%%%
converg_CSA                    = [];
hist_converg_CSA               = [];
hist_best_fitness_CSA          = [];

converg_PSO                    = [];
hist_converg_PSO               = [];
hist_best_fitness_PSO          = [];

converg_BA                     = [];
hist_converg_BA                = [];
hist_best_fitness_BA           = [];

converg_GA                     = [];
hist_converg_GA                = [];
hist_best_fitness_GA           = [];

converg_DA                     = [];
hist_converg_DA                = [];
hist_best_fitness_DA           = [];

converg_GWO                    = [];
hist_converg_GWO               = [];
hist_best_fitness_GWO          = [];

converg_CSA_Proposto           = [];
hist_converg_CSA_Proposto      = [];
hist_best_fitness_CSA_Proposto = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear CSA;
clear PSO;
clear GWO;
clear BA;
clear GA;
clear DA;
clear CSA_Proposto;

for rodada = 1:num_rodadas
    % CSA retorna um vetor com "max_iter" elementos. Cada elemento representa o 
    % melhor valor de fitness da iteração i das "max_iter" iterações.
    converg_CSA                = CSA(dim, populacao, max_iter, funcao);
    best_fitness_CSA           = converg_CSA(max_iter);
    hist_best_fitness_CSA      = [hist_best_fitness_CSA best_fitness_CSA];
    % Monta uma matriz de ordem "rodadas" x "max_iter"
    hist_converg_CSA           = [hist_converg_CSA; converg_CSA];
    
    converg_PSO                = PSO(dim, populacao, max_iter, funcao);
    best_fitness_PSO           = converg_PSO(max_iter);
    hist_best_fitness_PSO      = [hist_best_fitness_PSO best_fitness_PSO];
    hist_converg_PSO           = [hist_converg_PSO; converg_PSO];
    
    converg_BA                 = BA(dim, populacao, max_iter, funcao);
    best_fitness_BA            = converg_BA(max_iter);
    hist_best_fitness_BA       = [hist_best_fitness_BA best_fitness_BA];
    hist_converg_BA            = [hist_converg_BA; converg_BA];
    
    converg_GA                 = GA(dim, populacao, max_iter, funcao);
    best_fitness_GA            = converg_GA(max_iter);
    hist_best_fitness_GA       = [hist_best_fitness_GA best_fitness_GA];
    hist_converg_GA            = [hist_converg_GA; converg_GA];
    
    converg_DA                 = DA(dim, populacao, max_iter, funcao);
    best_fitness_DA            = converg_DA(max_iter);
    hist_best_fitness_DA       = [hist_best_fitness_DA best_fitness_DA];
    hist_converg_DA            = [hist_converg_DA; converg_DA];
    
    converg_GWO                = GWO(dim, populacao, max_iter, funcao);
    best_fitness_GWO           = converg_GWO(max_iter);
    hist_best_fitness_GWO      = [hist_best_fitness_GWO best_fitness_GWO];
    hist_converg_GWO           = [hist_converg_GWO; converg_GWO];
    
    converg_CSA_Proposto           = CSA_Proposto(dim, populacao, max_iter, funcao);
    best_fitness_CSA_Proposto      = converg_CSA_Proposto(max_iter);
    hist_best_fitness_CSA_Proposto = [hist_best_fitness_CSA_Proposto best_fitness_CSA_Proposto];
    hist_converg_CSA_Proposto      = [hist_converg_CSA_Proposto; converg_CSA_Proposto];
end;

%%% Calcula médias
media_converg_CSA          = mean(hist_converg_CSA); 
media_converg_PSO          = mean(hist_converg_PSO);
media_converg_BA           = mean(hist_converg_BA); 
media_converg_GA           = mean(hist_converg_GA); 
media_converg_DA           = mean(hist_converg_DA);
media_converg_GWO          = mean(hist_converg_GWO);
media_converg_CSA_Proposto = mean(hist_converg_CSA_Proposto);
%%%

disp(['CSA.....: ', num2str(mean(hist_best_fitness_CSA)), ' (', num2str(std(hist_best_fitness_CSA)), ')']);
disp(['PSO.....: ', num2str(mean(hist_best_fitness_PSO)), ' (', num2str(std(hist_best_fitness_PSO)), ')']);
disp(['BA......: ', num2str(mean(hist_best_fitness_BA)), ' (', num2str(std(hist_best_fitness_BA)), ')']);
disp(['GA......: ', num2str(mean(hist_best_fitness_GA)), ' (', num2str(std(hist_best_fitness_GA)), ')']);
disp(['DA......: ', num2str(mean(hist_best_fitness_DA)), ' (', num2str(std(hist_best_fitness_DA)), ')']);
disp(['GWO.....: ', num2str(mean(hist_best_fitness_GWO)), ' (', num2str(std(hist_best_fitness_GWO)), ')']);
disp(['Proposta: ', num2str(mean(hist_best_fitness_CSA_Proposto)), ' (', num2str(std(hist_best_fitness_CSA_Proposto)), ')']);


% Teste: Proposta x CSA
[pval, H] = signrank(hist_best_fitness_CSA_Proposto, hist_best_fitness_CSA, 'tail','both', 'method','exact','alpha', 0.05);
if H == 1
    if mean(hist_best_fitness_CSA_Proposto) < mean(hist_best_fitness_CSA)
        disp(['Proposta supera CSA com significância ', num2str(pval)]);
    elseif mean(hist_best_fitness_CSA_Proposto) > mean(hist_best_fitness_CSA)
        disp(['Proposta perde para CSA com significância ', num2str(pval)]);    
    end
else
    disp('Não há significância estatística entre CSA-Proposta x CSA.');
end

% Teste: Proposta x PSO
[pval, H] = signrank(hist_best_fitness_CSA_Proposto, hist_best_fitness_PSO, 'tail','both', 'method','exact','alpha', 0.05);
if H == 1
    if mean(hist_best_fitness_CSA_Proposto) < mean(hist_best_fitness_PSO)
        disp(['Proposta supera PSO com significância ', num2str(pval)]);
    elseif mean(hist_best_fitness_CSA_Proposto) > mean(hist_best_fitness_PSO)
        disp(['Proposta perde para PSO com significância ', num2str(pval)]);    
    end
else
    disp('Não há significância estatística entre CSA-Proposta x PSO.');
end

% Teste: Proposta x BA
[pval, H] = signrank(hist_best_fitness_CSA_Proposto, hist_best_fitness_BA, 'tail','both', 'method','exact','alpha', 0.05);
if H == 1
    if mean(hist_best_fitness_CSA_Proposto) < mean(hist_best_fitness_BA)
        disp(['Proposta supera BA com significância ', num2str(pval)]);
    elseif mean(hist_best_fitness_CSA_Proposto) > mean(hist_best_fitness_BA)
        disp(['Proposta perde para BA com significância ', num2str(pval)]);    
    end
else
    disp('Não há significância estatística entre CSA-Proposta x BA.');
end

% Teste: Proposta x GA
[pval, H] = signrank(hist_best_fitness_CSA_Proposto, hist_best_fitness_GA, 'tail','both', 'method','exact','alpha', 0.05);
if H == 1
    if mean(hist_best_fitness_CSA_Proposto) < mean(hist_best_fitness_GA)
        disp(['Proposta supera GA com significância ', num2str(pval)]);
    elseif mean(hist_best_fitness_CSA_Proposto) > mean(hist_best_fitness_GA)
        disp(['Proposta perde para GA com significância ', num2str(pval)]);    
    end
else
    disp('Não há significância estatística entre CSA-Proposta x GA.');
end

% Teste: Proposta x DA
[pval, H] = signrank(hist_best_fitness_CSA_Proposto, hist_best_fitness_DA, 'tail','both', 'method','exact','alpha', 0.05);
if H == 1
    if mean(hist_best_fitness_CSA_Proposto) < mean(hist_best_fitness_DA)
        disp(['Proposta supera DA com significância ', num2str(pval)]);
    elseif mean(hist_best_fitness_CSA_Proposto) > mean(hist_best_fitness_DA)
        disp(['Proposta perde para DA com significância ', num2str(pval)]);    
    end
else
    disp('Não há significância estatística entre CSA-Proposta x DA.');
end


% Teste: Proposta X GWO
[pval, H] = signrank(hist_best_fitness_CSA_Proposto, hist_best_fitness_GWO, 'tail','both', 'method','exact','alpha', 0.05);
if H == 1
    if mean(hist_best_fitness_CSA_Proposto) < mean(hist_best_fitness_GWO)
        disp(['Proposta supera GWO com significância ', num2str(pval)]);
    elseif mean(hist_best_fitness_CSA_Proposto) > mean(hist_best_fitness_GWO)
        disp(['Proposta perde para GWO com significância ', num2str(pval)]);    
    end
else
    disp('Não há significância estatística entre CSA-Proposta x GWO.');
end

if gera_grafico == 1
    step     = 0.04 * max_iter;
    % plot((1:step:max_iter), log(media_converg_CSA(1:step:end)),'--bp', ...                  
    %      (1:step:max_iter), log(media_converg_PSO(1:step:end)),'-go',...    
    %      (1:step:max_iter), log(media_converg_BA(1:step:end)),'-k',... 
    %      (1:step:max_iter), log(media_converg_GA(1:step:end)),'sm',... 
    %      (1:step:max_iter), log(media_converg_CSA_Proposto(1:step:end)), '-.r*','MarkerSize', 3.5);  
    % legend('CSA', 'PSO', 'BA', 'GA', 'CSA-Proposto', 'Location','northoutside', 'Orientation','horizontal');
    
    

    [qntdColum, local] = verifica_funcao_grafico(funcao);
    
    %Marcadores: https://www.mathworks.com/help/matlab/creating_plots/create-line-plot-with-markers.html
    semilogy((1:step:max_iter), media_converg_CSA(1:step:end),'--bp', ...                  
             (1:step:max_iter), media_converg_PSO(1:step:end),'-go',...    
             (1:step:max_iter), media_converg_BA(1:step:end),'-k',... 
             (1:step:max_iter), media_converg_GA(1:step:end),'-xm',...
             (1:step:max_iter), media_converg_GWO(1:step:end),'-h',...
             (1:step:max_iter), media_converg_DA(1:step:end),'-s',...
             (1:step:max_iter), media_converg_CSA_Proposto(1:step:end), '-.r*','MarkerSize', 4);  
    legend({'CSA', 'PSO', 'BA', 'GA', 'GWO', 'DA', 'Proposta'}, ...
            'Location', local,'NumColumns', qntdColum);
        
    
    ylabel('Fitness');
    xlabel('Iterações');
    set(gca,'XTick',[0:max_iter/10:max_iter]);
    curtick = get(gca, 'XTick');
    set(gca, 'XTickLabel', cellstr(num2str(curtick(:))));           

    fig_exportada = [funcao, int2str(populacao), 'corvos_', int2str(dim), 'D.png'];  
    set(gca, 'color', 'none');
    set(gcf, 'color', 'none');
    path = 'graficos\';
    export_fig([path, fig_exportada], '-png');
end;