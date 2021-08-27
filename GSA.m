function [Fbest, Lbest, BestChart, MeanChart] = GSA(N, max_it, dim, funcao)
 
ElitistCheck=1; 
Rpower=1;
min_flag=1; % 1: minimization, 0: maximization
F_index=1;

%V:   Velocity.
%a:   Acceleration.
%M:   Mass.  Ma=Mp=Mi=M;
%dim: Dimension of the test function.
%N:   Number of agents.
%X:   Position of agents. dim-by-N matrix.
%R:   Distance between agents in search space.
%[low-up]: Allowable range for search space.
%Rnorm:  Norm in eq.8.
%Rpower: Power of R in eq.7.
 
Rnorm = 2; 
 
[low, up] = get_espaco_busca(funcao);
funcao = str2func(funcao);

X = initialization(dim, N, up, low); 

% %Novo
% for i=1:N
%     fitness(i) = funcao(X(i, :));
% end 
% [best best_X] = min(fitness); %minimization.
% Fbest   = best;
% Lbest         = X(best_X,:);
% %Fim

%create the best so far chart and average fitnesses chart.
BestChart=[];
MeanChart=[];

V=zeros(N,dim);

for iteration=1:max_it
    X = space_bound(X, up, low); 

    for i=1:N
      fitness(i) = funcao(X(i, :));
    end 

    [best best_X] = min(fitness); %minimization.

    if iteration==1
      Fbest=best;
      Lbest=X(best_X,:);
    end
    
    if best < Fbest  %minimization.
        Fbest = best;
        Lbest = X(best_X,:);
    end
    
    BestChart=[BestChart Fbest];
    MeanChart=[MeanChart mean(fitness)];

    %Calculation of M. eq.14-20
    [M]=massCalculation(fitness,min_flag); 

    %Calculation of Gravitational constant. eq.13.
    G=Gconstant(iteration,max_it); 

    %Calculation of accelaration in gravitational field. eq.7-10,21.
    a=Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it);

    %Agent movement. eq.11-12
    [X,V]=move(X,a,V);

end %iteration
end

function [X]=initialization(dim, N, up, down)

    if size(up,2)==1
        X=rand(N,dim).*(up-down)+down;
    end
    if size(up,2)>1
        for i=1:dim
            high = up(i);
            low  = down(i);
            X(:,i) = rand(N,1).*(high-low)+low;            
        end
    end
end

function  X=space_bound(X,up,low);

[N,dim] = size(X);
for i=1:N 
    Tp=X(i,:)>up;Tm=X(i,:)<low;X(i,:)=(X(i,:).*(~(Tp+Tm)))+((rand(1,dim).*(up-low)+low).*(Tp+Tm));
end

end

function [M]=massCalculation(fit,min_flag);

 Fmax=max(fit); Fmin=min(fit); Fmean=mean(fit); 
 [i N]=size(fit);

if Fmax==Fmin
   M=ones(N,1);
else
    
   if min_flag==1 %for minimization
      best=Fmin;worst=Fmax; %eq.17-18.
   else %for maximization
      best=Fmax;worst=Fmin; %eq.19-20.
   end
  
   M=(fit-worst)./(best-worst); %eq.15,

end

M=M./sum(M); %eq. 16.
end

function G=Gconstant(iteration,max_it)
  alfa=20;G0=100;
  G=G0*exp(-alfa*iteration/max_it); %eq. 28.
end

function [X,V]=move(X,a,V)
    [N,dim]=size(X);
    V = rand(N,dim).*V+a; %eq. 11.
    X = X+V; %eq. 12
end

function a = Gfield(M, X, G, Rnorm, Rpower, ElitistCheck, iteration, max_it);
    [N,dim]=size(X);
    final_per=2; %In the last iteration, only 2 percent of agents apply force to the others.

    %%%%total force calculation
     if ElitistCheck==1
         kbest=final_per+(1-iteration/max_it)*(100-final_per); %kbest in eq. 21.
         kbest=round(N*kbest/100);
     else
         kbest=N; %eq.9.
     end
     [Ms ds]=sort(M,'descend');

     for i=1:N
         E(i,:)=zeros(1,dim);
         for ii=1:kbest
             j=ds(ii);
             if j~=i
                R=norm(X(i,:)-X(j,:),Rnorm); %Euclidian distanse.
             for k=1:dim 
                 E(i,k)=E(i,k)+rand*(M(j))*((X(j,k)-X(i,k))/(R^Rpower+eps));
                  %note that Mp(i)/Mi(i)=1
             end
             end
         end
     end

    %%acceleration
    a=E.*G; %note that Mp(i)/Mi(i)=1
end

