function convergencia = DA(dim, SearchAgents_no, Max_iteration, fobj)
% Definindo os limites 
[lb, ub] = get_espaco_busca(fobj);
fobj = str2func(fobj);
%Limite Inferior
lb = lb*ones(1, dim);
%Limite Superior
ub = ub*ones(1, dim);

cg_curve = zeros(1, Max_iteration);
if size(ub,2) == 1
    ub = ones(1, dim) * ub;
    lb = ones(1, dim) * lb;
end

%The initial radius of gragonflies' neighbourhoods
r         = (ub - lb) / 10;
Delta_max = (ub-lb) / 10;

Food_fitness = inf;
Food_pos     = zeros(dim,1);

Enemy_fitness = -inf;
Enemy_pos     = zeros(dim,1);

X       = inicializacao(SearchAgents_no, dim, ub, lb);
Fitness = zeros(1, SearchAgents_no);

DeltaX  = inicializacao(SearchAgents_no, dim, ub, lb);
convergencia = [];
for iter=1:Max_iteration   
    r = (ub-lb)/4+((ub-lb)*(iter/Max_iteration)*2);
    w = 0.9 - iter * ((0.9 - 0.2)/Max_iteration);
       
%     my_c = 0.1 - iter*((0.1-0)/(Max_iteration/2));
%     if my_c < 0
%        my_c=0;
%     end
%     
%     s = 2*rand*my_c; % Seperation weight
%     a = 2*rand*my_c; % Alignment weight
%     c = 2*rand*my_c; % Cohesion weight
%     f = 2*rand;      % Food attraction weight
%     e = my_c;        % Enemy distraction weight
    
    s = 0.1; % Seperation weight
    a = 0.1; % Alignment weight
    c = 0.7; % Cohesion weight
    f = 1;   % Food attraction weight
    e = 1;   % Enemy distraction weight
    
    for i=1:SearchAgents_no %Calculate all the objective values first
        Fitness(1,i)=fobj(X(:,i)');
        if Fitness(1,i)<Food_fitness
            Food_fitness=Fitness(1,i);
            Food_pos=X(:,i);
        end
        
        if Fitness(1,i) > Enemy_fitness
            if all(X(:,i)<ub') && all( X(:,i)>lb')
                Enemy_fitness=Fitness(1,i);
                Enemy_pos=X(:,i);
            end
        end
    end
    
    for i=1:SearchAgents_no
        index=0;
        neighbours_no=0;
        
        clear Neighbours_DeltaX
        clear Neighbours_X
        %find the neighbouring solutions
        for j=1:SearchAgents_no
            Dist2Enemy=distance(X(:,i),X(:,j));
            if (all(Dist2Enemy<=r) && all(Dist2Enemy~=0))
                index=index+1;
                neighbours_no=neighbours_no+1;
                Neighbours_DeltaX(:,index)=DeltaX(:,j);
                Neighbours_X(:,index)=X(:,j);
            end
        end
        
        % Seperation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.1)
        S=zeros(dim,1);
        if neighbours_no>1
            for k=1:neighbours_no
                S=S+(Neighbours_X(:,k)-X(:,i));
            end
            S=-S;
        else
            S=zeros(dim,1);
        end
        
        % Alignment%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.2)
        if neighbours_no>1
            A=(sum(Neighbours_DeltaX')')/neighbours_no;
        else
            A=DeltaX(:,i);
        end
        
        % Cohesion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.3)
        if neighbours_no>1
            C_temp=(sum(Neighbours_X')')/neighbours_no;
        else
            C_temp=X(:,i);
        end
        
        C=C_temp-X(:,i);
        
        % Attraction to food%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.4)
        Dist2Food=distance(X(:,i),Food_pos(:,1));
        if all(Dist2Food<=r)
            F=Food_pos-X(:,i);
        else
            F=0;
        end
        
        % Distraction from enemy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Eq. (3.5)
        Dist2Enemy=distance(X(:,i),Enemy_pos(:,1));
        if all(Dist2Enemy<=r)
            Enemy=Enemy_pos+X(:,i);
        else
            Enemy=zeros(dim,1);
        end
        
        for tt=1:dim
            if X(tt,i)>ub(tt)
                X(tt,i)=lb(tt);
                DeltaX(tt,i)=rand;
            end
            if X(tt,i)<lb(tt)
                X(tt,i)=ub(tt);
                DeltaX(tt,i)=rand;
            end
        end
        
        if any(Dist2Food>r)
            if neighbours_no>1
                for j=1:dim
                    DeltaX(j,i)=w*DeltaX(j,i)+rand*A(j,1)+rand*C(j,1)+rand*S(j,1);
                    if DeltaX(j,i)>Delta_max(j)
                        DeltaX(j,i)=Delta_max(j);
                    end
                    if DeltaX(j,i)<-Delta_max(j)
                        DeltaX(j,i)=-Delta_max(j);
                    end
                    X(j,i)=X(j,i)+DeltaX(j,i);
                end
            else
                % Eq. (3.8)
                X(:,i)=X(:,i)+Levy(dim)'.*X(:,i);
                DeltaX(:,i)=0;
            end
        else
            for j=1:dim
                % Eq. (3.6)
                DeltaX(j,i)=(a*A(j,1)+c*C(j,1)+s*S(j,1)+f*F(j,1)+e*Enemy(j,1)) + w*DeltaX(j,i);
                if DeltaX(j,i)>Delta_max(j)
                    DeltaX(j,i)=Delta_max(j);
                end
                if DeltaX(j,i)<-Delta_max(j)
                    DeltaX(j,i)=-Delta_max(j);
                end
                X(j,i)=X(j,i)+DeltaX(j,i);
            end 
        end
        
        Flag4ub=X(:,i)>ub';
        Flag4lb=X(:,i)<lb';
        X(:,i)=(X(:,i).*(~(Flag4ub+Flag4lb)))+ub'.*Flag4ub+lb'.*Flag4lb;
    end
    convergencia = [convergencia Food_fitness];
    %Best_score=Food_fitness;
    %Best_pos=Food_pos;
    
    %cg_curve(iter)=Best_score;
    end
end

function pop_inicial = inicializacao(SearchAgents_no, dim, ub, lb)

    Boundary_no= size(ub, 2); % numnber of boundaries

    % If the boundaries of all variables are equal and user enter a signle
    % number for both ub and lb
    if Boundary_no==1
        ub_new=ones(1,dim)*ub;
        lb_new=ones(1,dim)*lb;
    else
         ub_new=ub;
         lb_new=lb;   
    end

    % If each variable has a different lb and ub
        for i=1:dim
            ub_i=ub_new(i);
            lb_i=lb_new(i);
            Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
        end

    pop_inicial = Positions';
end

function o = distance(a,b)
    for i=1:size(a,1)
        o(1,i)=sqrt((a(i)-b(i))^2);
    end
end
