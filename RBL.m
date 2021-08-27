
function nova_solucao = RBL(solucao, l, u)
    nova_solucao = [];
    
    dim = length(solucao);
    for i = 1:dim
        z = solucao(1, i);
        u = z - (l + u) / 2;
        v = sqrt((z - l) * (u - z));
        r = (u - l) / 2;
        beta = 180 * (0.25 * randn + 1);
        u_star = u * cosd(beta) - v * sind(beta);
        z_star = (l + u) / 2 + u_star;
        nova_solucao = [nova_solucao z_star];
    end