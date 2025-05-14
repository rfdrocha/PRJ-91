function [lambda_i] = calculate_lambda_i_V2(lambda_i0,mu,lambda_c)
erro_max = 10^(-6);
erro = 1;
lambda_i = lambda_i0; % chute inicial sendo lambda_i0

if lambda_c < 0
    lambda_c = 0;
end

% Funcao para encontrar lambda_i de forma iterativa
    while erro > erro_max
        new_lambda_i = lambda_i0^2/(sqrt(mu^2 + (lambda_c + lambda_i)^2));
        erro = new_lambda_i - lambda_i;
        lambda_i = new_lambda_i;
    end
end

