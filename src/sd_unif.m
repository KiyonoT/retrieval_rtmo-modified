function sd = sd_unif(lb, ub, mu)
    sd = sqrt((lb^2 + lb*ub + ub^2)/3 - mu*(lb+ub) + mu^2);
end