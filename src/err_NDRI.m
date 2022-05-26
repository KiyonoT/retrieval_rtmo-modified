function result = err_NDRI(r1, r2, s1, s2)
    NDRI = (r1-r2)/(r1+r2);
    err2 = (2*r2/(r1+r2)^2)^2*s1^2 + (2*r1/(r1+r2)^2)^2*s2^2;
    result = [NDRI, sqrt(err2)];