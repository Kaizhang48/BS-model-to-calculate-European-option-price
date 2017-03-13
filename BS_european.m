function P=BS_european(id,K,T,S,sigma,q,r)
        d1=(log(S./K)+T*(r-q+power(sigma,2)/2))/(sigma*power(T,1/2));%d+
        d2=d1-sigma*power(T,1/2);%d-
switch id
    case 1
        P=S.*exp(-q.*T).*normcdf(d1,0,1)-K.*exp(-r.*T).*normcdf(d2,0,1);
    case -1
        P=K.*exp(-r.*T).*normcdf(-d2,0,1)-S.*exp(-q.*T).*normcdf(-d1,0,1);
    otherwise
        error('id is not satisfied')
end
end