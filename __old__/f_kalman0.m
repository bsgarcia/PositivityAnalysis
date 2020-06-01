function fx = f_kalman0(x,P,u,in)
a = u(1);
y = u(2);
a2 = exp(P(1));
s2 = exp(P(2));
m0 = x(1:2);
S0 = x(3:4);

% prediction
mp = m0;
Sp = S0 + a2;

% update
m = NaN(2,1);
S = NaN(2,1);
if a==1
    S(1) = s2./((s2./(S0(1)+a2))+1); 
    m(1) = m0(1) + (S(1)./(S0(1)+a2)).*(y-m0(1));
    S(2) = Sp(2);
    m(2) = mp(2);
else
    S(2) = s2./((s2./(S0(2)+a2))+1); 
    m(2) = m0(2) + (S(2)./(S0(2)+a2)).*(y-m0(2));
    S(1) = Sp(1);
    m(1) = mp(1);
end

% wrap-up
fx = [m;S];

if isweird(fx)
    keyboard
end

