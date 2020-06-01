function gx = g_kalman0(x,P,u,in)
m0 = x(1:2);

gx = sig((m0(1)-m0(2))./exp(P(1)));

