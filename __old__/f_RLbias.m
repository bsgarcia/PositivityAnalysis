function fx = f_RLbias(x,P,u,in)
a = u(1);
y = u(2);
ap = sig(P(1));
an = sig(P(2));
Q0 = x(1:2);

% update
Q = NaN(2,1);
if a==1
    PE = y-Q0(1);
    if PE >=0
        Q(1) = Q0(1) + ap.*(y-Q0(1));
        Q(2) = Q0(2);
    else
        Q(1) = Q0(1) + an.*(y-Q0(1));
        Q(2) = Q0(2);
    end
else
    PE = y-Q0(2);
    if PE >=0
        Q(2) = Q0(2) + ap.*(y-Q0(2));
        Q(1) = Q0(1);
    else
        Q(2) = Q0(2) + an.*(y-Q0(2));
        Q(1) = Q0(1);
    end
end

% wrap-up
fx = Q;

