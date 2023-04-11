function dWij = getdW(rij_norm,h)

x = rij_norm;
if(0 <= x && x<=1)
    dWij = (9/4)*x^2-3*x;
elseif(1<x && x<=2)
    dWij = -0.75*(2-x)^2;
else
    dWij = 0;
end
dWij = (1/(3.14*h^4)).*dWij;

end