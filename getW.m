function Wij = getW(rij_norm,h)

x = rij_norm;
if(0 <= x && x<=1)
    Wij = 1 - 1.5*x^2 + 0.75*x^3;
elseif(1<x && x<=2)
    Wij = 0.25*(2-x)^3;
else
    Wij = 0;
end
Wij = (1/(3.14*h^3)).*Wij;

end

