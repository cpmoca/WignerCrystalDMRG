

function E = energy(Ne, x, kappa, rs, p)
%dimensionless energy;
alpha = 0.3162;
clx = sort(x);

E_potential=1/4*sum(x.^4)+1/2*kappa*sum(x.^2)-p*sum(x);

%pure Coulomb
E_Coulomb=0;


for i=1:Ne
    for j=1:Ne
        if i==j
            continue
        end
        E_Coulomb=E_Coulomb+0.5*rs*1./abs((x(i)-x(j)));
    end
end

E=E_potential+E_Coulomb;
if Ne==1 
    E=E_potential;
end 
%Elliptic Coulomb;

%E=sum(y.^2)+sum( 2/pi*rs./sqrt(diff(y).^2+4*alpha^2).*ellipticK(4*alpha^2./(diff(y).^2+4*alpha^2)));


end
