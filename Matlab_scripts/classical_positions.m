%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates classical position of electrons %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ne = 2;       % Number of electrons
rs = 100.0000;% Strenght of the interaction
kappa = 0;    % Strength of the quadratix confinement
p = 0;        % External electric field 


L =10;        %in units of ld the oscillator length for x^4 potential.
steps = 10^6; %Number of Monte-Carlo steps

x0 = linspace(-L/2, L/2, Ne);
x0(2:Ne-1)= x0(2:Ne-1)+0.1*(rand-0.5)*L/(Ne);
x1=x0;

filename = sprintf('classical_positions_MC_Ne_%d_kappa_%.3f_rs_%.3f_p_%.3f.dat',Ne,kappa, rs,p );
file = fopen(filename, 'w');

%if kappa <-3 
%  kappa=-3;%to fix kappa; to -3
%end

for i=1:steps
    
    %now pick a random electron
    
    el_no = randi(Ne);
    
    % move it randomly
    
    if Ne>1 & el_no>1 & el_no<Ne
        
        x1(el_no)=x0(el_no)+(rand-0.5)*(L/(2*Ne));
        
        while x1(el_no) > x0(el_no+1) | x1(el_no)<x1(el_no-1)
            x1(el_no)=x0(el_no)+0.1*rand*(L*Ne);
        end
    end
    if el_no==Ne & Ne>1
        x1(Ne)=x0(Ne)+(rand-0.5)*(L/(2*Ne));
        while x1(Ne) < x0(Ne-1) | x1(Ne)>L/2
            x1(Ne)=x0(Ne)+0.1*(rand-0.5)*(L/(2*Ne));
        end
    end
    
    if el_no==1 & Ne>1
        x1(1)=x0(1)+(rand-0.5)*(L/(2*Ne));
        while x1(1) > x0(2) | x1(1)<-L/2
            x1(1)=x0(1)+0.1*(rand-0.5)*(L/(2*Ne));
        end
    end
    
    if Ne==1
        x1(1)=x0(1)+(rand-0.5)*(L/(2*Ne));
        while x1(1) > L/2 | x1(1)<-L/2
            x1(1)=x0(1)+0.5*(rand-0.5)*(L/(2*Ne));
        end
    end
    
    
    
    x1=sort (x1);
    dE = energy(Ne, x1, kappa, rs, p)-energy(Ne, x0, kappa, rs, p);
    
    if dE<0
        x0=x1;
    end
    
    x1(1:Ne)=-x1(1:Ne);
    x1=sort (x1);
    dE = energy(Ne, x1, kappa, rs, p)-energy(Ne, x0, kappa, rs, p);
    if dE<0
        x0=x1;
    end
    
    
    
    
end

for j=1:Ne
  fprintf(file, '%.3f\t', x0(j));
end


%if Ne==1
%    fprintf(file, '%.3f\n',  x0(1));
%end
%
%
%if Ne==2
%    fprintf(file, '%.3f\t%3f\n',  x0(1), x02));
%end
%
%if Ne==3
%    fprintf(file, '%.3f\t%.3f\t%.3f\n', x0(1), x0(2), x0(3));
%end
%
%if Ne==5
%    fprintf(file, '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', x0(1), x0(2), x0(3), x0(4), x0(5));
%end
%
%if Ne==7
%    fprintf(file, '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n', x0(1), x0(2), x0(3), x0(4), x0(5),x0(6), x0(7) );
%end





y= linspace(-L/2, L/2, 400);
V= 1/4*y.^4+1/2*kappa*y.^2-p*y;
V=V-min(V);

%
%             fprintf('i=%d\n',i);
%             V1=interp1(y, V, x0);
%             plot(y, V,'b-');
%             hold on
%             plot(x0 ,V1,  'ro','LineWidth',1,...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','g',...
%                 'MarkerSize',10);
%             ylim([0, 20]);
%             hold off;
%
%
%             keyboard;


fclose(file);








