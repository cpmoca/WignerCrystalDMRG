function plot_density(Ne,Norb, rs, kappa, p)


xmin = -6;
xmax = 6;


%Norb =3; %number of states kept at each positions.

N=Ne*Norb;

x = linspace(xmin, xmax, 101);
x_int=linspace(xmin, xmax, 1001);

dirname = sprintf('Ne_%d_Norb_%d_rs_%.3f_kappa_%.3f_p_%.3f',Ne, Norb, rs, kappa, p);

fprintf('reading data from %s\n', dirname);
%kappa1=-2;
filename = sprintf('%s/classical_positions_MC_Ne_%d_kappa_%.3f_rs_%.3f_p_%.3f.dat',dirname, Ne,kappa, rs,p );
%filename = sprintf('%s/classical_positions_MC_Ne_%d_kappa_%.3f_rs_%.3f_p_%.3f.dat',dirname, Ne(kk),kappa1, rs(kk),p(kk) );

clasic_pos = load(filename);
%ld = min(diff(clasic_pos))/2;
ld =0.6;
%clasic_pos = linspace(-3,3,9);
error = 1e-5;
intervals = 401;
x = linspace(-15, 15, intervals);
dx = (x(end)-x(1))/(intervals-1);

HO_WF=@(n, z, x0) interp1(x, 1./sqrt(2^(n-1)*factorial(n-1))*(1/(pi*ld^2))^0.25*exp(-(x-x0).^2/(2*ld^2)).*hermiteH(n-1,(x-x0)/ld), z, 'spline');

n=1;

WF = zeros(N, numel(x));

n=1;
for k=1:Norb
    for j=1:Ne
        WF(n,:)= HO_WF(k, x, clasic_pos(j));
        n=n+1;
    end
end

%Now we construct an orthogonal set using the Gram-Schmidt transformation;

WF = WF';
Q = zeros(numel(x), N);

[m,n] = size(WF);
% compute QR using Gram-Schmidt
for j = 1:n
    v = WF(:,j);
    for i=1:j-1
        R(i,j) = Q(:,i)'*WF(:,j);
        v = v - R(i,j)*Q(:,i);
    end
    R(j,j) = norm(v);
    Q(:,j) = v/R(j,j);
end

WF=WF';
WF1=Q';


%proper normalization;

for j=1:N
    norm1 = trapz(x, WF1(j,:).^2);
    WF1(j,:)=WF1(j,:)/sqrt(norm1);
end

psi=@(n, z) interp1(x, WF1(n,:), z, 'spline');



for i=1:N
    WaveFunc(:,i) = psi(i,x);
end


for i=1:N
    WaveFunc_int(:,i) = psi(i,x_int);
end

Eigval = 1;



rho_ind = build_1pt_density([dirname, '/'], Eigval);



rho_x = zeros(1, numel(x));


for i=1:N
    for j =1:N
        for k=1:numel(x)
            rho_x(k) = rho_x(k)+WaveFunc(k, i)*rho_ind(i,j).*WaveFunc(k, j)+...
                WaveFunc(k, i)*rho_ind(i+N,j+N).*WaveFunc(k, j);
        end
        
    end
end

rho_x= 0.5*(rho_x+flip(rho_x));


if 1
    
    % plotting the data
    
    
    lw = 1.5;
    fsize = 14;
    
    h=figure('color','white','units','inches','position',[1 1 8 6]);
    
    
    h1 = plot(x,rho_x,'b-','LineWidth',lw);
    hold on;
    
    plot(clasic_pos ,zeros(1, Ne),  'ro','LineWidth',1,...
                 'MarkerEdgeColor','k',...
                 'MarkerFaceColor','g',...
                 'MarkerSize',10);
     
    
    xlim([xmin, xmax]);
    %xlim([-1, 1]);    %ylim([0 6]);
    %xlim([1e-9, energy_max]);    %ylim([0 6]);
    %ylim([0 0.01]);
    xlabel(sprintf('$x/l_d$'),'Interpreter','latex','FontSize',30);
    ylabel(sprintf('$\\rho_{%d}(x)$', Ne), 'Interpreter', 'latex', 'FontSize',30);
    
    yl=get(gca,'ylim');
    
    xl=get(gca,'xlim');
    
    % text(0.3*xl(2), 0.8*yl(2),sprintf('$\\rho_{%d} -\\rho_{%d}$', Ne2, Ne1), 'Interpreter', 'latex', 'FontSize',15);
    text(xl(1), 1.02*yl(2),sprintf('$N_{orb}=%d\\;r_s=%.2f\\;\\kappa=%.2f\\; p = %.2f$',Norb, rs, kappa, p),'Interpreter', 'latex','FontSize',20);
    
    
    yl=get(gca,'ylim');
    
    xl=get(gca,'xlim');
    
    %             text(0.6*xl(2), 0.8*yl(2),sprintf('N_e=%d', Ne(kk)),'FontSize',20);
    %             text(0.6*xl(2), 0.7*yl(2),sprintf('r_s=%.3f', rs(kk)),'FontSize',20);
    %
    
    set(gca,'FontSize',fsize);
    set(gcf, 'PaperPositionMode', 'auto');
    fname1 = sprintf('density_Ne_%d_kappa_%.3f_rs_%.3f_p_%.3f.jpg', Ne,kappa, rs, p );
    
    
    set(gcf,'paperunits','in');
    set(gcf,'papersize',[8.2,6.2]) % Desired outer dimensions
    % of figure
    
    
    if 1
        hfig = gcf;
        
        print(hfig,'-djpeg',fname1);
        
        
    end
    
    
    
end

end



