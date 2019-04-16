Ne = 2;
Norb = 2;

N=Ne*Norb;%total number of orbitals.

rs = 100.0000;
kappa = 0;
p = 0;
alpha = 0.3162;
X=500; 


nfft = 2^15;%number of points in the fft; 
dxn = X/nfft; 
Fs = 1/dxn; 
df = Fs/nfft;

xn = linspace(-X/2, X/2-dxn,nfft);
fn = linspace( - Fs/2 , Fs/2 - df , nfft );
omegan = 2*pi*fn;

%Fourier transform of the kernel ;
K_x= 1./sqrt(xn.^2+alpha^2);
K_fft = 1/sqrt(2*pi)*fftshift(fft(ifftshift(K_x)))*dxn;

filename = sprintf('classical_positions_MC_Ne_%d_kappa_%.3f_rs_%.3f_p_%.3f.dat',Ne,kappa, rs,p )
clasic_pos = load(filename);


if p==0
%clasic_pos((end+1)/2)=0.0;
end

ld =1; 

error = 1e-5;
intervals = 20000; 
x = linspace(-X/2, X/2, intervals);
x1 = linspace(-X/2, X/2, intervals-1);
x2 = linspace(-X/2, X/2, intervals-2);
dx = (x(end)-x(1))/(intervals-1);
dx2=dx*dx;

V_pot = 1/4*x.^4 +1/2*(kappa)*x.^2-p*x;
%V_pot = 1/2*(kappa)*x.^2-p*x;

V_min = min(V_pot);
V_pot = V_pot-V_min;

HO_WF=@(n, z, x0) interp1(x, 1./sqrt(2^(n-1)*factorial(n-1))*(1/(pi*ld^2))^0.25*exp(-(x-x0).^2/(2*ld^2)).*hermiteH(n-1,(x-x0)/ld), z, 'spline');
n=1;

WF = zeros(N, numel(x));



n=1;
for k=1:Norb
  for j=1:Ne
        WF(n,:)=HO_WF(k, x, clasic_pos(j))  ;
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
%first derivative
psi1=@(n,z) interp1(x1, diff(WF1(n,:))./dx, z, 'spline');
%second derivative
psi2=@(n,z) interp1(x2, diff(diff(WF1(n,:)))./dx2, z, 'spline');



%two body Coulomb integrals.
V = zeros(N, N, N, N);

for i=1:N
    tic;
    fprintf('i=%d\n', i);
    for j=1:N
        for k=1:N
            for l=1:N
               % tic
                psi_ij = psi(i, xn).*psi(j, xn);
                psi_ij_fft = 1/sqrt(2*pi)*fftshift(fft(ifftshift(psi_ij)))*dxn;
                psi_ij_fft = fliplr(psi_ij_fft);
                
                psi_kl = psi(k, xn).*psi(l, xn);
                psi_kl_fft = 1/sqrt(2*pi)*fftshift(fft(ifftshift(psi_kl)))*dxn;
               
             
                V(i, j, k, l) = V(i, j, k, l)+trapz(omegan, psi_ij_fft.*K_fft.*psi_kl_fft);
                V(i,j,k,l) = V(i,j,k,l)*sqrt(2*pi);
              %  toc
              %  tic
                
%                 integrand = @(z1, z2)  psi(i, z1).*psi(k,z2).*1./sqrt((z1-z2).^2+alpha^2).*psi(j, z1).*psi(l,z2);
%                 V1 =quad2d(integrand, -15, 15, -15, 15);
%                 
%               
%                 toc
%                 fprintf('V(i,j,k,l)=%.3f\t%.3f\n',V(i,j,k,l), V1);
%                 keyboard;
            end             
        end
    end
    toc;
end

%on site energy corresponding to diagonal term;
t=zeros(N);
t0=zeros(N);
t1 = zeros(N);
t2 = zeros(N);
t4 = zeros(N);

for i=1:N
    for j=1:N
       % integrand=@(z) psi(i, z).*psi2(j,z);
       % int = quad(integrand, -X/2, X/2);
        t0(i, j)= trapz(x,psi(i, x).*psi2(j,x) );
    end
end
%how compute the matrix elements of the z

for i=1:N
    for j=1:N
        %integrand=@(z) psi(i, z).*(z).*psi(j,z);
        %int = quad(integrand, -X/2, X/2);
        t1(i, j)= trapz(x, psi(i, x).*(x).*psi(j,x));
    end
end


for i=1:N
    for j=1:N
        %integrand=@(z) psi(i, z).*(z.^2).*psi(j,z);
        %int = quad(integrand, -X/2, X/2);
        t2(i, j)= trapz(x,psi(i, x).*(x.^2).*psi(j,x) );
    end
end


for i=1:N
    for j=1:N
        %integrand=@(z) psi(i, z).*(z.^4).*psi(j,z);
        %int = quad(integrand, -X/2, X/2);
        t4(i, j)= trapz(x,psi(i, x).*(x.^4).*psi(j,x));
    end
end


%t= -0.5*t0-p*t1+0.5*(kappa)*t2;
t= -0.5*t0-p*t1+0.5*(kappa)*t2+0.25*t4;

t=0.5*(t+t');

%here we save the data to the file.

filename = sprintf('FCIDUMP.dat');
file = fopen(filename, 'w');
fprintf(file, '&FCI NORB= %d,NELEC= %d,MS2= 0,\n',2*N, Ne);
fprintf(file, 'ORBSYM=%s\n', repmat('1,', 1, 2*N));
fprintf(file, 'ISYM=1\n');
fprintf(file, '&END\n');

for i=1:N
    for j=1:N
        for k=1:N
            for l=1:N
                if abs(real(rs*V(i, j, k, l)))>error
                    fprintf(file,'%10.10f\t%d\t%d\t%d\t%d\n', real(rs*V(i, j, k, l)), i, j,k,l );
                end
            end
        end
    end
end

for i=1:N
    for j=1:N
        for k=1:N
            for l=1:N
                if abs(real(rs*V(i, j, k, l)))>error
                    fprintf(file,'%10.10f\t%d\t%d\t%d\t%d\n', real(rs*V(i, j, k, l)), i, j,k+N,l+N );
                end
            end
        end
    end
end



for i=1:N
    for j=1:N
        for k=1:N
            for l=1:N
                if abs(real(rs*V(i, j, k, l)))>error
                    fprintf(file,'%10.10f\t%d\t%d\t%d\t%d\n',real(rs*V(i, j, k, l)), i+N, j+N,k,l );
                end
            end
        end
    end
end

for i=1:N
    for j=1:N
        for k=1:N
            for l=1:N
                if abs(rs*V(i, j, k, l))>error
                    fprintf(file,'%10.10f\t%d\t%d\t%d\t%d\n', rs*V(i, j, k, l), i+N, j+N,k+N,l+N );
                end
            end
        end
    end
end



for i=1:N
    for j=1:N
        if abs(t(i, j))>error
            fprintf(file,'%10.10f\t%d\t%d\t%d\t%d\n',t(i,j), i, j, 0, 0 );
        end
    end
end

for i=1:N
    for j=1:N
        if abs(t(i, j))>error
            fprintf(file,'%10.10f\t%d\t%d\t%d\t%d\n',t(i,j), i+N, j+N, 0, 0 );
        end
    end
end


fprintf(file,'%10.10f\t%d\t%d\t%d\t%d\n',0, 0, 0, 0, 0 );
fclose(file);
