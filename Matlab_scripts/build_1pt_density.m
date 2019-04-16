
function [D] = build_1pt_density(Dirname, Eigval)
Aup   = load([Dirname, 'cp_uc_u_rho_1pt_eig', num2str(Eigval), 'eig', num2str(Eigval), '.mat.out'], '-mat');
Adown = load([Dirname, 'c_ucp_u_rho_1pt_eig', num2str(Eigval), 'eig', num2str(Eigval), '.mat.out'], '-mat');
Bup   = load([Dirname, 'cp_dc_d_rho_1pt_eig', num2str(Eigval), 'eig', num2str(Eigval), '.mat.out'], '-mat');
Bdown = load([Dirname, 'c_dcp_d_rho_1pt_eig', num2str(Eigval), 'eig', num2str(Eigval), '.mat.out'], '-mat');
%Adown = load([Dirname, 'c_ucp_u_rho_1pt_eig1eig1.mat.out'], '-mat')
%Bup   = load([Dirname, 'cp_dc_d_rho_1pt_eig1eig1.mat.out'], '-mat')
%Bdown = load([Dirname, 'c_dcp_d_rho_1pt_eig1eig1.mat.out'], '-mat')
%Aup.A_in
%Bup.A_in
%Adown.A_in
%Bdown.A_in
Cup   = Aup.A_in+Bup.A_in;
Cdown = Adown.A_in+Bdown.A_in;
%D=triu(C)+triu(C,1)'
D=triu(Cup)+Aup.Phasefactor*triu(Cdown,1)';
%save('D.mat', 'D');

%D-D' %this must be zero for a symmetric matrix
