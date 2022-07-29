% Two-phase: cytosol and actomyosin network
% Steady-state
% Unknowns:
% pc, vc, vn, theta_n, theta_c, cin, v0
% Velocities and equations are written in the frame of the moving cell
% Questions? Contact: Yizeng Li - liyizeng52@hotmail.com

clear
clc
close all

L0 = 200;         % (micron) channel length
b  = 3.5;           % (micron) channel width
% The above two parameters are used to calcuate the hydraulic resistence.
% dg = 12 * viscosity (Pa s) * (L0 - L) (micron) / b^2 (micron^2)
% this dg is in Pa d/micron, need to divide it by 1d3 to get Pa s/nm,
% the latter is the unit system used in this program

L = 50.d3;              % (nm) cell length
R = 8.31451;        % (J/mol K) Ideal gas constant
T = 310;            % (K) absolute temperature

p0f = 0*1d5;            % (Pa) external pressure at the front
p0b = 0*1d5;            % (Pa) external pressure at the back
c0f = 340;            % (mol/m^3 = mM) external ion concentration at the front
c0b = 340;            % (mol/m^3 = mM) external ion concentration at the back
fextf = 0d2;            % (Pa) external force per unit area at the front of the cell
fextb = 0d2;            % (Pa) external force per unit area at the back of the cell

Thetac = 0.1;           % (mM) reference value of G-actin
Thetan = 0.2;           % (mM) reference value of F-actin
Theta  = 0.3;           % (mM) reference total G- and F-actin, \int_0^L (thetan + thetac)dx/L
thetacc  = 0.2d-3;      % (mM) Critical value for actin polymerization

ksigman = 1d2;          % (Pa /mM) Coefficient of passive actin network pressure
gamma0 = 6.0d-4;          % (1/s) constant rate of actin depolymerization

eta   = 1d-9;           % (Pa s/nm^2/mM)

Dc  = 1.d8;             % (nm^2/s) diffusion constant for solutes
Dtc = 1.d7;             % (nm^2/s) diffusion constant for theta_c

alphaf = 1.d-1;     % (nm/Pa/s) coefficient of water permeation
alphab = 1.d-1;     % (nm/Pa/s) coefficient of water permeation
gf = 5d4;         % (nm/s) passive channel coefficient at the front
gb = 5d4;         % (nm/s) passive channel coefficient at the back

cinb = 340.3;
cinf = 340.8;
vc0 = 0;           % (nm) expected average cytosol velocity

%%
Iter = 21;
N  = 121;
dx = L/(N-1);
x  = linspace(0,L,N);

% number of variables of the model
n_var = 7; % pc, v_c, v_n, theta_n, theta_c, cin, v0
var_legh = [N,1,N,N,N,N,1]; % length of eavh variable
var_strt = zeros(1,n_var); % start position of each variable
var_strt(1) = 1;
for in = 2:n_var
    var_strt(in) = var_strt(in-1) + var_legh(in-1);
end
N_var = var_strt(n_var) + var_legh(n_var) - 1; % total size of the matrix

s_pc = var_strt(1); l_pc = var_legh(1);
s_vc = var_strt(2); l_vc = var_legh(2);
s_vn = var_strt(3); l_vn = var_legh(3);
s_tn = var_strt(4); l_tn = var_legh(4);
s_tc = var_strt(5); l_tc = var_legh(5);
s_c  = var_strt(6); l_c  = var_legh(6);
s_v0 = var_strt(7); l_v0 = var_legh(7);

%%

nLV = 2;
bLV = 0.3;

nHV = 4;
bHV = 0.25;

etast0_ref = 0.8d-4;           % (Pa s/nm^2/mM) drag coefficient between the actin and the substrate

Ntest = 4;
V0 = zeros(Ntest,1);
for ih = 1:Ntest
    
    if ih == 1
        Mode  = 'LV';
        Model = 3;
    elseif ih == 2
        Mode  = 'LV';
        Model = 1;
    elseif ih == 3
        Mode  = 'HV';
        Model = 3;
    elseif ih == 4
        Mode  = 'HV';
        Model = 1;
    end
    
    if strcmp(Mode,'LV')
        Polar = 'LV';
        ratio = 1.67;              % NHE polarization ratio
        L = 85d3;           % (nm) cell length
        eta_fluid = 0.77d-3;    % (Pa s) fluid dynamic viscosity
        dg    = 12*eta_fluid*(L0-L/1d3)/b^2*1d-3;   % (Pa s/nm) drag coefficient from the displaced water
        gamma0 = 6.0d-4;          % (1/s) constant rate of actin depolymerization
        etast0 = etast0_ref;
        kad = 1.8d-1;           % (Pa s/nm) adhesive force, Fad^b = kad*v0
    elseif strcmp(Mode,'HV')
        Polar = 'HV';
        ratio = 2.84;              % NHE polarization ratio
        L = 125d3;           % (nm) cell length
        eta_fluid = 10d-3;    % (Pa s) fluid dynamic viscosity
        dg    = 12*eta_fluid*(L0-L/1d3)/b^2*1d-3;   % (Pa s/nm) drag coefficient from the displaced water
        gamma0 = 6.0d-4;          % (1/s) constant rate of actin depolymerization
        etast0 = 10*etast0_ref;           % (Pa s/nm^2/mM) drag coefficient between the actin and the substrate
        kad = 9d-1;           % (Pa s/nm) adhesive force, Fad^b = kad*v0
    end
    
    Jactinf0 = 2.5;           % (nm mM/s) Jactinf = Jactinf0*thetac^f/(thetacc + thetac^f)
    
    if Model == 1
        Jcactiveb = 0;  % (mM nm/s) active flux at the front
    elseif Model == 2 || Model == 3
        Jcactiveb = -5.3d6; %1*(gf*(cinf-c0f)+Dcommon/L*(cinf-cinb)-vc0*(cinb+cinf)/2);
    end
    
    dx = L/(N-1);
    x  = linspace(0,L,N);
    
    dgf = dg/2;
    dgb = dg/2;
    
    nx = (0:N-1)'/(N-1);
    if strcmp(Polar,'LV')
        etast = etast0*(2^nLV*(1-bLV)*(abs(nx-1/2)).^nLV + bLV);
    elseif strcmp(Polar,'HV')
        etast = etast0*((1-bHV)*nx.^nHV + bHV);
    end
    
    Jcactivef = -ratio*Jcactiveb;
    
    % initial guess
    pc = ((cinb+cinf)/2-(c0b+c0f)/2)*R*T*ones(N,1);
    thetan = linspace(Thetan*1,Thetan*1,N)';
    thetac = linspace(Thetac*1,Thetac*1,N)';
    Jactinf = Jactinf0*thetac(N)/(thetacc+thetac(N));
    DJactinfDtcN = Jactinf0*thetacc/(thetacc+thetac(N))^2;
    Jwaterf = -alphaf*(pc(N)-p0f-R*T*(cinf-c0f));
    v0 = 0*(etast(1)*Jactinf+dg/L*Jwaterf)/(Thetan*etast(1)+kad/L+dg/L);
    vc = v0 - Jwaterf;
    vn = linspace(v0,v0-Jactinf,N)';
    cin = linspace(cinb,cinf,N)';
    X = [pc; vc; vn; thetan; thetac; cin; v0];
    
    iter = 0;
    ITER = true;
    while ITER
        iter = iter + 1;
        
        DF = zeros(N_var,N_var);
        Fn = zeros(N_var,1);
        
        sigma_n = ksigman*thetan;     % N by 1 vector
        sigma   = sigma_n;
        dsigmadtn =  ksigman;
        
        gamma = gamma0*logspace(-0,0,N)';
        dgammadtn = zeros(N,1);
        dgammadmn = zeros(N,1);
        
        Jactinf = Jactinf0*thetac(N)/(thetacc+thetac(N));
        DJactinfDtcN = Jactinf0*thetacc/(thetacc+thetac(N))^2;
        
        %% Equations for pc and the Derivatives
        Fn(s_pc) = -pc(2)+pc(1) + dx*eta*thetan(1)*(vn(1)-vc);
        DF(s_pc,[s_pc,s_pc+1]) = [1, -1];
        DF(s_pc,s_vc) = -dx*eta*thetan(1);
        DF(s_pc,s_vn) = dx*eta*thetan(1);
        DF(s_pc,s_tn) = dx*eta*(vn(1) - vc);
        
        Fn(s_pc+1:s_pc+N-2) = -pc(3:N)+pc(1:N-2) ...
            + 2*dx*eta*thetan(2:N-1).*(vn(2:N-1)-vc);
        for i = 2:N-1
            DF(s_pc+i-1,[s_pc+i-2, s_pc+i]) = [1, -1];
            DF(s_pc+i-1,s_vc)     = -2*dx*eta*thetan(i);
            DF(s_pc+i-1,s_vn+i-1) =  2*dx*eta*thetan(i);
            DF(s_pc+i-1,s_tn+i-1) = -2*dx*eta*(vc - vn(i));
        end
        
        Fn(s_pc+N-1) = -alphab/alphaf*(pc(1)-p0b) - (pc(N)-p0f)...
            + (dgf-alphab/alphaf*dgb)*(vc + v0) ...
            + alphab/alphaf*R*T*(cin(1)-c0b) + R*T*(cin(N)-c0f);
        DF(s_pc+N-1,[s_pc,s_pc+N-1]) = -[alphab/alphaf, 1];
        DF(s_pc+N-1,s_vc) = dgf - alphab/alphaf*dgb;
        DF(s_pc+N-1,[s_c, s_c+N-1])  = R*T*[alphab/alphaf, 1];
        DF(s_pc+N-1,s_v0) = dgf - alphab/alphaf*dgb;
        
        %% Equations for vc and the Derivatives
        Fn(s_vc) = -(pc(N)-p0f) + (dgf + 1/alphaf)*vc + dgf*v0 + R*T*(cin(N)-c0f);
        
        DF(s_vc,s_pc+N-1) = -1;
        DF(s_vc,s_vc)     = dgf + 1/alphaf;
        DF(s_vc,s_c+N-1)  = R*T;
        DF(s_vc,s_v0) = dgf;
        
        %% Equations for vn and the Derivatives
        Fn(s_vn) = - (sigma(2)-sigma(1)) + dx*eta*thetan(1)*(vc-vn(1)) ...
            - dx*etast(1)*thetan(1)*(vn(1) + v0);
        DF(s_vn,s_vc) = dx*eta*thetan(1);
        DF(s_vn,s_vn) = -dx*eta*thetan(1) - dx*etast(1)*thetan(1);
        DF(s_vn,s_tn) = dsigmadtn + dx*eta*(vc-vn(1)) - dx*etast(1)*(vn(1)+v0);
        DF(s_vn,s_tn+1) = -dsigmadtn;
        DF(s_vn,s_v0) = -dx*etast(1)*thetan(1);
        
        Fn(s_vn+1:s_vn+N-2) = -(sigma(3:N)-sigma(1:N-2)) ...
            + 2*dx*eta*thetan(2:N-1).*(vc-vn(2:N-1)) ...
            - 2*dx*etast(2:N-1).*thetan(2:N-1).*(vn(2:N-1)+v0);
        for i = 2:N-1
            DF(s_vn+i-1,s_vc) = 2*dx*eta*thetan(i);
            DF(s_vn+i-1,s_vn+i-1) = -2*dx*thetan(i)*(eta+etast(i));
            DF(s_vn+i-1,[s_tn+i-2,s_tn+i]) = [1,-1]*dsigmadtn;
            DF(s_vn+i-1,s_tn+i-1) = 2*dx*eta*(vc-vn(i)) - 2*dx*etast(i)*(vn(i)+v0);
            DF(s_vn+i-1,s_v0) = -2*dx*etast(i)*thetan(i);
        end
        
        Fn(s_vn+N-1) = thetan(N)*vn(N) + Jactinf;
        DF(s_vn+N-1,s_vn+N-1) = thetan(N);
        DF(s_vn+N-1,s_tn+N-1) = vn(N);
        DF(s_vn+N-1,s_tc+N-1) = DJactinfDtcN;
        
        %% Equations for thetan and the Derivatives
        Fn(s_tn) = thetan(1)*vn(1);
        DF(s_tn,s_vn) = thetan(1);
        DF(s_tn,s_tn) = vn(1);
        
        Fn(s_tn+1:s_tn+N-2) = (thetan(3:N).*vn(3:N)-thetan(1:N-2).*vn(1:N-2))...
            + 2*dx*gamma(2:N-1).*thetan(2:N-1);
        for i = 2:N-1
            DF(s_tn+i-1,[s_vn+i-2,s_vn+i]) = [-thetan(i-1), thetan(i+1)];
            DF(s_tn+i-1,s_tn+i-2) = -vn(i-1);
            DF(s_tn+i-1,s_tn+i-1) = 2*dx*gamma(i) + 2*dx*dgammadtn(i)*thetan(i);
            DF(s_tn+i-1,s_tn+i)   =  vn(i+1);
        end
        
        Fn(s_tn+N-1) = (thetan(N)*vn(N)-thetan(N-1)*vn(N-1)) ...
            + dx*gamma(N)*thetan(N);
        DF(s_tn+N-1,[s_vn+N-2,s_vn+N-1]) = [-thetan(N-1), thetan(N)];
        DF(s_tn+N-1,s_tn+N-2) = -vn(N-1);
        DF(s_tn+N-1,s_tn+N-1) =  vn(N) + dx*gamma(N) + dx*thetan(N)*dgammadtn(N);
        
        %% Equations for thetac and the Derivatives
        Fn(s_tc) = thetac(1)*vc - Dtc/dx*(thetac(2)-thetac(1));
        DF(s_tc,s_vc)   = thetac(1);
        DF(s_tc,s_tc)   = vc + Dtc/dx;
        DF(s_tc,s_tc+1) = - Dtc/dx;
        
        Fn(s_tc+1:s_tc+N-2) = vc * (thetac(3:N)-thetac(1:N-2)) ...
            - 2*Dtc/dx*(thetac(1:N-2) - 2*thetac(2:N-1) + thetac(3:N))...
            - 2*dx*gamma(2:N-1).*thetan(2:N-1);
        for i = 2:N-1
            DF(s_tc+i-1,s_vc) = (thetac(i+1) - thetac(i-1));
            DF(s_tc+i-1,s_tn+i-1) = -2*dx*gamma(i) - 2*dx*thetan(i)*dgammadtn(i);
            DF(s_tc+i-1,s_tc+i-2) = -vc - 2*Dtc/dx;
            DF(s_tc+i-1,s_tc+i-1) = 4*Dtc/dx;
            DF(s_tc+i-1,s_tc+i)   =  vc - 2*Dtc/dx;
        end
        
        Fn(s_tc+N-1) = 1/2*(sum(thetan(1:N-1)+thetac(1:N-1)) + sum(thetan(2:N)+thetac(2:N)))...
            - L*(Thetan + Thetac)/dx;
        DF(s_tc+N-1,[s_tn,s_tn+N-1]) = 1/2;
        DF(s_tc+N-1,s_tn+1:s_tn+N-2) = 1;
        DF(s_tc+N-1,[s_tc,s_tc+N-1]) = 1/2;
        DF(s_tc+N-1,s_tc+1:s_tc+N-2) = 1;
        
        %% Equations for cin and the Derivatives
        Fn(s_c) = vc*cin(1) - Dc/dx*(cin(2)-cin(1))...
            + gb*(cin(1)-c0b) - Jcactiveb;
        DF(s_c,s_vc)  =  cin(1);
        DF(s_c,s_c)   =  vc + Dc/dx + gb;
        DF(s_c,s_c+1) =     - Dc/dx;
        
        Fn(s_c+1:s_c+N-2) = vc*(cin(3:N)-cin(1:N-2)) ...
            - 2*Dc/dx.*(cin(3:N) - 2*cin(2:N-1) + cin(1:N-2));
        
        for i = 2:N-1
            DF(s_c+i-1,s_vc)    =   cin(i+1) - cin(i-1);
            DF(s_c+i-1,s_c+i-2) = - vc - 2*Dc/dx;
            DF(s_c+i-1,s_c+i-1) =        4*Dc/dx;
            DF(s_c+i-1,s_c+i)   =   vc - 2*Dc/dx;
        end
        
        Fn(s_c+N-1) = vc*cin(N) - Dc/dx*(cin(N)-cin(N-1))...
            - gf*(cin(N) - c0f) + Jcactivef;
        DF(s_c+N-1,s_vc)    =  cin(N);
        DF(s_c+N-1,s_c+N-2) =       Dc/dx;
        DF(s_c+N-1,s_c+N-1) =  vc - Dc/dx - gf;
        
        %% Equation for v0 and the derivatives
        Fn(s_v0) = (fextf-fextb) + (p0f-p0b) + (dgf + dgb)*(vc + v0) + kad*v0 ...
            + dx/2*(sum(etast(1:N-1).*thetan(1:N-1).*vn(1:N-1)) ...
            + sum(etast(2:N).*thetan(2:N).*vn(2:N)))...
            + dx/2*v0*(sum(etast(1:N-1).*thetan(1:N-1)) + sum(etast(2:N).*thetan(2:N)));
        DF(s_v0,s_vc) = dgf + dgb;
        DF(s_v0,[s_vn,s_vn+N-1]) = dx/2*[etast(1)*thetan(1),etast(N)*thetan(N)];
        DF(s_v0,s_vn+1:s_vn+N-2) = dx*etast(2:N-1).*thetan(2:N-1);
        DF(s_v0,[s_tn,s_tn+N-1]) = dx/2*[etast(1)*(vn(1)+v0),etast(N)*(vn(N)+v0)];
        DF(s_v0,s_tn+1:s_tn+N-2) = dx*etast(2:N-1).*(vn(2:N-1)+v0);
        DF(s_v0,s_v0) = (dgf + dgb) + kad ...
            + dx/2*(sum(etast(1:N-1).*thetan(1:N-1)) + sum(etast(2:N).*thetan(2:N)));
        
        temp_Fn = Fn;
        
        %% Solve for the matrix
        DF = sparse(DF);
        X = X - DF\Fn;
        
        if sum(isnan(X)) >= 1 || iter == Iter || norm(temp_Fn) < norm(Fn)
            Solution = 0;
        else
            Solution = 1;
        end
        
        pc = X(s_pc:s_pc+N-1);
        vc = X(s_vc);
        vn = X(s_vn:s_vn+N-1);
        thetan = X(s_tn:s_tn+N-1);
        thetac = X(s_tc:s_tc+N-1);
        cin = X(s_c:s_c+N-1);
        v0 = X(s_v0);
        
        if iter > 1
            error = abs((X-temp_X)./(X+eps));
            error = sum(error)/(N_var);
            if error < 1d-6 || iter == Iter
                ITER = false;
            end
        end
        temp_X = X;
    end
    
    p_star_f = p0f + dgf*(vc + v0);
    Jwaterf = -alphaf*(pc(N)-p_star_f-R*T*(cin(N)-c0f));
    
    sigma_n = ksigman*thetan;     % N by 1 vector
    sigma   = sigma_n;
    
    gamma = gamma0*logspace(-0,0,N)';
    Jactinf = Jactinf0*thetac(N)/(thetacc+thetac(N));
    
    V0(ih) = v0;
    
    if ih == 1 || ih == 3
        
        figure(1)
        plot(x/L,etast/etast0,'linewidth',2); hold on
        set(gca,'fontsize',15);
        xlabel('x (normalized)','fontsize',15)
        ylabel('\eta_{st} (normalized)','fontsize',15)
        axis([0 1 0 1])
        box off
        if ih == 3
            legend('0.77 cp','8 cp')
        end
        
        figure(2)
        plot(x/L,thetan/thetan(1),'linewidth',2); hold on
        set(gca,'fontsize',15);
        xlabel('x (normalized)','fontsize',15)
        ylabel('Normalized \theta_n ({\mu}M)','fontsize',15)
        box off
        if ih == 3
            legend('0.77 cp','8 cp')
        end
        
    end
    
end

figure(3)
X = categorical({'0.77 cp Ctrl','0.77 cp NHE^{-/-}','8 cp Ctrl','8 cp NHE^{-/-}'});
X = reordercats(X,{'0.77 cp Ctrl','0.77 cp NHE^{-/-}','8 cp Ctrl','8 cp NHE^{-/-}'});
bar(X, V0*3.6, 0.6)
set(gca,'fontsize',15);
ylabel('v_0 ({\mu}m/h)','fontsize',15)
ylim([0 100])
box off