function [Fluo,x,beta,R_psbs,R_VDE,delta_psi,delta_pH_part,pmf,N]=PulseVariedflueqnRK4_PSBS(k,t) 
h=0.0001;
samplinginterval=t(2)-t(1);
t_end=394;
t_integration=h:h:t_end;
n=length(t_integration);
x=zeros(14,n); 

% %dark adaptation initial
x(6,1)=9.04;
x(7,1)=7.8;
x(8,1)= -0.094;
x(9,1)=0.22;
x(10,1)=0.051;


%Light adaptation initial
% x(:,1)=[0.394445165069418;0.0140329656670052;0.0534078475477987;5.14442227068724e-05;9.03999953449232;8.88854571557125;0;0.135068875350609;0;0];

%PRBS initial
% x(:,1)=[0.204330320307304;0.0137872778534449;0.225882319825202;0.00126543019413446;8.97361778624062;4.07419778802545;0.982348796873557];
%gain=1;

% % light intensity is 2500*10^(-6)mol/m^2/s;
% % Projected area of chloroplast is 15.5 +- 1.2 um2;
area=15.5*10^(-12);
% % one chloroplast absorbs 0.3 of incident radiation regardless of the wavelength;
% % photons in one chloroplast are 6.9982*10^9 photons/s;
Photon=6.9982*10^9;
% % 10^9 chlorophyll in one chloroplast;
% % 400 chlorophyll in one reaction center;
% % 2.5*10^6 reaction center in one chloroplast;
% N_reactionCenter=1.2*10^11;
N_reactionCenter=2.5*10^6;
% N_photonRection=Photon/N_reactionCenter;
N_photonRection=1.12*10^3;
% N_reactionCenter=2.5*10^6;
% % photons in one reaction center are 2.7993*10^3 photons/s;
% % NADP+ concentration is 0.6 to 0.8 mM.
% u=[2.7993*10^3;0;1;1];
%convert m^3 of Volume to L; 
% V_lumen=2.1*10^(-18)*1000;
V_lumen=2.1*10^(-15);
constant=6.022*10^23;
e=1.602176634 *10^(-19);
area_thylakoid=(5.4*10^(-6))/2;
capacitance=0.6*10^(-6)*area_thylakoid;
R=8.3145;
Temp=293;
F=96485.33;
pH0=7.8;


%constant illumination
n_constant=length(h:h:300);
addtime=0;
% Pulse high level
U_Pulse=2.23944*10^3;

U_total=zeros(1,2*n);
U_total(1,1:2*n_constant)=N_photonRection;

Pulse=ones(1,2*length(h:h:1))*U_Pulse;
Constant=ones(1,2*length(h:h:30))*N_photonRection;
PulseConstant=[Pulse Constant];
U_total(1,2*n_constant+1:end-2*addtime*1/h)=[PulseConstant PulseConstant PulseConstant Pulse];
% U_total(1,2*n_constant+1:end-2*addtime*1/h)=10;
U_total(1,end-2*addtime*1/h+1:end)=N_photonRection;

kp_VDE=10^(-6.5);
np_VDE=4;

kp=10^(-6.5);
np_psbs=1.36;

pKa_buffer=5.5;
%autoregression filter

% np_VDE=4;
for i=1:n-1
    % PRBS
%     if U_total(:,2*i)==10
%         x(2,i)=0;
%     end
    [S(:,1),N1] = Variedflueqn(t_integration(i),x(:,i),k,U_total(:,2*i));
    [S(:,2),N2] = Variedflueqn(t_integration(i)+h/2,x(:,i)+S(:,1)*(1/2)*h,k,U_total(:,2*i+1));
    [S(:,3),N3] = Variedflueqn(t_integration(i)+h/2,x(:,i)+S(:,2)*(1/2)*h,k,U_total(:,2*i+1));
    [S(:,4),N4] = Variedflueqn(t_integration(i)+h,x(:,i)+S(:,3)*h,k,U_total(:,2*i+2));
    x(:,i+1) = x(:,i) + h*(S(:,1)+2*S(:,2)+2*S(:,3)+S(:,4))/6;
    N(:,i+1)=(N1+N2+N3+N4)/4;
%     if x(1,i+1)<0
%         x(1,i+1)=0;
%     end
%     if x(2,i+1)<0
%         x(2,i+1)=0;
%     end
%     if x(3,i+1)<0
%         x(3,i+1)=0;
%     end
%     if x(4,i+1)<0
%         x(4,i+1)=0;
%     end
    %dbstop if naninf
end
constant=6.022*10^23;
% Fluo=k(16)*(k(2)*X(1,:)*N_reactionCenter)/(area*constant/10^6);
w=round(1:samplinginterval/h:n);

% beta=zeros(1,w);
% R_psbs=zeros(1,w);
% R_VDE=zeros(1,w);
% quencher=zeros(1,w);
% delta_psi=zeros(1,w);
% delta_pH_part=zeros(1,w);
% pmf=zeros(1,w);

x=x(:,w);
N=N(:,w);
unit3=V_lumen*constant;
c_K=0.15;
c_Cl=0.075;
delta_psi=e*(x(8,:)+x(9,:)-x(10,:)-c_K+c_Cl)*constant*V_lumen/capacitance;
%     delta_psi(:,i+1)=0;
delta_pH_part=2.3*R*Temp/F*(7.8-x(7,:));
pmf=delta_psi+delta_pH_part;
%     beta(i+1)=0.03;

% beta=zeros(1,length(w));
% beta(find(beta==0))=0.3;

beta=2.303*10.^(-x(7,:))+(0.69078*10^(-pKa_buffer)*10.^(-x(7,:)))./(10^(-pKa_buffer)+10.^(-x(7,:))).^2;
    
R_psbs=10.^(-x(7,:)*np_psbs)./(kp^np_psbs+10.^(-x(7,:)*np_psbs));
R_VDE=10.^(-x(7,:)*np_VDE)./(kp_VDE^np_VDE+10.^(-x(7,:)*np_VDE));
%     if t_integration(i)>10
%         x(7,i+1)=k(15)*R_psbs(i+1)*(x(9,i+1)+x(10,i+1))+1;
%     else

% accumulated_H=10.^(-x(6,:))+0.3*10.^(-x(6,:))./(10^(-5.5)+10.^(-x(6,:)))-10^(-pH0)-0.3*10^(-pH0)/(10^(-5.5)+10^(-pH0));

Fluo=k(31)*(k(2)*x(2,:)*N_reactionCenter)/(area*constant/10^6)+k(32);