function [slope,A]=Variedflueqn(t,x,k,u)
A0=290;
PQ=9.04;
% if x(2)<0
%     x(2)=0;
% end
if x(8)<0
    x(8)=0;
end
if x(11)<0
    x(11)=0;
end
%k(3)=0;
km=k(11);
V_max_VA=k(21);
km_VA=k(22);
V_max_AZ=k(23);
km_AZ=k(24);
V_max_ZA=k(27);
km_ZA=k(28);
V_max_AV=k(25);
km_AV=k(26);

pKa_buffer=5.5;

% km1=k(22);
% N_reactionCenter=1.2*10^11;
N_reactionCenter=2.5*10^6;
%convert m^3 of Volume to L; 
% V_lumen=2.1*10^(-18)*1000;
V_lumen=2.1*10^(-15);
constant=6.022*10^23;
e=1.602176634 *10^(-19);
% area in unit cm2
area_thylakoid=(5.4*10^(-6))/2;
capacitance=0.6*10^(-6)*area_thylakoid;
R=8.3145;
Temp=293;
F=96485.33;
kp=10^(-6.5);
np_psbs=1.36;
kp_VDE=10^(-6.5);
np_VDE=4;
% u_max=10^4;
pH0=7.8;
%in stroma, 150mM K+, 50-90mM Cl-, 2-3 mM Mg2+ in light, 0.5-1 mM Mg2+ in
%dark, 40�70 mM Na+ depending on the plant species and salt condition. 
% in lumen, 40mM Cl-, 135 mM K+
%
c_K=0.15;
c_Cl=0.075;
% R_psbs=10^(-x(6)*np_psbs)/(kp^np_psbs+10^(-x(6)*np_psbs));
% A_initial=A0*k_a*0.3*(1-R_psbs)*(u_max-u)/u_max;

% R_psbs=10.^(-x(6)*np_psbs)./(kp^np_psbs+10.^(-x(6)*np_psbs));

%*R_psbs.

% A=k(22)*(1-R_psbs)*A0;

% A=Amin+(Amax-Amin)/(k(23)*R_psbs+1);
% A=A0;
A=A0-x(13);
% r1=k(1)*u(1)*(A-x(1))-k(2)*x(1)-k(4)*x(1)*(1-x(2))+k(5)*x(2)*(A-x(1));
R_psbs=10^(-x(6)*np_psbs)/(kp^np_psbs+10^(-x(6)*np_psbs));

%both pathway 1 and 2
r1=k(1)*u(1)*(A-x(1))-k(2)*x(1)-R_psbs*(k(15)*x(12)+k(16)*x(11)+k(17)*x(10))*x(1)-k(3)*x(1)*(1-x(2))+k(4)*x(2)*(A-x(1));

%without pathway 1
% r1=k(1)*u(1)*(A-x(1))-k(2)*x(1)-k(4)*x(1)*(1-x(2))+k(5)*x(2)*(A-x(1));

r2=k(3)*x(1)*(1-x(2))-k(4)*x(2)*(A-x(1))-k(5)*x(2)*(1-x(3)-x(4))+k(6)*x(3)*(1-x(2))-k(7)*x(2)*x(3)+k(8)*x(4)*(1-x(2));
r3=k(5)*x(2)*(1-x(3)-x(4))-k(6)*x(3)*(1-x(2))-k(7)*x(2)*x(3)+k(8)*x(4)*(1-x(2));
r4=k(7)*x(2)*x(3)-k(8)*x(4)*(1-x(2))-k(9)*x(4)*x(5);  
% beta=0.2;
beta=2.303*10^(-x(6))+0.69078*10^(-pKa_buffer)*10^(-x(6))/(10^(-pKa_buffer)+10^(-x(6)))^2;

% accumulated_H=10^(-x(6))+0.3*10^(-x(6))/(10^(-5.5)+10^(-x(6)))-10^(-pH0)-0.3*10^(-pH0)/(10^(-5.5)+10^(-pH0));
r5=-k(9)*x(4)*x(5)+k(10)*x(6)/7.8*(PQ-x(5))/(km+(PQ-x(5)));
%K+concentration in stroma is about 135mmol/L

delta_psi=e*(x(7)+x(8)-x(9)-c_K+c_Cl)*constant*V_lumen/capacitance;
% delta_psi=0;
r8=k(12)*(-delta_psi+R*Temp/F*log(c_K/x(8)));
r9=-k(13)*(-delta_psi-R*Temp/F*log(c_Cl/x(9)));
delta_pH_part=2.3*R*Temp/F*(8-x(6));
pmf=delta_psi+delta_pH_part;


r6=-1/beta*(N_reactionCenter/V_lumen/constant*(k(3)*x(1)*(1-x(2))+4*k(10)*x(6)/7.8*(PQ-x(5))/(km+(PQ-x(5))))-k(14)*pmf);
r7=N_reactionCenter/V_lumen/constant*(k(3)*x(1)*(1-x(2))+4*k(10)*x(6)/7.8*(PQ-x(5))/(km+(PQ-x(5))))-k(14)*pmf;
% proton_production=10^(-x(6))+0.3*10^(-x(6))/(10^(-5.5)+10^(-x(6)))-10^(-pH0)-0.3*10^(-pH0)/(10^(-5.5)+10^(-pH0));
% r6=-1/beta*(N_reactionCenter/V_lumen/constant*(k(4)*x(1)*(1-x(2))+4*k(12)*(PQ-x(5))/(km+(PQ-x(5))))-k(14)*pmf);

R_VDE=10^(-x(6)*np_VDE)/(kp_VDE^np_VDE+10^(-x(6)*np_VDE));
% if pmf<0.
%     R_ZE=k(33);
% else
    R_ZE=0;
% end



r10=V_max_VA*R_VDE*(1-x(10)-x(11))/(km_VA+1-x(10)-x(11))-V_max_AZ*R_VDE*x(10)/(km_AZ+x(10))+V_max_ZA*R_ZE*(x(11))/(km_ZA+x(11))-V_max_AV*R_ZE*(x(10))/(km_AV+x(10));
%without A complex
% r11=V_max_AZ*R_VDE*x(10)/(km_AZ+x(10))-V_max_ZA*R_ZE*(x(11))/(km_ZA+x(11));
% r12=k(23)*R_VDE*(1-x(12))/(k(24)+1-x(12));
% r13=0;

%with A complex
r11=V_max_AZ*R_VDE*x(10)/(km_AZ+x(10))-V_max_ZA*R_ZE*(x(11))/(km_ZA+x(11))-k(19)*x(11)*(A0-x(13))/A0;
r12=k(29)*R_VDE*(1-x(12))/(k(30)+1-x(12))-k(20)*x(12)*(A0-x(13))/A0;

r13=k(18)*R_psbs*(k(20)*x(12)+k(19)*x(11))*(A0-x(13));

slope=[r1;r2;r3;r4;r5;r6;r7;r8;r9;r10;r11;r12;r13];
