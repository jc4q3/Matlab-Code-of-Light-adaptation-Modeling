close all
clc
Y=csvread('SP4.csv');
v=1:10:1000;
dtt1 = 1e-3;
dtt2 = 1e-2;
tt=[dtt1:dtt1:1,1+dtt2:dtt2:300,300+dtt1:dtt1:301,301+dtt2:dtt2:331,331+dtt1:dtt1:332,332+dtt2:dtt2:362,362+dtt1:dtt1:363,363+dtt2:dtt2:393,393+dtt1:dtt1:394];
Ydata=[Y(v);Y(1001:30900);Y(30900+v);Y(31901:34900);Y(34900+v);Y(35901:38900);Y(38900+v);Y(39901:42900);Y(42900+v)]';
samplinginterval=0.01;
t=samplinginterval:samplinginterval:394;
k=[0.00167329
13.45464086
151.8544403
0.588945609
580.9436628
4196.281928
971.3031091
571.7923856
619.4650614
174.1330374
5.241173076
165.2368658
51.78179117
4.542231291
19.65396685
63.09051233
0.854439311
0.004136785
0.002584267
24.12541449
3.363765926
47.01892247
64.01218561
1299.629788
1.146396523
1.616404233
1.49104289
1.218063986
26.57900245
41.05731863
11.28942854
250.6944121
];

[Fluo,x,beta,R_psbs,R_VDE,delta_psi,delta_pH_part,pmf]=PulseVariedflueqnRK4_PSBS(k,t);
Yold=Fluo;
figure(1)
plot(t,Ydata);
% ylim([40 500])
hold on
plot(t,Yold);
legend('Fluorescence experimental data','Fluorescence estimated');
xlabel('time');
ylabel('Fluorescence');

figure(2)
semilogx(t,Ydata);
% ylim([40 500])
hold on
semilogx(t,Yold);
legend('Fluorescence experimental data','Fluorescence estimated');
xlabel('time');
ylabel('Fluorescence');


figure(3)
subplot(3,3,1)
plot(t,x(1,:))
xlabel('time');
ylabel('A^*');
subplot(3,3,2)
plot(t,x(2,:))
xlabel('time');
ylabel('QA^-');
subplot(3,3,3)
plot(t,x(3,:))
xlabel('time');
ylabel('QB^-');
subplot(3,3,4)
plot(t,x(4,:))
xlabel('time');
ylabel('QB^2^-');
subplot(3,3,5)
plot(t,x(5,:))
xlabel('time');
ylabel('PQ');
subplot(3,3,6)
plot(t,x(6,:))
xlabel('time');
ylabel('pH');
subplot(3,3,7)
plot(t,beta)
xlabel('time');
ylabel('Buffer Capacity');
subplot(3,3,8)
semilogx(t,R_psbs)
xlabel('time');
ylabel('R_p_s_b_s');
subplot(3,3,9)
plot(t,R_VDE)
xlabel('time');
ylabel('R_V_D_E');


figure(4)
subplot(3,3,1)
plot(t,x(8,:))
xlabel('time');
ylabel('K^+');
subplot(3,3,2)
plot(t,x(9,:))
xlabel('time');
ylabel('Cl^-');
subplot(3,3,3)
semilogx(t,x(10,:))
xlabel('time');
ylabel('Anth');
subplot(3,3,4)
semilogx(t,x(11,:))
xlabel('time');
ylabel('Zea');
subplot(3,3,5)
semilogx(t,x(12,:))
xlabel('time');
ylabel('Lut');
subplot(3,3,6)
plot(t,x(13,:))
xlabel('time');
ylabel('A_C');
subplot(3,3,7)
plot(t,delta_psi)
xlabel('time');
ylabel('\Delta\psi');
subplot(3,3,8)
plot(t,delta_pH_part)
xlabel('time');
ylabel('\DeltapH Component');
subplot(3,3,9)
plot(t,pmf)
xlabel('time');
ylabel('pmf');


n=length(Ydata);
MAPE=1/n*sum(abs((Ydata-Yold)./Yold))*100
r = Ydata-Yold;
normr = norm(r);
SSE = normr.^2;
SST = norm(Ydata-mean(Ydata))^2;
R2 = 1 - SSE/SST