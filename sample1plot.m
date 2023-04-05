close all
Y=csvread('NSP2_4.csv');
v=1:10:1000;
dtt1 = 1e-3;
dtt2 = 1e-2;
tt=[dtt1:dtt1:1,1+dtt2:dtt2:300,300+dtt1:dtt1:301,301+dtt2:dtt2:331,331+dtt1:dtt1:332,332+dtt2:dtt2:362,362+dtt1:dtt1:363,363+dtt2:dtt2:393,393+dtt1:dtt1:394];
Ydata=[Y(v);Y(1001:30900);Y(30900+v);Y(31901:34900);Y(34900+v);Y(35901:38900);Y(38900+v);Y(39901:42900);Y(42900+v)]';
samplinginterval=0.01;
t=samplinginterval:samplinginterval:394;
k=[0.00167329
24.6197192
508.7713415
1.106640599
851.7895347
4165.682974
1151.025832
624.3902019
49.07614128
124.1815196
5.241172031
162.2160191
46.63350533
6.279215376
34.45119816
74.50352
0.111671282
0.005613246
0.028946165
17.41037772
0.036100138
0.093053513
70.07071569
1299.50877
1.146396523
1.616404233
1.49104289
1.218063986
89.97367605
137.0296629
17.04171653
480.6944121
];
[Fluo,x,beta,R_psbs,R_VDE,delta_psi,delta_pH_part,pmf]=PulseVariedflueqnRK4_PSBS(k,t);
Yold=Fluo;
figure(1)
plot(t,Ydata);

hold on
plot(t,Yold);
legend('Fluorescence experimental data','Fluorescence estimated');
xlabel('time');
ylabel('Fluorescence');

figure(2)
semilogx(t,Ydata);

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
