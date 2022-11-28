%function Path_loss(bb)
R=0.2; %the cell radius(KM)
R0=0.02; %the closest distance the mobile can be from the the BS antenna
hB=0.01; % BS antenna height
hm=0.002; % mobile antenna height


fc=9*10^8; %carrier frequence
Ru_iteration=0; 
a=2;
b=2;
K=1;
St=1;
NI=6;
Ns =10;
B = 0.01;
W =1;
pa  =B^(1/Ns);
Number_of_i = 1:1:NI;
P_nI = binopdf(Number_of_i,6,pa);

lamdac=(3*10^5)/fc; %calculate lamdac
g=(4*hB*hm)/lamdac; %calculate the break point of the pass-loss curve

for Ru=2:0.2:10;  % this is the normalized Ru
D=Ru.*R; %reuse distance
Ru_iteration=Ru_iteration+1;

for i=1:1:10000; % repeat the process for 10,000 times
%%%%%%%%%% STEP 1%%%%%
u=rand(1,1);
r= R0+(R - R0)*sqrt(u); % Step 1, the Random Distance of the user
%%%%%%%%%%STEP 2%%%%%%%%%%%

for j = 1:1:NI;

ui = rand(j,1);  %Step 2, ,6 means NI = 6

vi = rand(j,1);  % Step 2,
xi = R0+(R - R0)*sqrt(ui);  %the distance of the of the interferer UI
thetai=2*pi*vi; %%the distance of the of the interferer UI

%%%%%%%%%%STEP 3%%%%%%%%%%%
ri=(D*D+(xi).^2+(2*D*xi).*sin(thetai)).^(1/2); %the actual distance of the interferer

%%%%%%%%%%STEP 4 & Step 5%%%%%%%%%%%
Sd=(K/((r^a)*(1+r/g)^b))*St;
%S_R_interfered  = K*S_transmitted./((ri.^a).*(1+ri/g).^b);
S_R_interfered =(K./((ri.^a).*(1+ri./g).^b))*St; 
SI=sum(S_R_interfered);
gammad(j)  = Sd/(sum(S_R_interfered));



end


for j = 1:1:NI;
    uii_b=1./((D+R)^a*(1+((D+R)/g))^b)*St;
    uii_w=1./((D-R)^a*(1+((D-R)/g))^b)*St;
    gammad_partial_best(j) = Sd/NI/uii_b;
    gammad_partial_worst(j) = Sd/NI/uii_w;
end


Ck = 1/Ns*log2(1+gammad).*P_nI;
Ckb = 1/Ns*log2(1+gammad_partial_best).*P_nI;
Ckw = 1/Ns*log2(1+gammad_partial_worst).*P_nI;



Ae_partial(Ru_iteration,i) = 4*Ns*B^(1/Ns)/(pi*W*(Ru.^2)*(R*R))*sum(Ck);
Ae_partial_b(Ru_iteration,i) = 4*Ns*B^(1/Ns)/(pi*W*(Ru.^2)*(R*R))*sum(Ckb);

Ae_partial_w(Ru_iteration,i) = 4*Ns*B^(1/Ns)/(pi*W*(Ru.^2)*(R*R))*sum(Ckw);

gammad  = Sd/(sum(S_R_interfered));
gammab=((((R*(Ru+1))/r)^a)*((g+(Ru+1)*R)/(g+r))^b)/NI;%the best case of gammad
gammaw=((((R*(Ru-1))/r)^a)*((g+(Ru-1)*R)/(g+r))^b)/NI;%the worst case of gammad
%%%%%%%%%% STETP 6 %%%%%%%%%%%
Ae(Ru_iteration,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammad)); %the ASE of genera simulations
Aeb(Ru_iteration,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammab)); % ASE of the best cas interference
Aew(Ru_iteration,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammaw)); % ASE ofthe worst cas interference

end
end
%%%%%%%% PLOT 
aAe=mean(Ae,2); %the average of Ae after running 10000 times
aAeb=mean(Aeb,2); %the average of the best case after running 10000 times
aAew=mean(Aew,2); %the average of the worst case after running 10000 times

aAe_partial=mean(Ae_partial,2); %the average of the worst case after running 10000 times
aAe_partial_b =mean(Ae_partial_b,2); %the average of the worst case after running 10000 times
aAe_partial_w =mean(Ae_partial_w,2); %the average of the worst case after running 10000 times

Ru=2:0.2:10;
figure()
plot(Ru,aAe,'k',Ru,aAeb,'g',Ru,aAew,'r',Ru,aAe_partial,'k--',Ru,aAe_partial_b,'g--',Ru,aAe_partial_w,'r--');%set the best case is green,the worst line is red
legend('B=1, Fully-loaded System ','B=1, Fully-loaded System','B=1, Fully-loaded System','B=0.01, Ns=10','B=0.01, Ns=10','B=0.01, Ns=10')
xlabel('Normalized Reuse Distance Ru');grid;
ylabel('ASE[Bits/Sec/Hz/Km^2]');
title('ASE in Uplink Transmission(a=2,b=2,R=200m),');
%end