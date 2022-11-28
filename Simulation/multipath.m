%function Path_loss(bb)
R=200; %the cell radius(KM)
R0=20; %the closest distance the mobile can be from the the BS antenna
hB=10; % BS antenna height
hm=2; % mobile antenna height


Ru_iteration=0; 
a=2;
b=2;
K=1;
St=1;
NI=6;
Fading_md = 3;
Fading_mi = 1;

fc=9*10^8; %carrier frequence
lamdac=(3*10^8)/fc; %calculate lamdac
g=(4*hB*hm)/lamdac; %calculate the break point of the pass-loss curve

rho = Fading_md +Fading_mi * NI;

for Ru=2:0.2:10;  % this is the normalized Ru
    D=Ru.*R; %reuse distance
    %
    Ru_iteration=Ru_iteration+1;
    
    for i=1:1:10000; % repeat the process for 10,000 times
        
        %%%%%%%%%% STEP 1%%%%%
        u=rand(1,1);
       % u_mean=rand(1000,1);
        r= R0+(R - R0)*sqrt(u); % Step 1, the Random Distance of the user
       % rmean = R0+(R - R0).*sqrt(u_mean);
        %%%%%%%%%%STEP 2%%%%%%%%%%%
        ui = rand(6,1);  %Step 2, ,6 means NI = 6
        vi = rand(6,1);  % Step 2,
        xi = R0+(R - R0)*sqrt(ui);  %the distance of the of the interferer UI
        thetai=2*pi*vi; %%the distance of the of the interferer UI
        
        %%%%%%%%%%STEP 3%%%%%%%%%%%
        ri=(D*D+(xi).^2+(2*D*xi).*sin(thetai)).^(1/2); %the actual distance of the interferer


        %%%%%%%%%%STEP 4 & Step 5%%%%%%%%%%%
        Sd=(K/((r^a)*(1+r/g)^b))*St;
      %  Sd_mean = mean((K./((rmean.^a).*(1+rmean./g).^b))*St);
        S_R_interfered = (K./((ri.^a).*(1+ri./g).^b))*St; 
        SI = sum(S_R_interfered);

        S_R_interfered_worst = K./((D-R)^a*(1+((D-R)/g))^b)*St;
        S_R_interfered_best = K./((D+R)^a*(1+((D+R)/g))^b)*St;

        gamma_a_d = Fading_md;
        gamma_a_i = Fading_mi;

        gamma_b_d = Sd/Fading_md;
        gamma_b_i = S_R_interfered/Fading_mi;

        PDF_SD = gamrnd(gamma_a_d,gamma_b_d,1,1);
        PDF_SI = gamrnd(gamma_a_i,gamma_b_i,6,1);

        
        gamma_b_i_best = S_R_interfered_best/Fading_mi;
        gamma_b_i_worst = S_R_interfered_worst/Fading_mi;




        PDF_SI_best = gamrnd(gamma_a_i,gamma_b_i_best,6,1);
        PDF_SI_worst = gamrnd(gamma_a_i,gamma_b_i_worst,6,1);


       % PDF_SD = (Fading_md/Sd_mean)^Fading_md*Sd^(Fading_md - 1)/gamma(Fading_md)*exp(-Fading_md*Sd/Sd_mean);
        %PDF_SI = (Fading_mi/SI/NI)^Fading_mi.*S_R_interfered.^(Fading_mi - 1)/gamma(Fading_mi).*exp(-Fading_mi*S_R_interfered/SI/NI);

     
      %  gammad_multi  = Sd*PDF_SD/(sum(S_R_interfered.*PDF_SI));
        gammad_multi = PDF_SD/sum(PDF_SI);
        gammad_multi_worst  = PDF_SD/sum(PDF_SI_worst);
        gammad_multi_best  = PDF_SD/sum(PDF_SI_best);
        
        gammad  = Sd/(sum(S_R_interfered));
        gammab=((((R*(Ru+1))/r)^a)*((g+(Ru+1)*R)/(g+r))^b)/NI;%the best case of gammad
        gammaw=((((R*(Ru-1))/r)^a)*((g+(Ru-1)*R)/(g+r))^b)/NI;%the worst case of gammad
        %%%%%%%%%% STETP 6 %%%%%%%%%%%
        Ae(Ru_iteration,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1 + gammad)); %the ASE of genera simulations
        Aeb(Ru_iteration,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammab)); % ASE of the best cas interference
        Aew(Ru_iteration,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammaw)); % ASE ofthe worst cas interference
         Ae_multi(Ru_iteration,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1 + gammad_multi)); %the ASE of genera simulations
        Ae_multi_worst(Ru_iteration,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1 +gammad_multi_worst)); % ASE of the best cas interference
        Ae_multi_best(Ru_iteration,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1 +gammad_multi_best)); % ASE ofthe worst cas interference
        
    end
end
%%%%%%%% PLOT 
aAe=mean(Ae,2)*10^6; %the average of Ae after running 10000 times%
aAeb=mean(Aeb,2)*10^6; %the average of the best case after running 10000 times
aAew=mean(Aew,2)*10^6; %the average of the worst case after running 10000 times

aAe_mul=mean(Ae_multi,2)*10^6; %the average of Ae after running 10000 times%
aAe_mulb=mean(Ae_multi_best,2)*10^6; %the average of the best case after running 10000 times
aAe_mulw=mean(Ae_multi_worst,2)*10^6; %the average of the worst case after running 10000 times

Ru=2:0.2:10;
figure()
plot(Ru,aAe,'k',Ru,aAeb,'g',Ru,aAew,'r',Ru,aAe_mul,'r--',Ru,aAe_mulb,'p--',Ru,aAe_mulw,'b--');%set the best case is green,the worst line is red

%plot(Ru,aAe,'k',Ru,aAeb,'g',Ru,aAew,'r');%set the best case is green,the worst line is red

%plot(Ru,aAe,'k',Ru,aAeb,'g',Ru,aAew,'r',Ru,aAe_mul,'r--');

legend('Simulations ','Best-case Interference','Worst-case Interference','Multipath ','Best-case Multipath','Worst-case Multipath')
%legend('Simulations ','Best-case Interference','Worst-case Interference','Multipath')%,'Multipath ','Best-case Multipath','Worst-case Multipath')

%legend('Simulations ')
xlabel('Normalized Reuse Distance Ru');grid;
ylabel('ASE[Bits/Sec/Hz/Km^2]');
title('ASE in Uplink Transmission(a=2,b=2,R=200m),');
%end