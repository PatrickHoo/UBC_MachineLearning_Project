R=200; %the cell radius(KM)
R0=20; %the closest distance the mobile can be from the the BS antenna
hB=10; % BS antenna height
hm=2; % mobile antenna height
fc=9*10^8; %carrier frequence
c=0;
a=2;
b=2;
K=1;
St=1;
NI=6;
Xi=10/log(10); %the constant number of xi
sigmad = 4;
sigmaI = 4;
W0=1;
for Ru=2:0.2:10;
    D=Ru.*R;
    c=c+1;
    for i=1:1:10000;
        u=rand(1,1);
        r=R0+(R-R0)*(u^(1/2)); %user's distance to the BS
        ui=rand(6,1); %interferer's random distance to the BSi
        vi=rand(6,1); %interferer's random angle to the BSi
        xi=R0+(R-R0)*(ui.^(1/2)); %interferer's position to the BSi
        thetai=2*pi*vi; %interferer'angle to the BSi
        ri=(D*D+(xi).^2+(2*D*xi).*sin(thetai)).^(1/2); %the distance from the ith interferer to theBS
        lamdac=(3*10^8)/fc; %calculate lamdac
        g=(4*hB*hm)/lamdac; %calculate the break point of the pass-loss curve
        Sd=(K/((r^a)*(1+r/g)^b))*St; % equation 1!(St)(K)
        Si=(K./((ri.^a).*(1+ri./g).^b))*St ; %6 of received power level from ith interfering mobile
        SI=sum(Si); %total interfering power
        Sib=(K/(((D+R)^a)*(1+(D+R)/g)^b))*St ; %6 of received power level from ith interfering mobile
        SIb=6*(Sib); %total interfering power
        Siw=(K/(((D-R)^a)*(1+(D-R)/g)^b))*St ; %6 of received power level from ith interfering mobile
        SIW=6*(Siw); %total interfering power
        udd=1/(r.^a.*(1+(r/g)).^b)*St;
        uii=1./(ri.^a.*(1+(ri/g)).^b)*St;
        
        udd=1/(r.^a.*(1+(r/g)).^b)*St;
        uii=1./(ri.^a.*(1+(ri/g)).^b)*St;
        mdd=log(udd);
        mii=log(uii);
        log_ud=lognrnd(mdd,sigmad/Xi,[1,1]);
        log_ui=lognrnd(mii,sigmaI/Xi,[6,1]);
        gammadi=(log_ud)/sum(log_ui);
        %Calculate log parameters best-case
        udd_b= (K/((r^a)*(1+r/g)^b))*St;
        uii_b=(K/(((D+R)^a)*(1+(D+R)/g)^b))*St;
        mdd_b=log(udd_b);
        mii_b=log(uii_b);
        %log function for desired user simulation   
        log_udb=lognrnd(mdd_b,sigmad/Xi,[1,1]);
        log_uib=lognrnd(mii_b,sigmaI/Xi,[6,1]);
        gammadib =log_udb/(sum(log_uib));



        %Calculate log parameters best-case
        udd_w=(K/((r^a)*(1+r/g)^b))*St;
        uii_w=(K/(((D-R)^a)*(1+(D-R)/g)^b))*St;
        mdd_w=log(udd_w);
        mii_w=log(uii_w);
        %log function for desired user simulation
        log_udw=lognrnd(mdd_w,sigmad/Xi,[1,1]);
        log_uiw=lognrnd(mii_w,sigmaI/Xi,[6,1]);
        gammadiw=log_udw/(sum(log_uiw)); 
        
        
        
        %Calculate sigmaSI
        sigmaSI=(Xi^2)*(log((NI-1+exp(sigmaI^2/(Xi^2)))/NI));%log mean variance of interference
        sigmard=sigmad^2+sigmaSI; %log variance of desired user
        Urdu=Xi.*log(((R*(Ru+1)/r)^a).*((g+(Ru+1)*R)/(g+r)).^b)-Xi*log(NI)+(sigmaSI-sigmaI^2)/(2*Xi);%the upper log mean
        Urdl=Xi.*log(((R*(Ru-1)/r).^a).*((g+(Ru-1)*R)/(g+r)).^b)-Xi*log(NI)+(sigmaSI-sigmaI^2)/(2*Xi);%the lower log mean
        
        %Calculate Cup+andCup-
        Cupu=W0*log2(exp(1))*(Urdu/Xi+exp(sigmard/(2*Xi^2)-Urdu/Xi)); %the upper bound +capacity on desired user
        Cupl=W0*log2(exp(1))*(Urdl/Xi+exp(sigmard/(2*Xi^2)-Urdl/Xi));%the upper bound -capacity on desired user
       

    
        %Calculate Clow+ and Clow-
        Au=Urdu/sqrt(sigmard); %calculate Au in order to use qfunc
        Al=Urdl/sqrt(sigmard);
        Bu=Urdu/sqrt(sigmard)+sqrt(sigmard)/Xi;%calculate Buin order to use qfunc
        Bl=Urdl/sqrt(sigmard)+sqrt(sigmard)/Xi;
        %Calculate qfunc
        C1u=qfunc(Au);
        C2u=qfunc(Bu);
        C1l=qfunc(Al);
        C2l=qfunc(Bl);

        Clowu=W0*log2(exp(1))*(Urdu/Xi+C1u-exp(Urdu/Xi+sigmard/(2*Xi^2))*C2u);%the lower bound+
        Clowl=W0*log2(exp(1))*(Urdl/Xi+C1l-exp(Urdl/Xi+sigmard/(2*Xi^2))*C2l);%the lower bound-

        ACupu(c,i)=4*Cupu./(pi*W0*Ru.^2*R^2);
        ACupl(c,i)=4*Cupl./(pi*W0*Ru.^2*R^2);
        AClowu(c,i)=4*Clowu/(pi*W0*Ru^2*R^2);
        AClowl(c,i)=4*Clowl/(pi*W0*Ru^2*R^2);

        Aeld(c,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammadi)); %the general simulations
        Aebldb(c,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammadib)); %the best case interference
        Aewldw(c,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammadiw)); %the worst case interference
        end
end
aAe=mean(Aeld,2)*10^6;%the average of Ae after running 10000 times
aAeb=mean(Aebldb,2)*10^6;%the average of the best of Ae after runnin 10000 times
aAew=mean(Aewldw,2)*10^6;%the average of the worst of Ae after running 10000 times
aCupu=mean(ACupu,2)*10^6;
aCupl=mean(ACupl,2)*10^6;
aClowu=mean(AClowu,2)*10^6;
aClowl=mean(AClowl,2)*10^6;
Ru=2:0.2:10;


Rb=3:0.2:10;

%aaClowu = aClowu(6:41);

aaClowl = aClowl(6:41);
aaCupl = aCupl(6:41);

%plot(Ru,aCupu,'--b',Ru,aCupl,'--m',Ru,aAe,'k',Ru,aAeb,'g',Ru,aAew,'r');

plot(Ru,aCupu,'--b',Rb,aaCupl,'--m',Ru,aClowu,':r',Rb,aaClowl,':k',Ru,aAe,'k',Ru,aAeb,'g',Ru,aAew,'r');

legend('Capacity up upper bound','Capacity low up bound','Capacity up lower bound','Capacity low lower bound');
xlabel('Normalized Reuse Distance Ru');grid;
ylabel('ASE[Bits/Sec/Hz/Km^2]');
title('Figure 8 Effect of Shadowing');