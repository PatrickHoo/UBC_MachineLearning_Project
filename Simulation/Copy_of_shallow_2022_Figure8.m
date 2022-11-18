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
sigmad = 16;
sigmaI = 16;
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
        
        mdd=Xi*log(udd);
        mii=Xi*log(uii);
        
        %log mean function for desired user simulation
        mudd = log((mdd^2)/sqrt(sigmad + mdd^2)); % rewrite mean power in log function
        muii = log((mii.^2)./sqrt(sigmad + mii.^2)); % rewrite mean power in log function
        
        %log variance function for interferer simulation
        sigmadd = sqrt(log(sigmad/(mdd^2)+1));
        sigmai = sqrt(log(sigmad./(mii.^2)+1));
        log_ud=lognrnd(mudd,sigmadd,[1,1]);
        log_ui=lognrnd(muii,sigmai,[6,1]);
        SII=sum(Si.*log_ui);
        gammadi=(Sd*log_ud)/SII;
        gammab=((((R*(Ru+1))/r)^a)*((g+(Ru+1)*R)/(g+r))^b)/NI;%the best case ofgammad
        
        %Calculate log parameters best-case
        udd_b=1/(r.^a.*(1+(r/g)).^b)*St;
        uii_b=1./((D+R)^a*(1+((D+R)/g))^b)*St;
        mdd_b=Xi*log(udd_b);
        mii_b=Xi*log(uii_b);
        %log function for desired user simulation
        mudd_b = log((mdd_b^2)/sqrt(sigmad+mdd_b^2));
        muii_b = log((mii_b^2)/sqrt(sigmad+mii_b^2));
        %log function for interferer simulation
        sigmadd_b = sqrt(log(sigmad/(mdd_b^2)+1));
        sigmai_b = sqrt(log(sigmad./(mii_b^2)+1));
        log_udb=lognrnd(mudd_b,sigmadd_b,[1,1]);
        log_uib=lognrnd(muii_b,sigmai_b,[6,1]);
        
        SIIb=sum(Sib*log_uib);
        gammadib=(Sd*log_udb)/SIIb;

      % gammadib = log_udb/sum(log_uib);


        %Calculate log parameters best-case
        udd_w=1/(r.^a.*(1+(r/g)).^b)*St;
        uii_w=1./((D-R)^a*(1+((D-R)/g))^b)*St;
        mdd_w=Xi*log(udd_w);
        mii_w=Xi*log(uii_w);
        %log function for desired user simulation
        mudd_w = log((mdd_w^2)/sqrt(sigmad+mdd_w^2));
        muii_w = log((mii_w^2)/sqrt(sigmad+mii_w^2));
        %log function for interferer simulation
        sigmadd_w = sqrt(log(sigmad/(mdd_w^2)+1));
        sigmai_w = sqrt(log(sigmad/(mii_w^2)+1));
        
        log_udw=lognrnd(mudd_w,sigmadd_w,[1,1]);
        log_uiw=lognrnd(muii_w,sigmai_w,[6,1]);
        
        SIIw=sum(Siw*log_uiw);
        gammadiw=(Sd*log_udw)/SIIw;

      %  gammadiw = log_udw/sum(log_uiw);
        
        
        
        %Calculate sigmaSI
        sigmaSI=(Xi^2)*(log((NI-1+exp(sigmad/(Xi^2)))/NI));%log mean variance of interference
        sigmard=sigmad+sigmaSI; %log variance of desired user
        Urdu=Xi.*log(((R*(Ru+1)/r)^a).*((g+(Ru+1)*R)/(g+r)).^b)-Xi*log(NI)+(sigmaSI-sigmaI)/(2*Xi);%the upper log mean
        Urdl=Xi.*log(((R*(Ru-1)/r).^a).*((g+(Ru-1)*R)/(g+r)).^b)-Xi*log(NI)+(sigmaSI-sigmaI)/(2*Xi);%the lower log mean
        
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