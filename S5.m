clear all
%% parameters setting
IC = 600000; %Installed capacity of coal-fired power plants (kW)
RT = 5500; %Annual running time of coal-fired power plants (h) 
EF = 0.000893; %CO2 emission factor (t CO2/kWh)
CE = 0.9; %CO2 capture efficiency (%)
UTC = 100; %Unit transportation cost for CO2 captured (CNY/t CO2)
USC = 50; %Unit storage cost for CO2 captured (CNY/t CO2)
UnitI = 2666.67;
I0 = UnitI*IC; %Initial cost of CCS retrofitting investment
UnitOM = 62.9;
IOM0 = UnitOM*IC; %Initial cost of O&M 
Alpha = 0.0202; %Experience parameter for technological learning rate (%)
Beta = 0.057; %Experience parameter for O&M cost (%)
Sclean=0;%Clean electricity tariff (CNY/kWh) 
lamda = 0;
Pai = 0.94;%Unit efficiency (%)
Eta = 0.179;
rr=0.05;%discounted rate
N=35;
M=10000;
dt = 1;
D=20;
EP0=0.35; 
CP0=88; 
r2=0.04; 
sigma2=0.03; 

%% uncertainty simulation

EP = ones(M,N).*Sclean;%Clean electricity subsidy
CP = fun_GBM(dt,r2,sigma2,N,M,CP0);%carbon price

%% cash flows
Eoutput=Pai*IC*RT;%the amount of electricity production in year t
QCO2t = Eoutput*EF*CE; %the amount of CO2 captured in year t
TCt = UTC * QCO2t; %total transportation costs in year t
SCt = USC * QCO2t; %total transportation costs in year t
EPA = repmat(Eoutput,M,N);
CRA = repmat(QCO2t,M,N);
II = zeros(1,N);
I=zeros(M,N);
for i=1:N
    II(1,i) = I0 * exp (-Alpha*i)*(1-lamda);
end
I=repmat(II,M,1);
IOMM = zeros(1,N);
IOM = zeros(M,N);
for i=1:N
    IOMM(1,i) = IOM0 * exp (-Beta*i);
end
IOM = repmat(IOMM,M,1);

%% ROA
CFO = 0;
for k=1:D    
    EPU=0;
    CPU=0;
    EPAU=0;
    CRAU = 0;
    IU=0;
    IOMU=0;
    EGGU=0; 
    TCU = 0;
    SCU = 0;
    CR = 0;
    OR = 0;
    ER = 0;
    EL = 0;
    CF = 0;
    EPU=EP(1:M,k+1:N); 
    CPU=CP(1:M,k+1:N); 
    EPAU = EPA(1:M,k+1:N); 
    CRAU = CRA(1:M,k+1:N);
    IOMU = IOM(1:M,k+1:N); 
    IU=I(1:M,k); 
    TCU = repmat(TCt,M,N-k);
    SCU = repmat(SCt,M,N-k);
    ER = EPAU.*EPU; 
    CR=CPU.*CRAU;
    EL=Eta*EPAU.*EP0; 
    CF=ER+CR-EL-IOMU-TCU-SCU;
    for i=1:M  
       A=pvvar(CF(i,1:N-k),rr);
       CFO(i,k)=A;
    end
end

NPV0=exp(-rr).*CFO-I(1:M,1:D);
npv0=max(NPV0,0);
TO=zeros(M,2);

for i = 1:M
    A = 0 ;
    for j = 1:D
        if exp(-rr*(j-1))*npv0 (i,j) > A
            A = exp(-rr*(j-1))*npv0 (i,j);
            TO(i,1) = j;
            TO(i,2) = A;
        end
    end
end
Table = tabulate(TO(:,1));
[Opt_F,Opt_T] = max(Table(:,2));
Opt_DeT = Table(Opt_T,1);% optimal timing
count1=0;
count2=0;
for i = 1:M
    if TO(i,1)== Opt_DeT
        count1 = count1 + TO(i,2);
        count2 = count2 + NPV0(i,1);
    end
end
Opt_DeV = count1/Opt_F;% optimal value for deferred option
npv=count2/Opt_F;
