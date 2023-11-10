tic

% This code test thether the extracted events, using stripes, are renewal or not (renewal experiment).


ta = 100 ; % length of the age
SmaLen = 10000 ;
StripeSize = 0.01 ;
                    % For renewal experiment we need long time series. So, we ran the EI.m code with Trials = 1e7 and for NS = 20.
                            Data = DATA9x(1:1e7, 1)  ; 
Len = length(Data) ;
P1 = zeros(Len, 1) ;
P1Age = zeros(Len, 1) ;
P1SH = zeros(Len, 1) ;
P1SHAge = zeros(Len, 1) ;
DicatomX2 = zeros(Len, 1) ;
DPhiY1 = zeros(Len,1) ;
DPhiYn1 = zeros(Len,1) ;

% FINDING tavs
tav1 = zeros(Len,1);

Data = Data - min(Data) ;
Data = Data ./ max(Data) ;
RoundedData = round(Data./StripeSize, 0); 
Tau = zeros(1,Len) ;
k = 1 ;
Tau(1) = 1 ;

for i = 2 : Len
    if RoundedData(i) == RoundedData(i-1)
        Tau(k) = Tau(k) + 1 ;
    else
        k = k + 1 ;
        Tau(k) = 1 ;
    end
end

% sums the Taus
Tau = Tau(1:k) ;
tav1 = Tau ;
                                                % RENEWAL EXPERIMENT

% Aging  tav1
Ltav1 = zeros(k, 1);

Ltav1(1) = tav1(1) ;

for uuu = 2 : k
    Ltav1(uuu) = Ltav1(uuu-1) + tav1(uuu) ;
end

tav1Age = zeros( 10 * k , 1) ;

% After first event
c1Age = 0 ;
XX = 1  ;
while XX  >= 0
    c1Age = c1Age + 1 ;
    tav1Age(c1Age) = tav1(c1Age) - ta ;

    if tav1Age(c1Age) <= 0

        AAA = tav1(c1Age) ;
        hhh = c1Age ;
        while tav1Age(c1Age) <= 0
            hhh = hhh + 1 ;
            AAA = AAA + tav1(hhh) ;

            tav1Age(c1Age) = AAA - ta ;

            XX = Ltav1( k ) - (Ltav1(hhh) + ta ) ;
            if  XX < 0
                tav1Age(c1Age) = 0 ;
                break
            end
        end

        XX = Ltav1( k ) - (Ltav1(c1Age) + ta ) ;

        if  XX < 0
            tav1Age(c1Age) = 0 ;
            break
        end

    end
    XX = Ltav1(k ) - (Ltav1(c1Age) + ta ) ;
end

% Shuffle of  tav1 s
SH1 = k ;
tav1SH = zeros(k ,1);    

for ttt = 1 : Len

    r1 = round((SH1 - 1)*rand + 1) ;  
    r2 = round((SH1 - 1)*rand + 1) ;

    AA = tav1(r1) ;
    BB = tav1(r2) ;

    tav1SH(r1) = BB ;
    tav1SH(r2) = AA ;

end

% Aging the shuffled tav1 s

% Lenths
Ltav1SH = zeros( k, 1) ;
Ltav1SH(1) = tav1SH(1) ;

for u1 = 2 :1:k
    Ltav1SH(u1) = Ltav1SH(u1-1) + tav1SH(u1) ;
end


tav1SHAge = zeros( 10*k , 1) ;

c1SHAge = 0 ;
XX = Ltav1SH( k ) - ta  ;

while XX  >= 0

    c1SHAge = c1SHAge + 1 ;
    tav1SHAge(c1SHAge) = tav1SH(c1SHAge) - ta ;

    if tav1SHAge(c1SHAge) <= 0

        AAA = tav1SH(c1SHAge) ;
        hhh = c1SHAge ;
        while tav1SHAge(c1SHAge) <= 0
            hhh = hhh + 1 ;
            AAA = AAA + tav1SH(hhh) ;

            tav1SHAge(c1SHAge) = AAA - ta ;

            XX = Ltav1SH( k ) - (Ltav1SH(hhh) + ta ) ;
            if  XX < 0
                tav1SHAge(c1SHAge) = 0 ;
                break
            end
        end

        XX = Ltav1SH( k ) - (Ltav1SH(c1SHAge) + ta ) ;
        if  XX < 0
            tav1SHAge(c1SHAge) = 0 ;
            break
        end
    end
    XX = Ltav1SH( k ) - (Ltav1SH(c1SHAge) + ta ) ;
end

% Histograms
%  tav1
for  ii = 1 : 10000
    for  jj = 1 : k
        if tav1(jj) == ii
            P1(ii) = P1(ii) + 1;
        end
    end
end

%  tav1 Age
for  ii = 1 : 10000
    for  jj = 1 : c1Age
        if tav1Age(jj) == ii
            P1Age(ii) = P1Age(ii) + 1;
        end
    end
end

%  tav1SH
for  ii = 1 : 10000
    for  jj = 1 : k
        if tav1SH(jj) == ii
            P1SH(ii) = P1SH(ii) + 1;
        end
    end
end

%  tav1SHAge
for  ii = 1 : 10000
    for  jj = 1 : c1SHAge
        if tav1SHAge(jj) == ii
            P1SHAge(ii) = P1SHAge(ii) + 1;
        end
    end
end

% P1
for y = 1 : Len
    P1(y) = P1(y)/1 ;
end

Sum = 0 ;
for o = 1: Len
    Sum = Sum + P1(o) ;
end

for zz = 1: SmaLen
    P1(zz) = P1(zz)/ Sum ;
end

SmaLenP1 = zeros(SmaLen - 1, 1);
for o = 1: SmaLen
    SmaLenP1(o) = P1(o);
end

% P1Age
for y = 1 : Len
    P1Age(y) = P1Age(y)/1 ;
end

Sum = 0 ;
for o = 1: Len
    Sum = Sum + P1Age(o) ;
end

for zz = 1: SmaLen
    P1Age(zz) = P1Age(zz)/ Sum ;
end

SmaLenP1Age = zeros(SmaLen - 1, 1);

for o = 1: SmaLen
    SmaLenP1Age(o) = P1Age(o);
end

% P1SH
for y = 1 : Len
    P1SH(y) = P1SH(y)/1 ;
end

Sum = 0 ;
for o = 1: Len
    Sum = Sum + P1SH(o) ;
end

for zz = 1: SmaLen
    P1SH(zz) = P1SH(zz)/ Sum ;
end

SmaLenP1SH = zeros(SmaLen - 1, 1);

for o = 1: SmaLen
    SmaLenP1SH(o) = P1SH(o);
end

% P1SHAge
for y = 1 : Len
    P1SHAge(y) = P1SHAge(y)/1 ;
end

Sum = 0 ;
for o = 1: Len
    Sum = Sum + P1SHAge(o) ;
end

for zz = 1: SmaLen
    P1SHAge(zz) = P1SHAge(zz)/ Sum ;
end


A11 = P1(1:1e5, 1) ;
A22 = P1Age(1:1e5, 1) ;
A33 = P1SHAge(1:1e5, 1) ;
loglog(A11,'DisplayName','\psi(\tau)');hold on; loglog(A22,'DisplayName','Aged \psi(\tau)');hold on; loglog(A33,'DisplayName','Shuffled-Aged \psi(\tau)');hold off
xlabel('log(\tau)'), ylabel('log\psi(\tau)');


toc