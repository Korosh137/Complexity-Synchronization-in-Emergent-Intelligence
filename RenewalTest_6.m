tic

% This code tests whether or not the extracted events from a signal (using stripes) are renewal (renewal experiment).
% Korosh Mahmoodi Nov. 2023
% Using this code cite:


ta = 100 ; % length of the age
StripeSize = 0.01 ; % Size of the stripes for extracting events from the signal (Data)

                    % For renewal experiment we need long time series. So, we ran the EI.m code with Trials = 1e7 (NS = 20.)
                            Data = DATA9x(1:1e7, 1)  ; % input data
Len = length(Data) ;

% Waiting-time probability density function of the:
P1 = zeros(Len, 1) ; %  events
P1Age = zeros(Len, 1) ; % aged events
P1SH = zeros(Len, 1) ; % shuffled events
P1SHAge = zeros(Len, 1) ; % shuffled-aged events


% Extracting the waiting-times (Tau) using stripes
Tau = zeros(Len,1);

Data = Data - min(Data) ;
Data = Data ./ max(Data) ;
RoundedData = round(Data./StripeSize, 0); 
Ta = zeros(1,Len) ;
k = 1 ;
Ta(1) = 1 ;

for i = 2 : Len
    if RoundedData(i) == RoundedData(i-1)
        Ta(k) = Ta(k) + 1 ;
    else
        k = k + 1 ;
        Ta(k) = 1 ;
    end
end

% sums the Taus
Ta = Ta(1:k) ;
Tau = Ta ;
                                                % RENEWAL EXPERIMENT

% Aging  Taus
LTau = zeros(k, 1);

LTau(1) = Tau(1) ;

for uuu = 2 : k
    LTau(uuu) = LTau(uuu-1) + Tau(uuu) ;
end

TauAge = zeros( 10 * k , 1) ;

c1Age = 0 ;
XX = 1  ;
while XX  >= 0
    c1Age = c1Age + 1 ;
    TauAge(c1Age) = Tau(c1Age) - ta ;

    if TauAge(c1Age) <= 0

        AAA = Tau(c1Age) ;
        hhh = c1Age ;
        while TauAge(c1Age) <= 0
            hhh = hhh + 1 ;
            AAA = AAA + Tau(hhh) ;

            TauAge(c1Age) = AAA - ta ;

            XX = LTau( k ) - (LTau(hhh) + ta ) ;
            if  XX < 0
                TauAge(c1Age) = 0 ;
                break
            end
        end

        XX = LTau( k ) - (LTau(c1Age) + ta ) ;

        if  XX < 0
            TauAge(c1Age) = 0 ;
            break
        end

    end
    XX = LTau(k ) - (LTau(c1Age) + ta ) ;
end

% Shuffling the  Taus
SH1 = k ;
TauSH = zeros(k ,1);    

for ttt = 1 : Len

    r1 = round((SH1 - 1)*rand + 1) ;  
    r2 = round((SH1 - 1)*rand + 1) ;

    AA = Tau(r1) ;
    BB = Tau(r2) ;

    TauSH(r1) = BB ;
    TauSH(r2) = AA ;

end

% Aging the shuffled Tau s

% Lenths
LTauSH = zeros( k, 1) ;
LTauSH(1) = TauSH(1) ;

for u1 = 2 :1:k
    LTauSH(u1) = LTauSH(u1-1) + TauSH(u1) ;
end


TauSHAge = zeros( 10*k , 1) ;

c1SHAge = 0 ;
XX = LTauSH( k ) - ta  ;

while XX  >= 0

    c1SHAge = c1SHAge + 1 ;
    TauSHAge(c1SHAge) = TauSH(c1SHAge) - ta ;

    if TauSHAge(c1SHAge) <= 0

        AAA = TauSH(c1SHAge) ;
        hhh = c1SHAge ;
        while TauSHAge(c1SHAge) <= 0
            hhh = hhh + 1 ;
            AAA = AAA + TauSH(hhh) ;

            TauSHAge(c1SHAge) = AAA - ta ;

            XX = LTauSH( k ) - (LTauSH(hhh) + ta ) ;
            if  XX < 0
                TauSHAge(c1SHAge) = 0 ;
                break
            end
        end

        XX = LTauSH( k ) - (LTauSH(c1SHAge) + ta ) ;
        if  XX < 0
            TauSHAge(c1SHAge) = 0 ;
            break
        end
    end
    XX = LTauSH( k ) - (LTauSH(c1SHAge) + ta ) ;
end

% Histograms
%  Tau
for  ii = 1 : 10000
    for  jj = 1 : k
        if Tau(jj) == ii
            P1(ii) = P1(ii) + 1;
        end
    end
end

%  Tau Age
for  ii = 1 : 10000
    for  jj = 1 : c1Age
        if TauAge(jj) == ii
            P1Age(ii) = P1Age(ii) + 1;
        end
    end
end

%  TauSH
for  ii = 1 : 10000
    for  jj = 1 : k
        if TauSH(jj) == ii
            P1SH(ii) = P1SH(ii) + 1;
        end
    end
end

%  TauSHAge
for  ii = 1 : 10000
    for  jj = 1 : c1SHAge
        if TauSHAge(jj) == ii
            P1SHAge(ii) = P1SHAge(ii) + 1;
        end
    end
end

% P1
 gg = (cumsum(P1));
  P1=P1/gg(Len') ;
% P1Age
 gg = (cumsum(P1Age));
  P1Age=P1Age/gg(Len') ;
% P1SH
 gg = (cumsum(P1SH));
  P1SH=P1SH/gg(Len') ;
% P1SHAge
 gg = (cumsum(P1SHAge));
  P1SHAge=P1SHAge/gg(Len') ;


A11 = P1(1:1e5, 1) ;
A22 = P1Age(1:1e5, 1) ;
A33 = P1SHAge(1:1e5, 1) ;
loglog(A11,'DisplayName','\psi(\tau)');hold on; loglog(A22,'DisplayName','Aged \psi(\tau)');hold on; loglog(A33,'DisplayName','Shuffled-Aged \psi(\tau)');hold off
xlabel('log(\tau)'), ylabel('log\psi(\tau)');
legend;

toc
