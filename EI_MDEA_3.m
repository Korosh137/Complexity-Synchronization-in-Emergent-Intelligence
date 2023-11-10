tic

% This code evaluates the MDEA on the data generated from EI.m code.

CHANNELS = 2 ; % two time series for MDEA
str = 0.001 ; %  size of the stripes
ST =  0.4 ;  % ST and EN are the start and end points on the diffusion entropy graph which 
EN =  0.6 ;  % the slope of its linear fit gives the scaling delta
Rule = 1 ;   % Rule for MDEA
PLOT = 14 ;  % if 1, it plots the diffusion entropy graph (Note that there will be many graphs, so you should define in the code which graphs you want to see)

Slice = 3e4 ; % size of the slice of data for evaluating its scaling
Correl= zeros( SizeIndex ,   10  ) ; % correlation between the scaling time series of the two selected signals
MeanDelF= zeros( SizeIndex ,  10  ) ; % average of the scaling time series for agent F
MeanDelS= zeros( SizeIndex ,  10  ) ;
VarDelF= zeros( SizeIndex ,   10  ) ; % variance of the scalings evaluated from ten experiments for agent F
VarDelS= zeros( SizeIndex ,   10  ) ;

for kk = 1 :  SizeIndex  
kk
    for jj = 1 : 10   % 10 experiment for each N
            
    % select two time series which you want to evaluate their scaling time series using MDEA
        data(:,1) =  DATAall( kk, 1e4:1e6 , 5 , jj) ; % firs time series for DEA
        data(:,2) =  DATAall( kk, 1e4:1e6 , 6 , jj) ; % second time series for DEA

        nn =  floor(  ( length(data(:, 1)) - 0)/(Slice/3)   ) -  2 ; % number of slices for DEA (sliding windows with overlap of 2/3)

        Deltann = zeros(nn , CHANNELS) ; % records the scaling time series of the two time series

        for gg11 = 1 :   nn
            sta =  floor((gg11 -1)*(Slice/3) ) ;
            DaTaa1 = zeros(Slice, 1) ;
            DaTaa2 = zeros(Slice, 1) ;

            for yytt  = 1 : Slice
                DaTaa1(yytt) = data( yytt  + sta, 1 ) ;
                DaTaa2(yytt) = data( yytt  + sta, 2 ) ;
            end
            Deltann(gg11, 1) = MDEA(DaTaa1, str, Rule, ST, EN, PLOT , gg11)  ;
            Deltann(gg11, 2) = MDEA(DaTaa2, str, Rule, ST, EN, PLOT , gg11)  ;
        end

        XX = Deltann(:, 1) ;
        YY = Deltann(:, 2) ;

        [xc, pxc, rlo, rup] = corrcoef(Deltann);
        Correl(kk , jj) =  xc(1,2)  ;

        MeanDelF(kk , jj) = mean(XX) ;
        MeanDelS(kk , jj) = mean(YY) ;

        VarDelF(kk , jj) = var(XX) ;
        VarDelS(kk , jj) = var(YY) ;

    end
end

Meancorr = mean(Correl') ;

plot(Correl,'DisplayName','Correl')

meanmeanScalingF = mean(MeanDelF') ;
meanmeanScalingS = mean(MeanDelS') ;


toc