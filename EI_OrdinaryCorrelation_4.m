
% This code evaluates cross-correlation between the two time sereis

tic

CHANNELS = 2 ; % two time series for MDEA
corrMeanOrdinary = zeros(SizeIndex , 10) ;
Slice = 3e4 ; % size of the slice of data for evaluating its scaling

for kk = 1 : SizeIndex
    for jj = 1 : 10  % for each N there are 10 simulations

        data = zeros(1e6-1e4+1, CHANNELS) ;

            % select two time series which you want to evaluate their scaling time series using MDEA
        data(:,1) = DATAall( kk, 1e4:1e6 , 5 , jj) ; 
        data(:,2) = DATAall( kk, 1e4:1e6 , 6 , jj) ; 

        nn =  floor(  ( length(data(:, 1)) - 0)/(Slice/3)   ) -  2 ;

        Deltann = zeros(nn , CHANNELS) ;

        TimeSlice = zeros(nn , 1) ;
        for hh11 = 1 :   nn
            TimeSlice(hh11) = (Slice/2) +  (hh11-1)*(Slice/3) ;
        end

        Correl = zeros(nn, 1) ;

        for gg11 = 1 :   nn
            %%% each slice of data has 2/3 of the previous one
            sta =  floor((gg11 -1)*(Slice/3) ) ;
            DaTaa1 = zeros(Slice, 1) ;
            DaTaa2 = zeros(Slice, 1) ;

            for yytt  = 1 : Slice
                DaTaa1(yytt) = data( yytt  + sta, 1 ) ;
                DaTaa2(yytt) = data( yytt  + sta, 2 ) ;
            end

            XXXX = zeros(Slice, 2) ;
            XXXX(:, 1) = DaTaa1 ;
            XXXX(:, 2) = DaTaa2 ;

            [xc, pxc, rlo, rup] = corrcoef(XXXX);
            Correl(gg11) =  xc(1,2)  ;

        end
        corrMeanOrdinary( kk,  jj) = mean(Correl) ;
    end

end  

Meancorr = mean(corrMeanOrdinary') ;

toc
