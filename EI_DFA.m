
% This code evaluates the DFA on the data generated from EI.m code.

tic

Order = 1 ; % Order of detrend in DFA
PLOT = 18 ; % if 1, it plots the diffusion entropy graph (Note that there will be many graphs, so you should define in the code which graphs you want to see)

             Maxplot = 1 ;  % Maximum number of DFA graphs to be ploted
             pts = 10:5:100 ; % length of linear fit of the DFA graph

CHANNELS = 2 ; % two time series for DFA

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
        data(:,1) =  diff( DATAall( kk, 1e4:1e6 , 5 , jj) ); % firs time series for DEA
        data(:,2) =  diff( DATAall( kk, 1e4:1e6 , 6 , jj) ); % second time series for DEA

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

            [A1,F1] = DFA_fun(DaTaa1 , pts , Order , PLOT,  Maxplot ); 
            [A2,F2] = DFA_fun(DaTaa2 , pts , Order , PLOT,  Maxplot ); 

            Deltann(gg11, 1) = A1(1,1);
            Deltann(gg11, 2) = A2(1,1) ;

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



%%%  https://www.mathworks.com/matlabcentral/fileexchange/67889-detrended-fluctuation-analysis-dfa
function[A,F] = DFA_fun(data, pts, order, PLOT, maxplot)
% -----------------------------------------------------
% DESCRIPTION:
% Function for the DFA analysis.
% INPUTS: 
% data: a one-dimensional data vector.
% pts: sizes of the windows/bins at which to evaluate the fluctuation
% order: (optional) order of the polynomial for the local trend correction.
% if not specified, order == 1;
% OUTPUTS: 
% A: a 2x1 vector. A(1) is the scaling coefficient "alpha",
% A(2) the intercept of the log-log regression, useful for plotting (see examples).
% F: A vector of size Nx1 containing the fluctuations corresponding to the
% windows specified in entries in pts.
% -----------------------------------------------------
% Checking the inputs
if nargin == 4
   order = 1; 
end
sz = size(data);
if sz(1)< sz(2)
    data = data';
end
exit = 0;
if min(pts) == order+1
        disp(['WARNING: The smallest window size is ' num2str(min(pts)) '. DFA order is ' num2str(order) '.'])
        disp('This severly affects the estimate of the scaling coefficient')
        disp('(If order == [] (so 1), the corresponding fluctuation is zero.)')
elseif min(pts) < (order+1)
        disp(['ERROR: The smallest window size is ' num2str(min(pts)) '. DFA order is ' num2str(order) ':'])
        disp(['Aborting. The smallest window size should be of ' num2str(order+1) ' points at least.'])
        exit = 1;
end
if exit == 1
    return
end
% DFA
npts = numel(pts);
F = zeros(npts,1);
N = length(data);
for h = 1:npts
    
    w = pts(h);
    
    n = floor(N/w);
    Nfloor = n*pts(h);
    D = data(1:Nfloor);
    
    y  = cumsum(D-mean(D)); %  (D-mean(D));
    
    bin = 0:w:(Nfloor-1);
    vec = 1:w;
    
    coeff = arrayfun(@(j) polyfit(vec',y(bin(j) + vec),order),1:n,'uni',0);
    y_hat = cell2mat(cellfun(@(y) polyval(y,vec),coeff,'uni',0));
    F(h)  = mean((y - y_hat').^2)^0.5;
    
end

A = polyfit(log10(pts),log10(F)',1);

if PLOT == 1  &&  maxplot < 5
figure;
scatter(log10(pts),log10(F))
plot_fun = @(xp,A,ord) polyval(A,log10(xp));

hold on
x =  pts ;
plot(log10(x),plot_fun(x,A),'--')
xlabel('log_{10} W'), ylabel('log_{10} F(W)');
legend(['\alpha = ' num2str( sprintf('%.3f',A(1) ))],'Location','northwest');

hold off

end



end




