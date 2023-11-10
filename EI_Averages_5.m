% This code evaluates the average of the payoffs

PAYmatixF = zeros(SizeIndex, 10) ;
PAYmatixS = zeros(SizeIndex, 10) ;


for hh = 1 :    SizeIndex
    for mn = 1 :10  % ensemble
        PAYmatixF(hh, mn) =  mean(DATAall(hh, 1e4:1e6, 3, mn)) ;  
        PAYmatixS(hh, mn) =  mean(DATAall(hh, 1e4:1e6, 4, mn)) ; 
    end
end

meanPayMatrixF = mean(PAYmatixF') ;
meanPayMatrixS = mean(PAYmatixS') ;

summmm = meanPayMatrixF + meanPayMatrixS ;

SUMMAtPay = PAYmatixF + PAYmatixS ;


