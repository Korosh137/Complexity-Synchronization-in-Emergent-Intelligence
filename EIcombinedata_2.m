%  upload the data of simulations of different size N=NS=NF from EI.m code and using this code 
% you can combine them in a single matrix DATAall which is the input for the rest of the analysis

    SizeIndex = 9 ; % nine experiments with different size N

    DATAall = zeros(SizeIndex , 1e6, 12, 10) ;
    DATAall(1, :, :, :) = DATAA10(: , : , :) ; 
    DATAall(2, :, :, :) = DATAA15(: , : , :) ;
    DATAall(3, :, :, :) = DATAA20(: , : , :) ;
    DATAall(4, :, :, :) = DATAA25(: , : , :) ;
    DATAall(5, :, :, :) = DATAA30(: , : , :) ;
    DATAall(6, :, :, :) = DATAA35(: , : , :) ;
    DATAall(7, :, :, :) = DATAA40(: , : , :) ;
    DATAall(8, :, :, :) = DATAA45(: , : , :) ;
    DATAall(9, :, :, :) = DATAA50(: , : , :) ;