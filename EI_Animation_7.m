% Korosh Mahmoodi
% 05/18/2023
% funccc
% There are two types of agents, named 'F' and 'S', as metaphor for Fish and Shark
% Cite the corresponding paper:

tic
clc
clear all
close all

mov = VideoWriter('CS.avi');   % Creates CS.avi file of F agents as circles and S agents by circles
% Red circles are the F agents which copy their neighbors with chance >= 0.3.
% Green circles are the F agents which make their own decision
open(mov)
PLOT = 1 ; % If = 1, creates and shows the animation

Trials = 1e4 ; % Number of trials

    deltaF = 0.2 ; % Constant relating payoff to the change of the threshold
    deltaS = 0.2 ;
    L =  1 ;  % Domain size
    dt =  1 ; % Time step
    %%%%%%%%%%%%%%% Parameters for agents S
    NS =  50 ;    % Number of agents S
    rS = 0.15 ;   % Vision radius of agents S
    velS = 0.01 ; % Magnitude of the velocity of agents S
    %%%%%%%%%%%%%%% Parameters for agents F
    NF =   NS ;   % Number of agents F
    rF = rS/1 ;   % Radius of vision of agent F
    velF = velS/1 ; % Velocity of agents F
    rG = 1*rF ;   % Radius in which agents F and S can play
    %%%%%%%% Initial conditions
    thetaF=2*pi*(rand(1,NF)-0.5); % Angles of agents F
    thetaS=2*pi*(rand(1,NS)-0.5); % Angles of agents S
    thetaFS = 2*pi*(rand(1,NF+NS)-0.5);  % Matrix of all the angles
    PayF0 = zeros(NF , 1) ;  % previous payoffs of agents F
    PayF = zeros(NF , 1) ;   % current payoffs of agents F
    PayS0 = zeros(NS , 1) ;  
    PayS = zeros(NS , 1) ;   

    RelianceinfoF = 0.5*ones(NF, 1)  ;   % Threshold for agent F to use information of its neighboring F
    RelianceinfoF_01 =  round(rand(NF, 1) ) ;     % Records whether agent F relied (1) or not (0)   
    RelianceinfoS = 0.5*ones(NS, 1)  ; 
    RelianceinfoS_01 = round(rand(NS, 1) ) ; 

    OF_FS = 0.5*ones(NF, 1)  ;           % Threshold for agent F to Follow or Oppose the predicted position of S agents
    OF_FS_01 =  round(rand(NF, 1) ) ;    % Records whether agent F Followed or Opposed
    OF_SF = 0.5*ones(NS, 1)  ; 
    OF_SF_01 =  round(rand(NS, 1) ) ; 

    aveRelianceinfoF = zeros(Trials, 1) ; % Average of the RelianceinfoF threshold
    aveRelianceinfoS = zeros(Trials, 1) ;

    aveOF_FS = zeros(Trials, 1) ;        % Average of the OF_FS threshold
    aveOF_SF = zeros(Trials, 1) ;

    thetaAllyF = zeros(NF , 1) ;         % Angle for agent F, evaluated from neighboring F agents
    thetaFoeF = zeros(NF , 1) ;          % Angle for agent F, evaluated from neighboring S agents

    thetaAllyS = zeros(NF , 1) ;
    thetaFoeS = zeros(NS , 1) ;

    AnyF= zeros(NF , 1) ;               % Records '13' if there isn't Any F in neighberhood
    AnyS= zeros(NS , 1) ;

    OrderF = zeros(Trials, 1) ;  % Order parameter for F
    OrderS = zeros(Trials, 1) ;  

    IndividualFVx = zeros(Trials, 1) ; % X component of the velocity of one of the F agents
    IndividualSVx = zeros(Trials, 1) ;
    IndividualFX = zeros(Trials, 1) ; % Position on x axis of the one of the F agents
    IndividualSX = zeros(Trials, 1) ;

    xF=L*rand(1,NF);  % x position of agents F
    yF=L*rand(1,NF);  % y position of agents F
    xS=L*rand(1,NS);  
    yS=L*rand(1,NS);  

    xFS = zeros(1, NS + NF) ; % x position of agents F and S (saved in one matrix)
    yFS = zeros(1, NS + NF) ;

    PayttF = zeros(Trials, 1) ; % Average payoff of F agents at each trial
    PayttS = zeros(Trials, 1) ;

    for tt= 2 : Trials  %%% tt is the trial number

        for ii = 1 : NF             
            xFS(1, ii) = xF(ii) ;
            yFS(1, ii) = yF(ii) ;
            thetaFS(ii) = thetaF(ii) ;
        end
        for ii = NF +1 : NF + NS
            xFS(1, ii) = xS(1, ii - NF) ;
            yFS(1, ii) = yS(1, ii - NF) ;
            thetaFS(ii) = thetaS(ii- NF) ;
        end

  % Decision (angle), if agent uses the information of the agents of other group
        [l1FS,l2FS]= Finddistance(xFS, yFS, rF, L) ; % S agents in vision of F agents

        for i = 1:NF        
            list = l1FS(l2FS==i);
            list = list(list > NF  ) ;

            if ~isempty(list)
                                  % predicting the mean location of S agents
                xSthatFpredicts = mean(xFS(list)) + velS*mean(cos(thetaFS(list)))*dt  ;
                ySthatFpredicts = mean(yFS(list)) + velS*mean(sin(thetaFS(list)))*dt  ;

                r = rand ;
                if r > OF_FS(i)  % Decision to go Away or Toward the predicted position

                    thetaFoeF(i) =   0 + atan2( (ySthatFpredicts - yFS(i)) , (xSthatFpredicts - xFS(i)) )  ;
                    OF_FS_01(i) = 1 ;
                else
                    thetaFoeF(i) =   pi + atan2( (ySthatFpredicts - yFS(i)) , (xSthatFpredicts - xFS(i)) )  ;
                    OF_FS_01(i) = 0 ;
                end
            else
                thetaFoeF(i) = thetaFS(i) ;
                OF_FS_01(i) = -1 ;
            end
        end

        [l1FS,l2FS]= Finddistance(xFS, yFS, rS, L) ; %%% F in vison of S

        for i = 1 + NF : NS + NF
            list = l1FS(l2FS==i);
            list = list(list <= NF  ) ;

            if  ~isempty(list)
                                  % predicting the mean location of F agents
                xFthatSpredicts = mean(xFS(list)) + velF*mean(cos(thetaFS(list)))*dt;
                yFthatSpredicts = mean(yFS(list)) +  velF*mean(sin(thetaFS(list)))*dt;

                r = rand ;
                if r > OF_SF(i- NF)  % Decision to go Away or Toward the predicted position

                    thetaFoeS(i - NF) =  atan2( (yFthatSpredicts - yFS(i)) , (xFthatSpredicts - xFS(i)) ) ;
                    OF_SF_01(i- NF) = 1 ;
                else
                    thetaFoeS(i - NF) =  pi+ atan2( (yFthatSpredicts - yFS(i)) , (xFthatSpredicts - xFS(i)) ) ;
                    OF_SF_01(i- NF) = 0 ;
                end
            else
                thetaFoeS(i-NF) = thetaFS(i) ;
                OF_SF_01(i- NF) = -1 ;
            end
        end

        % Calculation of average angle, if agent uses information of the same group
        [l1F,l2F]= Finddistance(xF, yF, rF, L) ;

        for i = 1:NF
            list = l1F(l2F==i);
            if ~isempty(list)
                thetaAllyF(i) = atan2(mean(sin(thetaF(list))),mean(cos(thetaF(list))));
            else
                thetaAllyF(i) = thetaF(i) ; 
                AnyF(i) = 13 ;
            end
        end

        for ii = 1 : NF
            r00 = rand ;
            if  r00            >     RelianceinfoF(ii)   
                thetaF(ii) = thetaAllyF(ii);
                RelianceinfoF_01(ii) = 1 ;
            else
                thetaF(ii) = thetaFoeF(ii) ;
                RelianceinfoF_01(ii) = 0 ;
            end
        end

        [l1,l2]= Finddistance(xS, yS, rS, L) ;
        for i = 1:NS
            list = l1(l2==i);
            if ~isempty(list)
                thetaAllyS(i) = atan2(mean(sin(thetaS(list))),mean(cos(thetaS(list))));
            else
                thetaAllyS(i) =  thetaS(i) ;
            end
        end

        for ii = 1 : NS
            r00 = rand ;
            if  r00          >  RelianceinfoS(ii) 
                thetaS(ii) = thetaAllyS(ii); 
                RelianceinfoS_01(ii) = 1 ;
            else
                thetaS(ii) = thetaFoeS(ii) ;
                RelianceinfoS_01(ii) = 0 ;
            end
        end

        % Update
        xF = xF + velF*cos(thetaF)*dt;
        yF = yF + velF*sin(thetaF)*dt;
        % Periodic condition
        xF(xF<0) = L + xF(xF<0);
        xF(L<xF) = xF(L<xF) - L;
        yF(yF<0) = L + yF(yF<0);
        yF(L<yF) = yF(L<yF) - L;

        xS = xS + velS*cos(thetaS)*dt;
        yS = yS + velS*sin(thetaS)*dt;
        
        xS(xS<0) = L + xS(xS<0);
        xS(L<xS) = xS(L<xS) - L;
        yS(yS<0) = L + yS(yS<0);
        yS(L<yS) = yS(L<yS) - L;

        OrderF(tt) = ( (mean(cos(thetaF)))^2 + (mean(sin(thetaF)))^2  )^0.5  ;
        OrderS(tt) = ( (mean(cos(thetaS)))^2 + (mean(sin(thetaS)))^2  )^0.5  ;

        IndividualFVx(tt) = cos(thetaF(1)) ;
        IndividualSVx(tt) = cos(thetaS(1)) ;
        IndividualFX(tt) = xF(1) ;
        IndividualSX(tt) = xS(1) ;

        % payoff 
        xFS = zeros(1, NS + NF) ;
        yFS = zeros(1, NS + NF) ;

        for ii = 1 : NF             
            xFS(1, ii) = xF(ii) ;
            yFS(1, ii) = yF(ii) ;
        end
        for ii = NF +1 : NF + NS
            xFS(1, ii) = xS(1, ii - NF) ;
            yFS(1, ii) = yS(1, ii - NF) ;
        end

        [l1FS,l2FS]= Finddistance(xFS, yFS, rG, L) ;
        for i = 1:NF        % payoff of F
            list = l1FS(l2FS==i);
            listofS = list(list > NF  ) ;  % # of S in rG neighborihood
            LS= length(listofS) ;
           
            if LS == 0
                PayF(i) =  1 ; 
            else
                PayF(i) =  (-1)*LS ; 
            end
        end

        for i = 1 + NF : NS + NF        % payoff of S
            list = l1FS(l2FS==i);
            listofF = list(list <= NF ) ; % # of F in rG neighborihood
            LF= length(listofF) ;
      
            if LF == 0
                PayS(i-NF) = -1 ; 
            else
                PayS(i-NF) =  (1)*LF ;
            end
        end

        % updating Thresholds

        for jj = 1 : NF
            RelianceinfoF(jj) = UpdateThreshold(RelianceinfoF(jj) , RelianceinfoF_01(jj)  , PayF(jj) , PayF0(jj) , deltaF) ;

            if  RelianceinfoF_01(jj) == 0  && OF_FS_01(jj) ~= -1
                OF_FS(jj) = UpdateThreshold(OF_FS(jj) , OF_FS_01(jj) , PayF(jj) , PayF0(jj) , deltaF) ;
            end
        end

        for jj = 1 : NS
            RelianceinfoS(jj) = UpdateThreshold(RelianceinfoS(jj) ,     RelianceinfoS_01(jj)  , PayS(jj) , PayS0(jj) , deltaS) ;

            if  RelianceinfoS_01(jj) == 0  && OF_SF_01(jj) ~= -1
                OF_SF(jj) = UpdateThreshold(OF_SF(jj) , OF_SF_01(jj) , PayS(jj) , PayS0(jj) , deltaS) ;
            end
        end

        aveRelianceinfoF(tt) = mean(RelianceinfoF) ;
        aveRelianceinfoS(tt) = mean(RelianceinfoS) ;

        PayttF(tt) = mean(PayF) ;
        PayttS(tt) = mean(PayS) ;

        aveOF_FS(tt) = mean(OF_FS) ;
        aveOF_SF(tt) = mean(OF_SF) ;

        PayF0 = PayF;
        PayS0 = PayS;



    if PLOT==1  && tt/ Trials > 0.9
        %%%%%%%%%%% plotinig
        XFlead = zeros(1, NF) ;
        YFlead = zeros(1, NF) ;

        for yy = 1 : NF
            if RelianceinfoF(yy)           <= 0.25  %% Solo  %%%
                XFlead(yy) = xF(yy) ;
                YFlead(yy) = yF(yy) ;

            end

        end
        XFlead = XFlead(XFlead~=0);
        YFlead = YFlead(YFlead~=0);


        %%%%%%%% Figure
        plot(xF,yF,'r.','MarkerSize',15)
        xlim([0 L]);
        ylim([0 L]);
        axis square

        %     pause(0.1)

        hold on ;


        %%%% leaders
        plot(XFlead,YFlead,'g.','MarkerSize',15)
        xlim([0 L]);
        ylim([0 L]);
        axis square

        hold on ;

        %%%%%%%% Figure Shark
        plot(xS,yS,'bv','MarkerSize',7)
        xlim([0 L]);
        ylim([0 L]);
        axis square

        pause(0.1)

        hold off ;


        frame = getframe(gcf);
        writeVideo(mov,frame)

    end


    end


close(mov)
toc


function [A,B] = Finddistance(x, y, r, L)

D = pdist([x' y'],'euclidean') ;

% Periodic boundary 
x0( x < r ) = L + x( x < r ) ;
x0( x > L-r ) = x( x > L-r ) -L ;
x0( r <= x & x <= L-r ) = x(r <= x & x <= L-r ) ;
y0( y < r) = L + y( y < r ) ;
y0( y > L-r) = y( y > L-r )-L ;
y0( r <= y & y <= L-r ) = y( r<= y & y <= L-r) ;

D00 = pdist([x0' y0'],'euclidean') ;
D = min([D; D00]) ;

M = squareform(D) ; 
[A,B]= find( 0 < M & M < r  ) ;  

end


function aa = UpdateThreshold(pi0 , Decision0 , Pay , Paybefore , ChangeThreshold)

if Decision0 == 1   
    dp = -1*ChangeThreshold ;
else
    dp =  1*ChangeThreshold ;
end

DeltaPay =  (Pay - Paybefore)  / (abs(Pay) + abs(Paybefore) ) ;
aa =  pi0  + dp * DeltaPay ;

% boundary conditions
if  aa  < 0
    aa = 0 ;
end

if aa  > 1
    aa  = 1 ;
end

end
