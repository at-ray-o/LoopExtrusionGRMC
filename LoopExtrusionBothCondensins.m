%% Load the cdf_inv of actual distribution
delete('LoopListsCond1/*')
delete('LoopListsCond2/*')
cdf_inv_pofl = load("pofL_cdf_inv.mat");
cdf_inv_pofl = cdf_inv_pofl.dummyVar;
delX = cdf_inv_pofl(2,1)-cdf_inv_pofl(1,1);
%%


%Intialise Condensin at different sites
%rng(445)
monomerSize = 100; %bps
delT = 1; %seconds
pauseTime = 7; %seconds
residenceTimeRealCond1 = 2*60; %seconds
residenceTimeRealCond2 = 6*60; %seconds
nSites = 1000000; %Each monomer is 100 bp
probSwitch = 1;

nCond1Total = 2550;
nCond1Bound = 1210;

nCond2Total = 343;
nCond2Bound = 211;

extrusionRateReal = 1500; %bp/s
extrusionRateSim = extrusionRateReal/monomerSize;
residenceTimeSimCond1 = residenceTimeRealCond1/delT;
residenceTimeSimCond2 = residenceTimeRealCond2/delT;
%stepSzLenOnce=extrusionRateSim;

nCond1Unbound = nCond1Total-nCond1Bound;
nCond2Unbound = nCond2Total-nCond2Bound;

probCond1Ini = nCond1Bound/nSites;
probOffCond1 = 1/(residenceTimeSimCond1);
probOnCond1 = nCond1Bound*probOffCond1/nCond1Unbound;

probOffCond2 = 1/(residenceTimeSimCond2);
probOnCond2 = nCond2Bound*probOffCond2/nCond2Unbound;

%Simulation parameters
freq = 360;
tMaxSim = 60*60/delT;
avgLoopLen = zeros(tMaxSim,3);
runMax = 50;
dum = 0;
ndum = 0;
%flBound = fopen('boundTime.txt','a');
for run = 1:runMax
    head1LocsCond1 = zeros(nSites,1); %static head
    head2LocsCond1 = zeros(nSites,1); %mobile head
    pausedTimeAtEachLocCond1 = zeros(nSites,1); %keep track of how long it has been paused
    %and if chosen where is the block location
    cond1Locs = zeros(nCond1Total,6); %LocHead1 LocHead2 Orientation Pausedtime Willswitch Switchlocation
    
    head1LocsCond2 = zeros(nSites,1); %static head
    head2LocsCond2 = zeros(nSites,1); %mobile head
    pausedTimeAtEachLocCond2 = zeros(nSites,1); %keep track of how long it has been paused
    %and if chosen where is the block location
    cond2Locs = zeros(nCond2Total,6); %LocHead1 LocHead2 Orientation Pausedtime Willswitch Switchlocation
    
    
    nLoops = 0;

    for i=1:1%only put one condensin
        %Initialise Condesin locations
        loc = floor((rand*nSites)+1);
        if rand<0.5
            cond1Locs(i,3)=1;
        else
            cond1Locs(i,3)=-1;
        end
        
%         if i==1
%             loc = 5000;
%             cond2Locs(i,3)=-1;
%         elseif i==2
%             loc = 4500;
%             cond2Locs(i,3)=1;
%         end
        cond1Locs(i,1)=loc;
        cond1Locs(i,2)=loc;
        
        head1LocsCond1(loc)=head1LocsCond1(loc)+1; %Keep track of where the head1 are
        head2LocsCond1(loc)=head2LocsCond1(loc)+1; %Keep track of where the head2 are
%         
%         pausedTimeAtEachLoc(loc)=0; %Keep track of how long the Condensin has been paused
        nLoops = nLoops+1;
    end
    %disp(sum(head1Locs))
    %boundTime = 0;
    for t=1:tMaxSim
        %boundTime = boundTime+delT;
        %Change locations of head2
        %For all non-zero entries in cond2Locs change head2
        nonzeroIndicesCond1 = find(cond1Locs(:,1)~=0);
        zeroIndicesCond1 = find(cond1Locs(:,1)==0);
        
        nonzeroIndicesCond2 = find(cond2Locs(:,1)~=0);
        zeroIndicesCond2 = find(cond2Locs(:,1)==0);
        
        %Loop length change
        %calculations for condensin1
        for i=1:length(nonzeroIndicesCond1)
            currLoopIndex = nonzeroIndicesCond1(i);
            currLoopH2 = cond1Locs(currLoopIndex,2);
            currLoopH1 = cond1Locs(currLoopIndex,1);
            currLoopOri = cond1Locs(currLoopIndex,3);
            
            currLoopWillSwitch = cond1Locs(currLoopIndex,5);
            %Multiply orientation with step size length to get step-size
%             if currLoopWillSwitch==1
%                 disp(currLoopWillSwitch)
%             end
            stepSz = currLoopOri*round((SamplePofL(cdf_inv_pofl,rand,delX))/(0.34*monomerSize));%+(SamplePofL(cdf_inv_pofl,rand,delX))/(0.34*monomerSize)); % /3.4*10/monomerSize
            %stepSz = currLoopOri*round(extrusionRateSim*delT);
            
            %Check if next step is possible
            if currLoopH2+stepSz < nSites && currLoopH2+stepSz > 1
                %Find the nearest block
                %iBlock is distance of block from current head2
                if currLoopWillSwitch==0
                    iBlock=0;
                    for iterBlock = currLoopOri:currLoopOri:stepSz    %All possible stepping locations
                        iBlock = iterBlock;
                        if head1LocsCond1(currLoopH2+iterBlock) ~= 0 || head2LocsCond1(currLoopH2+iterBlock) ~= 0 || head1LocsCond2(currLoopH2+iterBlock) ~= 0 || head2LocsCond2(currLoopH2+iterBlock) ~= 0
                            break
                        end
                    end

                    numBlockades = head1LocsCond1(currLoopH2+iBlock)+head2LocsCond1(currLoopH2+iBlock)+head1LocsCond2(currLoopH2+iBlock)+head2LocsCond2(currLoopH2+iBlock); %Will be zero if no blocks
                    waitingTime = numBlockades*pauseTime;
                    if iBlock==stepSz && waitingTime==0     %No block to be considered
                        %Take step
%                         if stepSz<0
%                             disp(stepSz);
%                         end
                        head2LocsCond1(currLoopH2+stepSz) = head2LocsCond1(currLoopH2+stepSz)+1;
                        head2LocsCond1(currLoopH2) = head2LocsCond1(currLoopH2)-1;
                        cond1Locs(currLoopIndex,2) = cond1Locs(currLoopIndex,2)+stepSz;
                        cond1Locs(currLoopIndex,4) = 0; %Pause reset to zero
                        %disp(cond2Locs(2,2))
                    else        %Block to be considered
                        
                        if waitingTime==0
                            %Take step
                            head2LocsCond1(currLoopH2+stepSz) = head2LocsCond1(currLoopH2+stepSz)+1;
                            head2LocsCond1(currLoopH2) = head2LocsCond1(currLoopH2)-1;
                            cond1Locs(currLoopIndex,2) = cond1Locs(currLoopIndex,2)+stepSz;
                            cond1Locs(currLoopIndex,4) = 0; %Pause reset to zero
                        elseif abs(iBlock) ~= 1 && waitingTime~=0      %Block not adjacent move next to block
                            
                            %First check if switching should be considered
                            %If block location is a head2
                            if head2LocsCond1(currLoopH2+iBlock)~=0  %The block is a head2
                                %find corresponding head1
                                blockHead2Ind = currLoopH2+iBlock;
                                correspondingH1Locs = cond1Locs(cond1Locs(:,2)==blockHead2Ind,1);
                                %Pick the maximum loop
                                distVec = abs(correspondingH1Locs-blockHead2Ind);
                                maxDistLoop = max(distVec);
                                maxH1 = correspondingH1Locs(distVec==maxDistLoop);
                                %check if curr loop within this loop
                                if (currLoopH1-maxH1)*currLoopOri>=0 %Current loop nested
                                    if rand>probSwitch
                                        cond1Locs(currLoopIndex,5) = 1;
                                        cond1Locs(currLoopIndex,6) = maxH1;
                                    end
                                end
                                %disp(switchArray(currLoopH2,1))
                            end
                            %Check random switch
                            %if true find correspoding head2
                            %perform stop and wait subroutine
                            %else continue

                            blockedStep = iBlock-currLoopOri; %-1 for 1 and vice versa
                            head2LocsCond1(currLoopH2+blockedStep) = head2LocsCond1(currLoopH2+blockedStep)+1;
                            head2LocsCond1(currLoopH2) = head2LocsCond1(currLoopH2)-1;
                            cond1Locs(currLoopIndex,2) = cond1Locs(currLoopIndex,2)+blockedStep;
                            cond1Locs(currLoopIndex,4) = 0; %Pause reset to zero   

                        elseif abs(iBlock) == 1 && cond1Locs(currLoopIndex,4) < waitingTime %Block adjacent but pause less
                            %Wait here
                            cond1Locs(currLoopIndex,4) = cond1Locs(currLoopIndex,4) + delT; %keep track of how long waiting

                        elseif abs(iBlock) == 1 && cond1Locs(currLoopIndex,4) >= waitingTime %Block adjacent but pause exceeded
                            %Move on block
                            %disp("here")
                            %disp("here")
                            blockedStep = iBlock;
                            head2LocsCond1(currLoopH2+blockedStep) = head2LocsCond1(currLoopH2+blockedStep)+1;
                            head2LocsCond1(currLoopH2) = head2LocsCond1(currLoopH2)-1;
                            cond1Locs(currLoopIndex,2) = cond1Locs(currLoopIndex,2)+blockedStep;
                            cond1Locs(currLoopIndex,4) = 0; %Pause reset to zero
    %                         if currLoopIndex==9
    %                             disp(t)
    %                             disp([currLoopIndex,currLoopH2,currLoopOri,iBlock,stepSz])
    %                         end
                        else
                            disp("here")
                            input('')
                        end
                    end
                elseif currLoopWillSwitch==1 %%Is setup for switching, make calculations accordingly
                    %Compute blocking time and such
                    disp("switching direction")
                    iBlock = cond1Locs(currLoopIndex,6)-currLoopH2;
                    numBlockades = head1LocsCond1(currLoopH2+iBlock)+head2LocsCond1(currLoopH2+iBlock); %Will be zero if no blocks
                    waitingTime = numBlockades*pauseTime;
                    %if pause not exceeded, wait
                    if cond1Locs(currLoopIndex,4) < waitingTime %Block adjacent but pause less
                        %Wait here
                        cond1Locs(currLoopIndex,4) = cond1Locs(currLoopIndex,4) + delT;
                    elseif cond1Locs(currLoopIndex,4) >= waitingTime
                        %if pause exceeded, pass through block, change
                        %orientation, change h2 location, change switchArray
                        %values
                        blockedStep = iBlock;
                        head2LocsCond1(currLoopH2+blockedStep) = head2LocsCond1(currLoopH2+blockedStep)+1;
                        head2LocsCond1(currLoopH2) = head2LocsCond1(currLoopH2)-1;
                        cond1Locs(currLoopIndex,2) = cond1Locs(currLoopIndex,2)+blockedStep;
                        cond1Locs(currLoopIndex,3) = -cond1Locs(currLoopIndex,3); %Reverse orientation
                        cond1Locs(currLoopIndex,4) = 0; %Pause reset to zero
                        cond1Locs(currLoopIndex,5) = 0;
                        cond1Locs(currLoopIndex,6) = 0;
                    end
                end
            end
        end
        
        %Make calculations for condensin 2
        for i=1:length(nonzeroIndicesCond2)
            currLoopIndex = nonzeroIndicesCond2(i);
            currLoopH2 = cond2Locs(currLoopIndex,2);
            currLoopH1 = cond2Locs(currLoopIndex,1);
            currLoopOri = cond2Locs(currLoopIndex,3);
            
            currLoopWillSwitch = cond2Locs(currLoopIndex,5);
            %Multiply orientation with step size length to get step-size
%             if currLoopWillSwitch==1
%                 disp(currLoopWillSwitch)
%             end
            stepSz = currLoopOri*round((SamplePofL(cdf_inv_pofl,rand,delX))/(0.34*monomerSize));%+(SamplePofL(cdf_inv_pofl,rand,delX))/(0.34*monomerSize)); % /3.4*10/monomerSize
            %stepSz = currLoopOri*round(extrusionRateSim*delT);
            
            %Check if next step is possible
            if currLoopH2+stepSz < nSites && currLoopH2+stepSz > 1
                %Find the nearest block
                %iBlock is distance of block from current head2
                if currLoopWillSwitch==0
                    iBlock=0;
                    for iterBlock = currLoopOri:currLoopOri:stepSz    %All possible stepping locations
                        iBlock = iterBlock;
                        if head1LocsCond1(currLoopH2+iterBlock) ~= 0 || head2LocsCond1(currLoopH2+iterBlock) ~= 0 || head1LocsCond2(currLoopH2+iterBlock) ~= 0 || head2LocsCond2(currLoopH2+iterBlock) ~= 0
                            break
                        end
                    end

                    numBlockades = head1LocsCond1(currLoopH2+iBlock)+head2LocsCond1(currLoopH2+iBlock)+head1LocsCond2(currLoopH2+iBlock)+head2LocsCond2(currLoopH2+iBlock); %Will be zero if no blocks
                    waitingTime = numBlockades*pauseTime;
                    if iBlock==stepSz && waitingTime==0     %No block to be considered
                        %Take step
%                         if stepSz<0
%                             disp(stepSz);
%                         end
                        head2LocsCond2(currLoopH2+stepSz) = head2LocsCond2(currLoopH2+stepSz)+1;
                        head2LocsCond2(currLoopH2) = head2LocsCond2(currLoopH2)-1;
                        cond2Locs(currLoopIndex,2) = cond2Locs(currLoopIndex,2)+stepSz;
                        cond2Locs(currLoopIndex,4) = 0; %Pause reset to zero
                        %disp(cond2Locs(2,2))
                    else        %Block to be considered
                        
                        if waitingTime==0
                            %Take step
                            head2LocsCond2(currLoopH2+stepSz) = head2LocsCond2(currLoopH2+stepSz)+1;
                            head2LocsCond2(currLoopH2) = head2LocsCond2(currLoopH2)-1;
                            cond2Locs(currLoopIndex,2) = cond2Locs(currLoopIndex,2)+stepSz;
                            cond2Locs(currLoopIndex,4) = 0; %Pause reset to zero
                        elseif abs(iBlock) ~= 1 && waitingTime~=0      %Block not adjacent move next to block
                            
                            %First check if switching should be considered
                            %If block location is a head2
                            if head2LocsCond2(currLoopH2+iBlock)~=0  %The block is a head2
                                %find corresponding head1
                                blockHead2Ind = currLoopH2+iBlock;
                                correspondingH1Locs = cond2Locs(cond2Locs(:,2)==blockHead2Ind,1);
                                %Pick the maximum loop
                                distVec = abs(correspondingH1Locs-blockHead2Ind);
                                maxDistLoop = max(distVec);
                                maxH1 = correspondingH1Locs(distVec==maxDistLoop);
                                %check if curr loop within this loop
                                if (currLoopH1-maxH1)*currLoopOri>=0 %Current loop nested
                                    if rand>probSwitch
                                        cond2Locs(currLoopIndex,5) = 1;
                                        cond2Locs(currLoopIndex,6) = maxH1;
                                    end
                                end
                                %disp(switchArray(currLoopH2,1))
                            end
                            %Check random switch
                            %if true find correspoding head2
                            %perform stop and wait subroutine
                            %else continue

                            blockedStep = iBlock-currLoopOri; %-1 for 1 and vice versa
                            head2LocsCond2(currLoopH2+blockedStep) = head2LocsCond2(currLoopH2+blockedStep)+1;
                            head2LocsCond2(currLoopH2) = head2LocsCond2(currLoopH2)-1;
                            cond2Locs(currLoopIndex,2) = cond2Locs(currLoopIndex,2)+blockedStep;
                            cond2Locs(currLoopIndex,4) = 0; %Pause reset to zero   

                        elseif abs(iBlock) == 1 && cond2Locs(currLoopIndex,4) < waitingTime %Block adjacent but pause less
                            %Wait here
                            cond2Locs(currLoopIndex,4) = cond2Locs(currLoopIndex,4) + delT; %keep track of how long waiting

                        elseif abs(iBlock) == 1 && cond2Locs(currLoopIndex,4) >= waitingTime %Block adjacent but pause exceeded
                            %Move on block
                            %disp("here")
                            %disp("here")
                            blockedStep = iBlock;
                            head2LocsCond2(currLoopH2+blockedStep) = head2LocsCond2(currLoopH2+blockedStep)+1;
                            head2LocsCond2(currLoopH2) = head2LocsCond2(currLoopH2)-1;
                            cond2Locs(currLoopIndex,2) = cond2Locs(currLoopIndex,2)+blockedStep;
                            cond2Locs(currLoopIndex,4) = 0; %Pause reset to zero
    %                         if currLoopIndex==9
    %                             disp(t)
    %                             disp([currLoopIndex,currLoopH2,currLoopOri,iBlock,stepSz])
    %                             input('error')
    %                         end
                        else
                            disp("here")
                            input('')
                        end
                    end
                elseif currLoopWillSwitch==1 %%Is setup for switching, make calculations accordingly
                    %Compute blocking time and such
                    disp("switching direction")
                    iBlock = cond2Locs(currLoopIndex,6)-currLoopH2;
                    numBlockades = head1LocsCond2(currLoopH2+iBlock)+head2LocsCond2(currLoopH2+iBlock); %Will be zero if no blocks
                    waitingTime = numBlockades*pauseTime;
                    %if pause not exceeded, wait
                    if cond2Locs(currLoopIndex,4) < waitingTime %Block adjacent but pause less
                        %Wait here
                        cond2Locs(currLoopIndex,4) = cond2Locs(currLoopIndex,4) + delT;
                    elseif cond2Locs(currLoopIndex,4) >= waitingTime
                        %if pause exceeded, pass through block, change
                        %orientation, change h2 location, change switchArray
                        %values
                        blockedStep = iBlock;
                        head2LocsCond2(currLoopH2+blockedStep) = head2LocsCond2(currLoopH2+blockedStep)+1;
                        head2LocsCond2(currLoopH2) = head2LocsCond2(currLoopH2)-1;
                        cond2Locs(currLoopIndex,2) = cond2Locs(currLoopIndex,2)+blockedStep;
                        cond2Locs(currLoopIndex,3) = -cond2Locs(currLoopIndex,3); %Reverse orientation
                        cond2Locs(currLoopIndex,4) = 0; %Pause reset to zero
                        cond2Locs(currLoopIndex,5) = 0;
                        cond2Locs(currLoopIndex,6) = 0;
                    end
                end
            end
        end
        
        %input('')
        %Remove condensin
        %For all non-zero entries in cond2Locs remove randomly and make zero
        for i=1:length(nonzeroIndicesCond1)
            if rand<probOffCond1
                %fprintf(flBound,"%f\n",boundTime);
                %boundTime = 0;
               
                currLoopIndex = nonzeroIndicesCond1(i);
                currLoopH1 = cond1Locs(currLoopIndex,1);
                currLoopH2 = cond1Locs(currLoopIndex,2);
                head1LocsCond1(currLoopH1) = head1LocsCond1(currLoopH1)-1;
                head2LocsCond1(currLoopH2) = head2LocsCond1(currLoopH2)-1;
                cond1Locs(currLoopIndex,:) = 0;
            end
        end
        
        for i=1:length(nonzeroIndicesCond2)
            if rand<probOffCond2
                %fprintf(flBound,"%f\n",boundTime);
                %boundTime = 0;
               
                currLoopIndex = nonzeroIndicesCond2(i);
                currLoopH1 = cond2Locs(currLoopIndex,1);
                currLoopH2 = cond2Locs(currLoopIndex,2);
                head1LocsCond2(currLoopH1) = head1LocsCond2(currLoopH1)-1;
                head2LocsCond2(currLoopH2) = head2LocsCond2(currLoopH2)-1;
                cond2Locs(currLoopIndex,:) = 0;
            end
        end
        
        %Add condensin
        %For all zero entries in cond2Locs randomly put in some location
        for i=1:length(zeroIndicesCond1)
            if rand<probOnCond1
                loc = floor((rand*nSites)+1);
                currLoopIndex = zeroIndicesCond1(i);
                cond1Locs(currLoopIndex,1) = loc;
                cond1Locs(currLoopIndex,2) = loc;
                head1LocsCond1(loc) = head1LocsCond1(loc)+1;
                head2LocsCond1(loc) = head2LocsCond1(loc)+1;
                if rand<0.5
                    cond1Locs(currLoopIndex,3)=1;
                else
                    cond1Locs(currLoopIndex,3)=-1;
                end
            end
        end
        
        for i=1:length(zeroIndicesCond2)
            if rand<probOnCond2
                loc = floor((rand*nSites)+1);
                currLoopIndex = zeroIndicesCond2(i);
                cond2Locs(currLoopIndex,1) = loc;
                cond2Locs(currLoopIndex,2) = loc;
                head1LocsCond2(loc) = head1LocsCond2(loc)+1;
                head2LocsCond2(loc) = head2LocsCond2(loc)+1;
                if rand<0.5
                    cond2Locs(currLoopIndex,3)=1;
                else
                    cond2Locs(currLoopIndex,3)=-1;
                end
            end
        end
        %Print the loopList
        if mod(t,freq)==0
            save('LoopListsCond1/LoopList_'+string(run)+'.txt', 'cond1Locs', '-ASCII','-append');
            save('LoopListsCond2/LoopList_'+string(run)+'.txt', 'cond2Locs', '-ASCII','-append');
            
            dum = dum + (mean(abs(nonzeros(cond2Locs(:,1))-nonzeros(cond2Locs(:,2)))));
            ndum = ndum+1;
        end
        %Save avg loopLen
        avgLoopLen(t,2) = avgLoopLen(t,2) + mean(abs(nonzeros(cond1Locs(:,1))-nonzeros(cond1Locs(:,2))));
        avgLoopLen(t,3) = avgLoopLen(t,3) + std(abs(nonzeros(cond2Locs(:,1))-nonzeros(cond2Locs(:,2))))^2;
        avgLoopLen(t,1) = t*delT;
    end
    disp([dum/ndum,length(nonzeros(cond2Locs(:,1)))])
end
avgLoopLen(:,2) = avgLoopLen(:,2)./runMax;
avgLoopLen(:,3) = avgLoopLen(:,3)./runMax;
avgLoopLen(:,3) = avgLoopLen(:,3).^0.5;
% hold on m
box on
plot(avgLoopLen(:,1),avgLoopLen(:,2)*monomerSize/10^4,'LineWidth',4)
set(gca, 'fontsize', 20)
title('Condensin I P(s) - 60 min','FontSize',20)
xlabel('time (seconds)','FontSize',20)
ylabel('Average Loop Size (10 Kbp)','FontSize',20)
% 
%save('AvgLoopSizeForPlotting/AvgLoopSizeCond2.txt', 'avgLoopLen', '-ASCII');
