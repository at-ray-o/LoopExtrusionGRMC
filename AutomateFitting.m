%% Code for computing P(s) at different times

% Load the cdf_inv of actual distribution

cdf_inv_pofl = load("../PofLSim/PofLNumericalGenerate/pofL_cdf_inv.mat");
cdf_inv_pofl = cdf_inv_pofl.dummyVar;
delX = cdf_inv_pofl(2,1)-cdf_inv_pofl(1,1);
%%


%Intialise Condensin at different sites
%rng(445)
monomerSize = 100; %bps
delT = 1; %seconds
pauseTime = 7; %seconds
residenceTimeReal = 6*60; %seconds
nSites = 1000000; 
probSwitch = 1;
%Condensin II
nCond2Total = 343;
nCond2Bound = 210;
stepProb = 0.98;
extrusionRateReal = 1500; %bp/s
extrusionRateSim = extrusionRateReal/monomerSize;
residenceTimeSim = residenceTimeReal/delT;
%stepSzLenOnce=extrusionRateSim;

nCond2Unbound = nCond2Total-nCond2Bound;

probCond2Ini = nCond2Bound/nSites;
probOff = 1/(residenceTimeSim);
probOn = nCond2Bound*probOff/nCond2Unbound;

%Simulation parameters
freq = 360;
tMaxSim = 60*60/delT;
avgLoopLen = zeros(tMaxSim,3);
runMax = 5;
dum = 0;
ndum = 0;
%flBound = fopen('boundTime.txt','a');
for run = 1:runMax
    head1Locs = zeros(nSites,1); %static head
    head2Locs = zeros(nSites,1); %mobile head
    pausedTimeAtEachLoc = zeros(nSites,1); %keep track of how long it has been paused
    %and if chosen where is the block location
    cond2Locs = zeros(nCond2Total,6); %LocHead1 LocHead2 Orientation Pausedtime Willswitch Switchlocation
    nLoops = 0;

    for i=1:1%only put one condensin
%         if i==1
%             cond2Locs(i,1)=5700;
%             cond2Locs(i,2)=4500;
%             cond2Locs(i,3)=-1;
%             head1Locs(5700)=head1Locs(5700)+1; 
%             head2Locs(4500)=head2Locs(4500)+1;
%         elseif i==2
%             cond2Locs(i,1)=5250;
%             cond2Locs(i,2)=5050;
%             cond2Locs(i,3)=-1;
%             head1Locs(5250)=head1Locs(5250)+1; 
%             head2Locs(5050)=head2Locs(5050)+1;
%         end
        %Initialise Condesin locations
        loc = floor((rand*nSites)+1);
        if rand<0.5
            cond2Locs(i,3)=1;
        else
            cond2Locs(i,3)=-1;
        end
        
%         if i==1
%             loc = 5000;
%             cond2Locs(i,3)=-1;
%         elseif i==2
%             loc = 4500;
%             cond2Locs(i,3)=1;
%         end
        cond2Locs(i,1)=loc;
        cond2Locs(i,2)=loc;
        
        head1Locs(loc)=head1Locs(loc)+1; %Keep track of where the head1 are
        head2Locs(loc)=head2Locs(loc)+1; %Keep track of where the head2 are
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
        nonzeroIndices = find(cond2Locs(:,1)~=0);
        zeroIndices = find(cond2Locs(:,1)==0);
        
        %Loop length change
        
        for i=1:length(nonzeroIndices)
            currLoopIndex = nonzeroIndices(i);
            currLoopH2 = cond2Locs(currLoopIndex,2);
            currLoopH1 = cond2Locs(currLoopIndex,1);
            currLoopOri = cond2Locs(currLoopIndex,3);
            
            currLoopWillSwitch = cond2Locs(currLoopIndex,5);
            %Multiply orientation with step size length to get step-size
%             if currLoopWillSwitch==1
%                 disp(currLoopWillSwitch)
%             end
            if rand<=stepProb
                willStep = 1;
            else
                willStep = 0;
            end
            stepSz = willStep*currLoopOri*round((SamplePofL(cdf_inv_pofl,rand,delX))/(0.34*monomerSize));
            %stepSz = currLoopOri*round(extrusionRateSim*delT);
            
            %Check if next step is possible
            if currLoopH2+stepSz < nSites && currLoopH2+stepSz > 1
                %Find the nearest block
                %iBlock is distance of block from current head2
                if currLoopWillSwitch==0
                    iBlock=0;
                    for iterBlock = currLoopOri:currLoopOri:stepSz    %All possible stepping locations
                        iBlock = iterBlock;
                        if head1Locs(currLoopH2+iterBlock) ~= 0 || head2Locs(currLoopH2+iterBlock) ~= 0
                            break
                        end
                    end

                    numBlockades = head1Locs(currLoopH2+iBlock)+head2Locs(currLoopH2+iBlock); %Will be zero if no blocks
                    waitingTime = numBlockades*pauseTime;
                    if iBlock==stepSz && waitingTime==0     %No block to be considered
                        %Take step
%                         if stepSz<0
%                             disp(stepSz);
%                         end
                        head2Locs(currLoopH2+stepSz) = head2Locs(currLoopH2+stepSz)+1;
                        head2Locs(currLoopH2) = head2Locs(currLoopH2)-1;
                        cond2Locs(currLoopIndex,2) = cond2Locs(currLoopIndex,2)+stepSz;
                        cond2Locs(currLoopIndex,4) = 0; %Pause reset to zero
                        %disp(cond2Locs(2,2))
                    else        %Block to be considered
                        
                        if waitingTime==0
                            %Take step
                            head2Locs(currLoopH2+stepSz) = head2Locs(currLoopH2+stepSz)+1;
                            head2Locs(currLoopH2) = head2Locs(currLoopH2)-1;
                            cond2Locs(currLoopIndex,2) = cond2Locs(currLoopIndex,2)+stepSz;
                            cond2Locs(currLoopIndex,4) = 0; %Pause reset to zero
                        elseif abs(iBlock) ~= 1 && waitingTime~=0      %Block not adjacent move next to block
                            
                            %First check if switching should be considered
                            %If block location is a head2
                            if head2Locs(currLoopH2+iBlock)~=0  %The block is a head2
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
                            head2Locs(currLoopH2+blockedStep) = head2Locs(currLoopH2+blockedStep)+1;
                            head2Locs(currLoopH2) = head2Locs(currLoopH2)-1;
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
                            head2Locs(currLoopH2+blockedStep) = head2Locs(currLoopH2+blockedStep)+1;
                            head2Locs(currLoopH2) = head2Locs(currLoopH2)-1;
                            cond2Locs(currLoopIndex,2) = cond2Locs(currLoopIndex,2)+blockedStep;
                            cond2Locs(currLoopIndex,4) = 0; %Pause reset to zero
    %                         if currLoopIndex==9
    %                             disp(t)
    %                             disp([currLoopIndex,currLoopH2,currLoopOri,iBlock,stepSz])
    %                             input('poop')
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
                    numBlockades = head1Locs(currLoopH2+iBlock)+head2Locs(currLoopH2+iBlock); %Will be zero if no blocks
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
                        head2Locs(currLoopH2+blockedStep) = head2Locs(currLoopH2+blockedStep)+1;
                        head2Locs(currLoopH2) = head2Locs(currLoopH2)-1;
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
        for i=1:length(nonzeroIndices)
            if rand<probOff
                %fprintf(flBound,"%f\n",boundTime);
                %boundTime = 0;
               
                currLoopIndex = nonzeroIndices(i);
                currLoopH1 = cond2Locs(currLoopIndex,1);
                currLoopH2 = cond2Locs(currLoopIndex,2);
                head1Locs(currLoopH1) = head1Locs(currLoopH1)-1;
                head2Locs(currLoopH2) = head2Locs(currLoopH2)-1;
                cond2Locs(currLoopIndex,:) = 0;
            end
        end
        %Add condensin
        %For all zero entries in cond2Locs randomly put in some location
        for i=1:length(zeroIndices)
            if rand<probOn
                loc = floor((rand*nSites)+1);
                currLoopIndex = zeroIndices(i);
                cond2Locs(currLoopIndex,1) = loc;
                cond2Locs(currLoopIndex,2) = loc;
                head1Locs(loc) = head1Locs(loc)+1;
                head2Locs(loc) = head2Locs(loc)+1;
                if rand<0.5
                    cond2Locs(currLoopIndex,3)=1;
                else
                    cond2Locs(currLoopIndex,3)=-1;
                end
            end
        end
        %Print the loopList
        if mod(t,freq)==0
            save('LoopListsCond1/LoopList_'+string(run)+'.txt', 'cond2Locs', '-ASCII','-append');
            
            dum = dum + (mean(abs(nonzeros(cond2Locs(:,1))-nonzeros(cond2Locs(:,2)))));
            ndum = ndum+1;
        end
        %Save avg loopLen
        avgLoopLen(t,2) = avgLoopLen(t,2) + mean(abs(nonzeros(cond2Locs(:,1))-nonzeros(cond2Locs(:,2))));
        avgLoopLen(t,3) = avgLoopLen(t,3) + std(abs(nonzeros(cond2Locs(:,1))-nonzeros(cond2Locs(:,2))))^2;
        avgLoopLen(t,1) = t*delT;
    end
    disp([dum/ndum,length(nonzeros(cond2Locs(:,1)))])
end
avgLoopLen(:,2) = avgLoopLen(:,2)./runMax;
avgLoopLen(:,3) = avgLoopLen(:,3)./runMax;
avgLoopLen(:,3) = avgLoopLen(:,3).^0.5;
% hold on m
% box on
% plot(avgLoopLen(:,1),avgLoopLen(:,2)*monomerSize/10^4,'LineWidth',4)
% set(gca, 'fontsize', 20)
% title('Condensin I P(s) - 60 min','FontSize',20)
% xlabel('time (seconds)','FontSize',20)
% ylabel('Average Loop Size (10 Kbp)','FontSize',20)
% 
%save('AvgLoopSizeForPlotting/AvgLoopSizeCond2.txt', 'avgLoopLen', '-ASCII');

%% Compute P(s)

%Code to compute average end to end distance and it's variance
inp = [7,1];
%expPoS = csead("CapH2maid30minUpto1000.csv");
len = 10000; %1 bead = 1 kbp;
b = 1;
rc = 1;
kBT=1;
sqrtkBT = sqrt(kBT);
loopSepParam = inp(1);
relLoopStrength = inp(2);
k = 3/b^2*kBT;
contacts = zeros(len,len);
nSamples = 1;

nTimes = 1;
nLoops = 343;
nTimes = 9;
runMax = 5;
pArrAllTimes = zeros(nTimes,len-1);
contacts = zeros(len,len);
for t=9:nTimes
    %disp(t);
    sigma = zeros(len,len);
    
    for run=1:runMax
        disp(run);
        mat = zeros(len,len);
        crossLink = zeros(len,len);
        bias = zeros(len,1);
        isLoop = zeros(1,len);
        for i = 1:len-1
            mat(i,i+1) = -k;
            mat(i+1,i) = -k;
            mat(i,i) = 2*k;
        end
        mat(len,len)=k;
        %relLoopStrength = loopStrength;
        
        %load loop list
        loopList = load("LoopListsCond1/LoopList_"+string(run)+".txt","-ascii");
        
        %Save it once
%         save("bondlist.dat","loopList","-ascii");
%         input("poop");
        currLoopList = loopList(1+t*nLoops:(t+1)*nLoops,1:2);
        kLoop = relLoopStrength*k;
        nonZeroIndex = find(currLoopList(:,1)~=0);
        for i=1:length(nonZeroIndex)
            startLen = round(currLoopList(nonZeroIndex(i),1)/100);
            stopLen = round(currLoopList(nonZeroIndex(i),2)/100);
            if startLen==0
                startLen=1;
            end
            if stopLen==0
                stopLen=1;
            end
            %disp([startLen,stopLen])
            isLoop(startLen)=1;
            isLoop(stopLen)=1;
            mat(startLen,stopLen) = mat(startLen,stopLen) - kLoop;
            mat(stopLen,startLen) = mat(stopLen,startLen) - kLoop;
            mat(startLen,startLen) = mat(startLen,startLen)+kLoop;
            mat(stopLen,stopLen) = mat(stopLen,stopLen)+kLoop;
            
        end

        matInv = inv(mat);
        for i=1:len-1
            for j=i+1:len
                %Compute the sigma value for i and j
                si_2 = matInv(i,i);
                sj_2 = matInv(j,j);
                sij = matInv(i,j);
                sigma(i,j) = sqrt(si_2+sj_2-2*sij);
                sigma(j,i) = sigma(i,j);
                val = erf(rc/(sqrt(2)*sigma(i,j)))-sqrt(2/pi)*exp(-rc^2/(2*sigma(i,j)^2))*rc/sigma(i,j);
                contacts(i,j) = contacts(i,j)+val;
                contacts(j,i) = contacts(i,j);
            end
        end
    
    end
    %We only want the final term of the inverse
    contacts = contacts/runMax;
    
    sArr = 1:1:len-1;
    pArr = zeros(1,len-1);
    counts = zeros(1,len-1);
    for i = 1:len-1
        %Compute avg prob from i to i+s
        for j = i+1:len
            s = abs(j-i);
            counts(s) = counts(s)+1;
            pArr(s) = pArr(s)+contacts(i,j);
        end
    end
    calib = 10;
    pArr = pArr./counts;
    pArrSkip6 = pArr(2:len-1);
    pArrSkip6Min10 = pArr(1:len-1);
    mag = pArr(calib);
    %mag = 1;
    pArrSkip6 = pArrSkip6 ./mag;
    pArrSkip6Min10 = pArrSkip6Min10 ./mag;
    pArrAllTimes(t+1,:) = pArrSkip6Min10;
end

hold on
box on
set(gca, 'xscale','log')
set(gca, 'yscale','log')
set(gca, 'fontsize', 20)
set(gca, 'linewidth', 4)

calib = 10;
mag = pArr(calib);
pArr = pArr ./ mag;
plot((sArr),(pArr),'LineWidth',4);

mitoDat = csvread("ExpCSVs/Condensin2_30min.csv");
mitoDat(:,1) = mitoDat(:,1); %Setting step-size to 10kbp
interpDat = interp1(mitoDat(:,1),mitoDat(:,2),sArr,'linear');
mag = interpDat(calib);
interpDatSkip6 = interpDat(2:len-1) ./ mag;
interpDatSkip6Min10 = interpDat(1:len-1) ./ mag;
interpDatSkip6Min10(1) = mitoDat(1,2)/mag;
plot((sArr),(interpDatSkip6Min10),'LineWidth',4);
diff = norm(log(interpDatSkip6Min10)-log(pArr))


%% Delete the loop list files
for run=1:20
    delete("LoopListsCond1/LoopList_"+string(run)+".txt");
end