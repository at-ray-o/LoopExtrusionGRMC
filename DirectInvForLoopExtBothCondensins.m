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

nLoopsCond1 = 2550;
nLoopsCond2 = 343;
nTimes = 9;
runMax = 5;
pArrAllTimes = zeros(nTimes,len-1);
contacts = zeros(len,len);
for t=9:nTimes
    disp(t);
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
        loopListCond1 = load("LoopListsCond1/LoopList_"+string(run)+".txt","-ascii");
        loopListCond2 = load("LoopListsCond2/LoopList_"+string(run)+".txt","-ascii");
        %Save it once
%         save("bondlist.dat","loopList","-ascii");
%         input("poop");
        currLoopList = loopListCond1(1+t*nLoopsCond1:(t+1)*nLoopsCond1,1:2);
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
        
        currLoopList = loopListCond2(1+t*nLoopsCond2:(t+1)*nLoopsCond2,1:2);
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
        a = input("prompt")
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

%offset = mitoDat(1,1);
%mitoDat(:,1) = mitoDat(:,1)-offset;
% sArrSkip6 = sArr(2:len-1); %To prevent boundary effects
% sArrSkip6Min10 = sArr(1:len-1); %To prevent boundary effects

% logSarr = log(sArrSkip6Min10);
% xVals = linspace(min(logSarr),max(logSarr),100);
% yValsProb = interp1(logSarr,log(pArrSkip6Min10),xVals,'linear');
% yValsInterp = interp1(logSarr,log(interpDatSkip6Min10),xVals,'linear');
% diff = norm(yValsInterp-yValsProb);
% disp(inp)
% hold on
% %set(gca,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
% set(gca, 'xscale','log')
%  set(gca, 'yscale','log')
% 
% plot((sArrSkip6Min10),(pArrSkip6Min10),'LineWidth',4,'Color',[1,0.1,0.1]);

% hold off
% title('Loop length changes P(s) curve','FontSize',20)
% xlabel('s (in 10kb)','FontSize',14)
% ylabel('P(s)','FontSize',14)
% box off;
% 
% xlims = get(gca,'XLim');
% 
% ylims = get(gca,'YLim');
% line(xlims,[ylims(2) ylims(2)],'Color','black','LineWidth',7);
% legend('A','B','C','D','E')
%csvwrite('DvaryGamma20.csv',[(sArrSkip6Min10);(pArrSkip6Min10)]');


%%

hold on
box on
set(gca, 'xscale','log')
set(gca, 'yscale','log')
set(gca, 'fontsize', 14)
for t=9:nTimes
    plot((sArr),(pArrAllTimes(t+1,:)),'LineWidth',4)
end
mitoDat = csvread("60min.csv");
mitoDat(:,1) = mitoDat(:,1); %Setting step-size to 10kbp

% mitoDat = csvread("Cond2_60min.csv");
% mitoDat(:,1) = mitoDat(:,1); %Setting step-size to 10kbp


interpDat = interp1(mitoDat(:,1),mitoDat(:,2),sArr,'linear');
mag = interpDat(calib);
interpDatSkip6 = interpDat(2:len-1) ./ mag;
interpDatSkip6Min10 = interpDat(1:len-1) ./ mag;
interpDatSkip6Min10(1) = mitoDat(1,2)/mag;

plot((sArr),(interpDatSkip6Min10),'LineStyle',':','LineWidth',4,'Color',[0,0,0]);

title('Condensin I P(s) - 60 min','FontSize',20)
xlabel('s (in 10kb)','FontSize',14)
ylabel('P(s)','FontSize',14)

data = [sArr',pArrAllTimes(t+1,:)'];

save('P(s)_both.txt', 'data', '-ASCII');