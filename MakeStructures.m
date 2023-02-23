%% Code to compute structures, given loop configurations.
len = 10000; %1 bead = 1 kbp;
b = 1;
rc = 1;
kBT=1;
sqrtkBT = sqrt(kBT);
relLoopStrength = 1;
k = 3/b^2*kBT;
nSamples = 100;

nLoopsCond1 = 2550;
nLoopsCond2 = 343;
nTimes = 9;
runMax = 1;
fileID = fopen('structureTest.xyz','w');

meanEndToendDist = zeros(1,runMax);
varEndToendDist = zeros(1,runMax);

isLoop = zeros(len,1);

for run=1:runMax
    disp(run);
    for t=9:nTimes
        disp(t);
        mat = zeros(len,len);
        crossLink = zeros(len,len);
        
        isLoop = zeros(1,len);
        for i = 1:len-1
            mat(i,i+1) = -k;
            mat(i+1,i) = -k;
            mat(i,i) = 2*k;
        end
        mat(len,len)=k;
        kLoop = relLoopStrength*k;
        start = 0;
        stop = 0;
        loopListCond1 = load("LoopListsCond1/LoopList_"+string(run)+".txt","-ascii");
        currLoopList = loopListCond1(1+t*nLoopsCond1:(t+1)*nLoopsCond1,1:2);
        nonZeroIndex = find(currLoopList(:,1)~=0);
        for i=1:length(nonZeroIndex)
            start = round(currLoopList(nonZeroIndex(i),1)/100);
            stop = round(currLoopList(nonZeroIndex(i),2)/100);
            if start==0
                start=1;
            end
            if stop==0
                stop=1;
            end
            isLoop(start)=1;
            isLoop(stop)=1;
            mat(start,stop) = mat(start,stop) - kLoop;
            mat(stop,start) = mat(stop,start) - kLoop;
            mat(start,start) = mat(start,start)+kLoop;
            mat(stop,stop) = mat(stop,stop)+kLoop;

        end
        
        loopListCond2 = load("LoopListsCond2/LoopList_"+string(run)+".txt","-ascii");
        currLoopList = loopListCond2(1+t*nLoopsCond2:(t+1)*nLoopsCond2,1:2);
        nonZeroIndex = find(currLoopList(:,1)~=0);
        for i=1:length(nonZeroIndex)
            start = round(currLoopList(nonZeroIndex(i),1)/100);
            stop = round(currLoopList(nonZeroIndex(i),2)/100);
            if start==0
                start=1;
            end
            if stop==0
                stop=1;
            end
            isLoop(start)=1;
            isLoop(stop)=1;
            mat(start,stop) = mat(start,stop) - kLoop;
            mat(stop,start) = mat(stop,start) - kLoop;
            mat(start,start) = mat(start,start)+kLoop;
            mat(stop,stop) = mat(stop,stop)+kLoop;

        end
        
        sigma = inv(mat);
        cholSigma = chol(sigma);
        transform = cholSigma';
        xCoord = zeros(len,1);
        yCoord = zeros(len,1);
        zCoord = zeros(len,1);
        lambdaSamples = zeros(len,1);
        for i=1:nSamples
            for j=1:len-1
                lambdaSamples(j) = normrnd(0,1);
            end
            xCoord = transform*lambdaSamples;
            for j=1:len-1
                lambdaSamples(j) = normrnd(0,1);
            end
            yCoord = transform*lambdaSamples;
            for j=1:len-1
                lambdaSamples(j) = normrnd(0,1);
            end
            zCoord = transform*lambdaSamples;        
            xCoord = xCoord - mean(xCoord);
            yCoord = yCoord - mean(yCoord);
            zCoord = zCoord - mean(zCoord);            
      
            %Print them
            fprintf(fileID,'%d\n\n',len);
            for j=1:len

                if isLoop(j)==1
                    fprintf(fileID,'O\t%12f\t%12f\t%12f\n',xCoord(j),yCoord(j),zCoord(j));
                else
                    fprintf(fileID,'C\t%12f\t%12f\t%12f\n',xCoord(j),yCoord(j),zCoord(j));
                end
            end
            
        end
    end
end
% c = 1:len;
% scatter3(xCoord,yCoord,zCoord,2,c);
