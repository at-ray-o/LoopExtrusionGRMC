%Given loop configuration draw the full chromosome configuration

%Load configuration file
run = 1;
nLoops = 2;
lenChromosome = 1000;
loopList = load("LoopListsCond1/LoopList_"+string(run)+".txt","-ascii");
nTimes = length(loopList)/nLoops;
t = 1;
printTime = 40*60; %40 minutes
angVals = 0:0.05:pi;
Frames(nTimes) = struct('cdata',[],'colormap',[]);
L = 1;
for t=1:nTimes
    currLoopList = loopList(1+t*nLoops:(t+1)*nLoops,1:2);
    nonZeroIndex = find(currLoopList(:,1)~=0);
    
    for i=1:length(nonZeroIndex)
        disp(length(nonZeroIndex))
        startLoc = round(currLoopList(nonZeroIndex(i),1)/10);
        stopLoc = round(currLoopList(nonZeroIndex(i),2)/10);
        hold on
        line(1:lenChromosome,zeros(1,lenChromosome),'LineWidth',10)
        ylim([0,10])
        %xlim([400,600])
        %Draw rectangles
        rectangle('Position',[startLoc-0.05 0 10 0.5],'FaceColor','y')
        rectangle('Position',[stopLoc-0.05 0 10 0.5],'FaceColor','r','Curvature',[1 1])
        %Draw circles between them
        R = abs(startLoc-stopLoc)/2;  %or whatever radius you want
        x = R*cos(angVals) + (startLoc+stopLoc)/2;
        y = L*sin(angVals);
        plot(x,y,'g')
        drawnow
    end
    txt = {'Time:',string(t)};
    text(500,5,txt)
    Frames(t) = getframe;
    hold off
    clf('reset')
end
movie(Frames)