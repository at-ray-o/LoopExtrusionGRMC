function lVal = SamplePofL(cdf_inv_pofl,unifVal,delX)
    %Find the corresponding bin to given unifVal
    if unifVal~=0 && unifVal~=1
        residue=unifVal/delX;
        lbInd = floor(residue)+1;
        ubInd = lbInd+1;
        xlb = cdf_inv_pofl(lbInd,1);
        ylb = cdf_inv_pofl(lbInd,2);
        xub = cdf_inv_pofl(ubInd,1);
        yub = cdf_inv_pofl(ubInd,2);
        lVal = ylb+(unifVal-xlb)*(yub-ylb)/(xub-xlb);
    elseif unifVal==0
        lVal = cdf_inv_pofl(1,2);
    elseif unifVal==1
        lVal = cdf_inv_pofl(end,2);   
    end
end
