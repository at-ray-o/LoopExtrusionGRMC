#include<iostream>
#include<fstream>
#include<string>
#include<time.h>
#include<math.h>


using namespace std;

int main()
{
	ifstream fTraj;
	fTraj.open("Abramo.xyz");
	long nFrames = 1000;
	int startFrame = 0;
	std::string cdum;
	double *rx, *ry, *rz;
    int ndum, nAtoms = 1563;
	rx = new double[nAtoms];
	ry = new double[nAtoms];
	rz = new double[nAtoms];
    int d_parameter = 32;
    int maxBonds = nAtoms/2-d_parameter-1;
    
    double *bondCorrArr, *meanBondCorrArr, *stdBondCorrArr;
    bondCorrArr = new double[nAtoms/2];
    meanBondCorrArr = new double[nAtoms/2];
    stdBondCorrArr = new double[nAtoms/2];
    for (int s = 0; s<nAtoms/2; s++)
    {
        meanBondCorrArr[s] = 0;
        stdBondCorrArr[s] = 0;
    }
    
	for (long iFrame=0;iFrame<nFrames;iFrame++)
	{
        cout<<iFrame<<endl;
		fTraj >> ndum;
        //cout<<ndum<<endl;
		//fTraj >> cdum;
		//fTraj >> cdum;
		//fTraj >> ndum;

		for (long iatom=0;iatom<nAtoms;iatom++)
		{
			fTraj >> cdum >> rx[iatom] >> ry[iatom] >> rz[iatom];
            //cout<<atType[iatom]<<endl;
			//cout<<iatom<<"\t"<<rx[iatom]<<"\t"<<ry[iatom]<<"\t"<<rz[iatom]<<"\t"<<endl;
        }
        //getchar();
        //Compute bond vectors
        for (int s = 0; s<nAtoms/2; s++)
        {
            bondCorrArr[s] = 0;
            for (int iBond=0;iBond<maxBonds;iBond++)
            {
                double dx_iBond = rx[iBond+d_parameter]-rx[iBond];
                double dy_iBond = ry[iBond+d_parameter]-ry[iBond];
                double dz_iBond = rz[iBond+d_parameter]-rz[iBond];
                
                int jBond = iBond+s;
                double dx_jBond = rx[jBond+d_parameter]-rx[jBond];
                double dy_jBond = ry[jBond+d_parameter]-ry[jBond];
                double dz_jBond = rz[jBond+d_parameter]-rz[jBond];
                double bondCorr = (dx_iBond*dx_jBond+dy_iBond*dy_jBond+dz_iBond*dz_jBond)/sqrt((dx_iBond*dx_iBond+dy_iBond*dy_iBond+dz_iBond*dz_iBond)*(dx_jBond*dx_jBond+dy_jBond*dy_jBond+dz_jBond*dz_jBond));
                bondCorrArr[s] += bondCorr;
                //cout<<iBond<<"\t"<<jBond<<"\t"<<s<<"\t"<<bondCorr<<endl;
                //getchar();
            }
            bondCorrArr[s] /= maxBonds;
            meanBondCorrArr[s] += bondCorrArr[s];
            stdBondCorrArr[s] += bondCorrArr[s]*bondCorrArr[s];
        }
    }
    
    ofstream fBondCorr;
    fBondCorr.open("BondCorr.dat");
    for (int s = 0; s<nAtoms/2; s++)
    {
        meanBondCorrArr[s] /= nFrames;
        stdBondCorrArr[s] /= nFrames;
        stdBondCorrArr[s] = stdBondCorrArr[s]-meanBondCorrArr[s]*meanBondCorrArr[s];
        fBondCorr<<s<<"\t"<<meanBondCorrArr[s]<<"\t"<<stdBondCorrArr[s]<<endl;
    }
    fBondCorr.close();
    //Print the matrices
    return 0;
}
