/* MySimulationBodyMethod.c */
#include "mex.h"
#include "math.h"
#include "time.h"

double StayInBound(double val)
{
    if(val>1)
    {
        val = 1.0;
    }
    else if(val<-1)
    {
        val = -1.0;
    }
    return val;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double U;
    double *x;
	double *lastOpinions;
    double *A;
    double *meanMajOpininos;
	double *followerRatio;
	double *distribution;
	double alpha;
	size_t N,InformedAgentsSize,it,i,j,neighbor,index,randomness,MaximumSimulationSteps;
    double *adj;
    double *leng;
    double cnt,mean,uncertaintyI,uncertaintyNeighbor,si,sj,maxDeg,mu;
	
    if(nrhs != 7)  
    {
        mexErrMsgIdAndTxt("MyToolbox:MyArrayProduct:nrhs","Seven inputs for the function is required not %d .", nrhs);
    }
    if(nlhs != 4)  
    {
        mexErrMsgIdAndTxt("MyToolbox:MyArrayProduct:nlhs","One output for the function is required not %d .", nlhs);      
    }
	
    U = 1;
    x = mxGetPr(prhs[0]);
    N = mxGetM(prhs[0]);
    InformedAgentsSize = mxGetScalar(prhs[1]);
    MaximumSimulationSteps = mxGetScalar(prhs[2]);
    mu = mxGetScalar(prhs[3]);
    A = mxGetPr(prhs[4]);
    randomness = mxGetScalar(prhs[5]);
    alpha = mxGetScalar(prhs[6]);
    
    plhs[0] = mxCreateDoubleMatrix(1,MaximumSimulationSteps,mxREAL);
    meanMajOpininos = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(1,N,mxREAL);
	lastOpinions = mxGetPr(plhs[1]);
	plhs[2] = mxCreateDoubleMatrix(1,N,mxREAL);
	distribution = mxGetPr(plhs[2]);
	plhs[3] = mxCreateDoubleMatrix(1,MaximumSimulationSteps,mxREAL);
    followerRatio = mxGetPr(plhs[3]);
	
    srand(time(0));
    adj = (double *) mxMalloc ((N*N) * sizeof(double));
    leng = (double *) mxMalloc(N * sizeof(double));
    for(i=0;i<N;i++)
    {
        leng[i] = 0;
        for(j=0;j<N;j++)
        {
            if(A[i*N+j] != 0)
            {
                adj[(size_t)(i*N+leng[i])] = (double)j;
                leng[i]++;
            }
        }
    }
	
    maxDeg = 0;
    for(i=0;i<N;i++)
    {
        if(leng[i] > maxDeg)
        {
            maxDeg = leng[i];
        }
		distribution[i] = 0;
    }
	
    for(it=0;it<MaximumSimulationSteps;it++)
    {
        for(i=0;i<N-InformedAgentsSize;i++)
        {            
            if(leng[i] == 0)
            {
                continue;
            }
			
            index = rand() % (size_t)leng[i];
            neighbor = adj[i*N + index];
            si = pow(leng[i],alpha);
            sj = pow(leng[neighbor],alpha);

            uncertaintyI = U;
            /* uncertaintyI = (maxDeg - leng[i])/maxDeg + ((rand() % 1000)/1000) * randomness; */ 
            if(abs(x[i] - x[neighbor]) <= uncertaintyI)
            {
                x[i] = StayInBound(x[i] + mu * (sj/(si+sj)) * (x[neighbor] - x[i]));
            }

            if(neighbor<N)
            {
                uncertaintyNeighbor = U;
                /* uncertaintyNeighbor = (maxDeg - leng[neighbor])/maxDeg + ((rand() % 1000)/1000) * randomness; */  
                if(abs(x[i] - x[neighbor]) <= uncertaintyNeighbor && neighbor < N-InformedAgentsSize)
                {
                    x[neighbor] = StayInBound(x[neighbor] + mu * (si/(si+sj)) * (x[i] - x[neighbor]));
                }
            }
        }
		
		cnt = 0;
        for(i=0;i<N-InformedAgentsSize;i++)
        {
            if(x[i] > 0)
			{
				cnt++;
			}
        }
        cnt/=(N-InformedAgentsSize);
        followerRatio[it] = cnt;
        
        mean = 0;
        for(i=0;i<N-InformedAgentsSize;i++)
        {
            mean += x[i];
        }
        mean/=(N-InformedAgentsSize);
        meanMajOpininos[it] = mean;
		
		for(i=0;i<N;i++)
		{
			distribution[i] += (x[i]/MaximumSimulationSteps);
		}
    }
	
	for(i=0;i<N;i++)
	{
		lastOpinions[i] = x[i];
	}
	
	/* mexWarnMsgIdAndTxt("MyToolbox:MyArrayProduct:nlhs","completed ..."); */
}


 