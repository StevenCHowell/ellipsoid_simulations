//
// Hailiang Zhang
// california institute of technology
//

#include "cudaKernels.h"


/// @par Main functionality
/// Based on the pre-generated random number, and the pre-computed good/bad-move flag, accept/reject a random move.<br>
/// For the accepted random move, the candidate @b M and <i>llk</i>'s will copied to the to-be-updated @b M and <i>llk</i>'s 
/// @param [in] Ns number of samples
/// @param [in] NMparam number of parameters in the @b M matrix
/// @param [in] gR_select vector for the pre-generated random numbers for Metropolis screening (Dimension: Ns)
/// @param [in] rejects vector for the pre-evaluated flags for good/bad moves (Dimension: Ns)
/// @param [in, out] gpnacc the total number of accepted samples
/// @param [in, out] gM the to-be-updated sample (@b M) matrix on GPU (Dimension: NMparam * Ns, Ns is the leading index)
/// @param [in, out] gllkprior the to-be-updated <i>prior-llk</i> vector (Ns) on GPU (Dimension: Ns)
/// @param [in, out] gllkdata the to-be-updated <i>data-llk</i> vector (Ns) on GPU (Dimension: Ns)
/// @param [in, out] gllkpost the to-be-updated <i>posterior-llk</i> vector (Ns) on GPU (Dimension: Ns)
/// @param [in, out] gM_candidate the candidate sample (@b M) matrix on GPU (Dimension: NMparam * Ns)
/// @param [in, out] gllkprior_candidate the candidate <i>prior-llk</i> vector (on GPU (Dimension: Ns)
/// @param [in, out] gllkdata_candidate the candidate <i>data-llk</i> vector on GPU (Dimension: Ns)
/// @param [in, out] gllkpost_candidate the candidate <i>posterior-llk</i> vector on GPU (Dimension: Ns)
/// @note
/// <i>notes to developers:</i><br>
/// <b>cuda kernel layout:</b> gridDim--Ns/BLOCKDIM, blockDim--BLOCKDIM
__global__ void
altar::bayesian::cuda::cudaKernels_Metropolis::
cudaUpdateSample(const int Ns, const int NMparam, TYPE * const gR_select, const int * const rejects, int * gpnacc,
				TYPE * const gM, TYPE * const gllkprior, TYPE * const gllkdata, TYPE * const gllkpost,
				const TYPE * const gM_candidate, const TYPE * const gllkprior_candidate, const TYPE * const gllkdata_candidate, const TYPE * const gllkpost_candidate)
{
    int sample = blockIdx.x*blockDim.x + threadIdx.x;
    if (sample >= Ns) return;
	//if (sample <10 || rejects[sample]==1) printf("Sample: %d llk: %8.5f gLLKpost[sample]: %8.5f difference: %8.5f rnd: %8.5f reject? %d\n",sample,gllkpost_candidate[sample],gllkpost[sample],(gllkpost_candidate[sample]-gllkpost[sample]),log(gR_select[sample]), rejects[sample]);
	int i;
	//if (sample==0&&rejects[sample]==1) { printf("rejected:\n"); for (i=0; i<NMparam; i++) printf("%f ",gM_candidate[i*Ns+sample]); printf("\n");}
    if (rejects[sample]==0 && gllkpost_candidate[sample]-gllkpost[sample]>=log(gR_select[sample]))
    {   
        //printf("Accept sample: %d llk: %8.5f gLLKpost[sample]: %8.5f difference: %8.5f rnd: %8.5f NMparam: %d\n",sample,gllkpost_candidate[sample],gllkpost[sample],(gllkpost_candidate[sample]-gllkpost[sample]),log(gR_select[sample]), NMparam);
        gllkprior[sample] = gllkprior_candidate[sample];
        gllkdata[sample] = gllkdata_candidate[sample];
        gllkpost[sample] = gllkpost_candidate[sample];
        for (i=0; i<NMparam; i++) gM[i*Ns+sample] = gM_candidate[i*Ns+sample];
        atomicAdd(gpnacc,1);
    }   
}
