
#ifndef LIST_PLAIN_MUL
#define LIST_PLAIN_MUL
#include "subproduct_tree.h"
const sfixn Tmul = 512;
//__global__ void listPlainMulGpu2( sfixn *Mgpu1, sfixn *Mgpu2 , sfixn length_poly, sfixn poly_on_layer, sfixn threadsForAmul, sfixn mulInThreadBlock, sfixn p);
__global__ void listPlainMulGpu(sfixn *M_dev, sfixn start_offset, sfixn length_poly, sfixn num_poly, sfixn threadsForAmul, sfixn mulInThreadBlock, sfixn p);
#endif
