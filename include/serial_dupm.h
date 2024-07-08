#ifndef _SERIAL_DUPM_H_
#define _SERIAL_DUPM_H_

__host__ __device__ void serial_dupm(
        sfixn *polynomial1,
        sfixn *polynomial2,
        sfixn polynomial1_degree,
        sfixn polynomial2_degree,
        sfixn *polynomial_product
);

#endif
