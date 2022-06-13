 **AlphaSparse** 



AlphaSparse sparse matrix algorithm library interface description for DCU



This document provides a detailed stage function interface description of the AlphaSparse sparse matrix algorithm library for DCU, including mv, mm, etc. The function will perform the corresponding algorithm processing according to the sparse matrix format. The specific format of the sparse matrix will not be reflected in the interface name and the reception parameters.



1. **Level 1**



- - **AXPYI**

  - - axpyi(ALPHA_INT nnz,
    - const ALPHA_Number alpha,
    - const ALPHA_Number *x_val,
    - const ALPHA_INT *x_ind,
    - ALPHA_Number *y



The axpyi function multiplies the sparse vector 𝑥 with scalar 𝛼 and adds the result to the dense vector 𝑦, such that 

​	𝑦 := 𝑦 + 𝛼 · 𝑥 

| Input parameters | Description                                                  |
| ---------------- | ------------------------------------------------------------ |
| nnz              | number of non-zero entries of vector 𝑥                       |
| alpha            | Scalar value alpha                                           |
| x_val            | array of nnz elements containing the values of 𝑥             |
| x_indx           | array of nnz elements containing the indices of the non-zero values of 𝑥 |
| y                | array of values in dense format                              |







- - **AXPBY**



- - - axpyi(ALPHA_INT nnz,
    - const ALPHA_Number alpha,
    - const ALPHA_Number beta,
    - const ALPHA_Number *x_val,
    - const ALPHA_INT *x_ind,
    - ALPHA_Number *y





The aspby multiplies the sparse vector 𝑥 with scalar 𝛼 and adds the result to the dense vector 𝑦 that is multiplied with scalar 𝛽, such that 

​	𝑦 := 𝛼 · 𝑥 + 𝛽 · 𝑦 

| Input parameters | Description                                      |
| ---------------- | ------------------------------------------------ |
| alpha            | Scalar value alpha                               |
| x                | array of nnz elements containing the values of 𝑥 |
| Beta             | Scalar value beta                                |
| y                | array of values in dense format                  |





- - **DOTI**



- - - doti(ALPHA_INT nnz,
    - const ALPHA_Number *x_val,
    - const ALPHA_INT *x_ind,
    - const ALPHA_Number *y,
    - ALPHA_Number *result



The doti executes dot product operation of compressed real number vector and full storage real number vector and return the result value:



​	res = x[0]*y[indx[0]] + x[1]*y[indx[1]] +... + x[nz-1]*y[indx[nz-1]]



| Input parameters | Description                                                  |
| ---------------- | ------------------------------------------------------------ |
| nnz              | number of non-zero entries of vector 𝑥                       |
| Result           | pointer to the result, can be host or device memory          |
| x_val            | array of nnz elements containing the values of 𝑥             |
| x_indx           | array of nnz elements containing the indices of the non-zero values of 𝑥 |
| y                | array of values in dense format                              |





- - **DOTCI**



- - - dotci(ALPHA_INT nnz,
    - const ALPHA_Number *x_val,
    - const ALPHA_INT *x_ind,
    - const ALPHA_Number *y,
    - ALPHA_Number *result



dotci computes the dot product of the complex conjugate sparse vector 𝑥 with the dense vector 𝑦, such that

​	result += conj(x_val[i]) * y[x_ind[i]];



| Input parameters | Description                                                  |
| ---------------- | ------------------------------------------------------------ |
| nnz              | number of non-zero entries of vector 𝑥                       |
| Result           | pointer to the result, can be host or device memory          |
| x_val            | array of nnz elements containing the values of 𝑥             |
| x_indx           | array of nnz elements containing the indices of the non-zero values of 𝑥 |
| y                | array of values in dense format                              |



- - **GTHR**



- - - gthr(ALPHA_INT nnz,
    - const ALPHA_Number *y,
    - ALPHA_Number *x_val,
    - const ALPHA_INT *x_ind



The gthr gathers the elements that are listed in x_ind from the dense vector 𝑦 and stores them in the sparse vector 𝑥. 

​	x_val[i] = y[x_ind[i]]





| Input parameters | Description                                                  |
| ---------------- | ------------------------------------------------------------ |
| nnz              | number of non-zero entries of vector 𝑥                       |
| x_val            | array of nnz elements containing the values of 𝑥             |
| x_indx           | array of nnz elements containing the indices of the non-zero values of 𝑥 |
| y                | array of values in dense format                              |





- - **GTHRZ**



- - - gthrz(ALPHA_INT nnz,
    - ALPHA_Number *y,
    - ALPHA_Number *x_val,
    - const ALPHA_INT *x_ind)



The gthrz executes by index of gathering the elements of a full storage vector into the compressed vector format, and zeroing the elements at the corresponding positions in the original vector:



​	for(i = 0; i < nnz; ++i)

​	{

​		x_val[i]  = y[x_ind[i]];

  	y[x_ind[i]] = 0;

​	}



| Input parameters | Description                                                  |
| ---------------- | ------------------------------------------------------------ |
| nnz              | number of non-zero entries of vector 𝑥                       |
| x_val            | array of nnz elements containing the values of 𝑥             |
| x_indx           | array of nnz elements containing the indices of the non-zero values of 𝑥 |
| y                | array of values in dense format                              |



- - **ROTI**



- - - roti(ALPHA_INT nnz,
    - ALPHA_Number *x_val,
    - const ALPHA_INT *x_ind,
    - ALPHA_Number *y,
    - ALPHA_Number c,
    - ALPHA_Number s

The roti applies the Givens rotation matrix 𝐺 to the sparse vector 𝑥 and the dense vector 𝑦, where: 

​	



​	

​	for(i = 0; i < nnz; ++i) { 

  	x_tmp = x_val[i];

  	y_tmp = y[x_ind[i]];

  	x_val[i]  = c * x_tmp + s * y_tmp;

  	y[x_ind[i]] = c * y_tmp - s * x_tmp;

   }





| Input parameters | Description                                                  |
| ---------------- | ------------------------------------------------------------ |
| nnz              | Number of elements in vectors x and indx                     |
| x_val            | Array of nnz elements containing the non-zero values of 𝑥    |
| x_ind            | array of nnz elements containing the indices of the non-zero values of 𝑥 |
| y                | Store as an array, the length is at least max(indx[i])       |
| c                | pointer to the cosine element of 𝐺, can be on host or device |
| s                | pointer to the sine element of 𝐺, can be on host or device   |



- - **SCTR**





- - - gthr(ALPHA_INT nnz,
    - const ALPHA_Number *x_val,
    - const ALPHA_INT *x_ind,
    - ALPHA_Number *y



The sctr catters the elements that are listed in x_ind from the sparse vector 𝑥 into the dense vector 𝑦. Indices of 𝑦 that are not listed in x_ind remain unchanged. 



​	for(i = 0; i < nnz; ++i) { 

  	y[x_ind[i]] = x_val[i];

   }





| Input parameters | Description                                                  |
| ---------------- | ------------------------------------------------------------ |
| nz               | Number of elements in vectors x and indx                     |
| x                | Store as an array,,length is at least nz, contains the vector converted to full storage |
| indx             | Given the element index of x that will be scattered,Store as an array, length is at least nz |
| y                | Store as an array, length is at least max(indx[i]), Contains the updated vector element value |





1. **Level 2**



1. **dcu_gemv**

- - alphasparse_layout_t layout,
  - ALPHA_INT mb,
  - ALPHA_INT nb,
  -  ALPHA_INT nnzb,
  - const ALPHA_Number alpha,
  - const ALPHA_Number *bsr_val,
  - const ALPHA_INT *bsr_row_ptr,
  - const ALPHA_INT *bsr_col_ind,
  - ALPHA_INT bs,
  - const ALPHA_Number *x,
  - const ALPHA_Number beta,
  - ALPHA_Number *y



The dcu_gemv function performs the operation of multiplying a sparse matrix and a dense vector:



y := alpha * op(A) * x + beta * y



Alpha and beta are scalar values, A is a sparse matrix with k rows and m columns, x and y are vectors. "?" indicates the data format, which corresponds to the DATATYPE in the interface, s corresponds to float, d corresponds to double, and c corresponds to float complex, which is a single-precision complex number, and z corresponds to a double complex, which is a double-precision complex number. This function stores the output result in the vector y. The input parameters of the function are shown in the table:







| Input  parameters | Description                                                  |
| ----------------- | ------------------------------------------------------------ |
| operation         | For specific operations on the input matrix, there are the following options: OPERATION_NON_TRANSPOSE, no transposition, op(A) = A OPERATION_TRANSPOSE, transpose, op(A) = AT OPERATION_CONJUGATE_TRANSPOSE, ConjugationTranspose, op(A) = AH |
| mb                | number of block rows of the sparse BSR matrix                |
| nb                | number of block columns of the sparse BSR matrix             |
| nnzb              | number of non-zero blocks of the sparse BSR matrix           |
| bsr_val           | array of nnzb blocks of the sparse BSR matrix                |
| bsr_row_ptr       | array of mb+1 elements that point to the start of every block row of the sparse BSR matrix. |
| bsr_col_ind       | array of nnzb elements containing the block column indices of the sparse BSR matrix. |
| block_dim         | block dimension of the sparse BSR matrix.                    |
| alpha             | Scalar value alpha                                           |
| A                 | Data structure of sparse matrix                              |
| x                 | Dense vector x, stored as an array, if no transpose operation is performed on matrix A, the length is at least the number of columns of matrix A |
| beta              | Scalar value beta                                            |
| y                 | Dense vector y, stored as an array, if no transpose operation is performed on matrix A, the length is at least the number of rows of matrix A |







1. **dcu_gemv_?_coo**

2. 1. ALPHA_INT nnz,
   2. ​    ALPHA_INT loops,
   3.  ALPHA_Number alpha,
   4. const ALPHA_INT *__restrict__ coo_row_ind,
   5.  const ALPHA_INT *__restrict__ coo_col_ind,
   6.  const ALPHA_Number *__restrict__ coo_val,
   7. const ALPHA_Number *__restrict__ x,
   8.  ALPHA_Number *__restrict__ y,
   9. ALPHA_INT *__restrict__ row_block_red,
   10. ​	10.ALPHA_Number *__restrict__ val_block_red



The dcu_gemv_?_bsr function performs the operation of multiplying a sparse matrix and a dense vector:



**y := alpha \* op(A) \* x + beta \* y**



Alpha and beta are scalar values, A is a sparse matrix with k rows and m columns, x and y are vectors. "?" indicates the data format, which corresponds to the DATATYPE in the interface, s corresponds to float, d corresponds to double, and c corresponds to float complex, which is a single-precision complex number, and z corresponds to a double complex, which is a double-precision complex number. This function stores the output result in the vector y. The input parameters of the function are shown:



| Input  parameters | Description                                                  |
| ----------------- | ------------------------------------------------------------ |
| operation         | For specific operations on the input matrix, there are the following options: OPERATION_NON_TRANSPOSE, no transposition, op(A) = A OPERATION_TRANSPOSE, transpose, op(A) = AT OPERATION_CONJUGATE_TRANSPOSE, ConjugationTranspose, op(A) = AH |
| mb                | number of block rows of the sparse BSR matrix                |
| nb                | number of block columns of the sparse BSR matrix             |
| nnzb              | number of non-zero blocks of the sparse BSR matrix           |
| bsr_val           | array of nnzb blocks of the sparse BSR matrix                |
| bsr_row_ptr       | array of mb+1 elements that point to the start of every block row of the sparse BSR matrix. |
| bsr_col_ind       | array of nnzb elements containing the block column indices of the sparse BSR matrix. |
| block_dim         | block dimension of the sparse BSR matrix.                    |
| alpha             | Scalar value alpha                                           |
| A                 | Data structure of sparse matrix                              |
| x                 | Dense vector x, stored as an array, if no transpose operation is performed on matrix A, the length is at least the number of columns of matrix A |
| beta              | Scalar value beta                                            |
| y                 | Dense vector y, stored as an array, if no transpose operation is performed on matrix A, the length is at least the number of rows of matrix A |



**3.** **level 3**





- - ​       alphasparse_layout_t dir,
  - ​        ALPHA_INT mb,
  - ​        ALPHA_INT n,
  - ​        ALPHA_INT kb,
  - ​        ALPHA_INT nnzb,
  - ​        const ALPHA_Number alpha,
  - ​        const ALPHA_Number *bsr_val,
  - ​        const ALPHA_INT *bsr_row_ptr,
  - ​        const ALPHA_INT *bsr_col_ind,
  - ​        ALPHA_INT block_row_dim,
  - ​        ALPHA_INT block_col_dim,
  - ​        const ALPHA_Number *x,
  - ​        ALPHA_INT ldx,
  - ​        const ALPHA_Number beta,
  - ​        ALPHA_Number *y,
  - ​        ALPHA_INT ldy





multiplies the scalar 𝛼 with a sparse 𝑚𝑏 × 𝑘𝑏 matrix 𝐴, defined in BSR storage format, and the dense 𝑘 × 𝑛 matrix 𝐵 (where 𝑘 = 𝑐𝑜𝑙𝑏𝑙𝑜𝑐𝑘_𝑑𝑖𝑚 × 𝑘𝑏) and adds the result to the dense 𝑚 × 𝑛 matrix 𝐶 (where 𝑚 = 𝑟𝑜𝑤𝑏𝑙𝑜𝑐𝑘_𝑑𝑖𝑚 × 𝑚𝑏) that is multiplied by the scalar 𝛽, such that 



​	𝐶 **:=**𝛼**·**𝑜𝑝**(**𝐴**)·**𝑜𝑝**(**𝐵**)+**𝛽**·**𝐶**,** 







| Input  parameters | Description                                                  |
| ----------------- | ------------------------------------------------------------ |
| operation         | For specific operations on the input matrix, there are the following options: OPERATION_NON_TRANSPOSE, no transposition, op(A) = A OPERATION_TRANSPOSE, transpose, op(A) = AT OPERATION_CONJUGATE_TRANSPOSE, ConjugationTranspose, op(A) = AH |
| mb                | number of block rows of the sparse GEneral BSR matrix 𝐴      |
| Kb                | number of columns of the dense matrix 𝑜𝑝(𝐵) and 𝐶            |
| n                 | number of block columns of the sparse matrix                 |
| nnzb              | number of non-zero blocks of the sparse matrix               |
| bsr_val           | array of nnzb blocks of the sparse matrix                    |
| bsr_row_ptr       | array of mb+1 elements that point to the start of every block row of the sparse BSR matrix. |
| bsr_col_ind       | array of nnzb elements containing the block column indices of the sparse BSR matrix. |
| block_dim         | block dimension of the sparse BSR matrix.                    |
| alpha             | Scalar value alpha                                           |
| A                 | Data structure of sparse matrix                              |
| ldx               | leading dimension of 𝐵, must be at least max(1,𝑘) ( 𝑜𝑝(𝐵) == 𝐵) where 𝑘 = 𝑏𝑙𝑜𝑐𝑘_𝑑𝑖𝑚 × 𝑘𝑏, max (1, 𝑛) otherwise. |
| x                 | Dense vector x, stored as an array, if no transpose operation is performed on matrix A, the length is at least the number of columns of matrix A |
| beta              | Scalar value beta                                            |
| y                 | Dense vector y, stored as an array, if no transpose operation is performed on matrix A, the length is at least the number of rows of matrix A |
| ldy               | leading dimension of 𝐶, must be at least max(1,𝑚) ( 𝑜𝑝(𝐴) == 𝐴) where 𝑚 = 𝑏𝑙𝑜𝑐𝑘_𝑑𝑖𝑚 × 𝑚𝑏, max (1, 𝑘) where 𝑘 = 𝑏𝑙𝑜𝑐𝑘_𝑑𝑖𝑚 × 𝑘𝑏 otherwise |