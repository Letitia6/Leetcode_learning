# How to add kernel to AlphaSparse:


## Using sp2m as an example, the steps are as follows:

### 1. Add CSR、CSC、COO and BSR to these four def_x.h

../kernel/def_c.h

../kernel/def_d.h

../kernel/def_s.h

../kernel/def_z.h 



### 2. Add CSR、CSC、COO and BSR to these kernel_x_x.h. Using csr as an example:

../kernel/kernel_csr_c.h

../kernel/kernel_csr_d.h

../kernel/kernel_csr_s.h

../kernel/kernel_csr_z.h

 

### 3. Add CSR、CSC、COO and BSR  to ../src/core/op_ /alphasparse_sparse_sp2m.c.

 

### 4. Add kernel code in /kernel/x86_64/level3/sp2m.

For instsnce /src/kernel/x86_64/level3/sp2m/csr/generel/general/A/B/sp2m_csr_gen_gen_A_B_x_FM.c
