In do_setup.m
N: discretization level of sample W
numThetan: number of angles

To solve the system of equations by backslash: x=L\b(:);
norm(x-W(:)) gives the reconstruction error comparing to ground truth W; 

XTM_Tensor.m shows how L is constructed. 
