read "determinant.input":

condensation := proc(M1, n,  P)
local det1;

	det1 := LinearAlgebra:-Determinant(M1) mod P;
	

return(det1);
end proc:

condensation(mat, n, P);

