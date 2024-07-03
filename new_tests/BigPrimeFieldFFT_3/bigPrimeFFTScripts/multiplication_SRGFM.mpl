/* This file is part of the CUMODP library

    CUMODP is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CUMODP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CUMODP.  If not, see <http://www.gnu.org/licenses/>.

    Copyright: Sardar Anisul Haque <shaque4@uwo.ca>
               Xin Li <xli.software@gmail.com>
               Farnam Mansouri <mansouri.farnam@gmail.com>
               Davood Mohajerani <dmohajer@uwo.ca>
               Marc Moreno Maza  <moreno@csd.uwo.ca>
               Wei Pan <wei.pan@intel.com>
               Ning Xie <nxie6@csd.uwo.ca>
*/


number_of_plus := proc(f)
    local n, operands, results;
    if f=0 then return (0); end if;
    if type(f,symbol) then return (0); end if;
    if type(f,numeric) then return (0); end if;
    if type(f, `+`) then n := (nops(f) -1) else n := 0; end if;
    operands := [op(f)];
    results := map(number_of_plus, operands);
    return (n + sum(results[i], i=1..nops(results)));
end proc;

number_of_times := proc(f)
    local n, operands, results;
    if f=0 then return (0); end if;
    if type(f,symbol) then return (0); end if;
    if type(f,numeric) then return (0); end if;
    if type(f, `*`) then n := (nops(f) -1) else n := 0; end if;
    operands := [op(f)];
    results := map(number_of_times, operands);
    return (n + sum(results[i], i=1..nops(results)));
end proc;


number_of_additions := proc(poly, r)
     local d, i, f, adds;
     d := degree(poly, r);
     adds := 0;
     for i from 0 to d do
         f := coeff(poly, r, i);
         adds := adds + number_of_plus(f);
     end do;
     return adds;
end proc;

number_of_multiplications := proc(poly, r)
     local d, i, f, muls;
     d := degree(poly, r);
     muls := 0;
     for i from 0 to d do
         f := coeff(poly, r, i);
         muls := muls + number_of_times(f);
     end do;
     return muls;
end proc;



x := x_0 + x_1 * r + x_2 * r^2 + x_3 * r^3 + x_4 * r^4 + x_5 * r^5 + x_6 * r^6 + x_7 * r^7;
y := y_0 + y_1 * r + y_2 * r^2 + y_3 * r^3 + y_4 * r^4 + y_5 * r^5 + y_6 * r^6 + y_7 * r^7;


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%% FIRST VIEWING r AS A VARIABLE %%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%
#%% Long multiplication %%
#%%


xy_long := expand(x*y);

[seq(coeff(xy_long, r, i), i=0..14)];

number_of_additions_xy_long := number_of_additions(xy_long, r);
number_of_multiplications_xy_long := number_of_multiplications(xy_long, r);

number_of_arithmetic_ops_xy_long := number_of_additions_xy_long + number_of_multiplications_xy_long;


#%%
#%% Karatsuba %%
#%%

X_0 := x_0 + x_1 * r + x_2 * r^2 + x_3 * r^3; X_1 := x_4 + x_5 * r + x_6 * r^2 + x_7 * r^3;
Y_0 := y_0 + y_1 * r + y_2 * r^2 + y_3 * r^3; Y_1 := y_4 + y_5 * r + y_6 * r^2 + y_7 * r^3;

A := expand(X_0 * Y_0); 
B := expand(X_1 * Y_1); 
C := expand((X_0 + X_1) * (Y_0 + Y_1));
E := C - A - B;
K := A  + E * r^4 + B * r^8;

expand((X_0 + r^4 * X_1) * (Y_0 + r^4 * Y_1) - xy_long);
expand(K - xy_long);

number_of_additions_xy_Karatsuba := number_of_additions(A, r) + 
                                    number_of_additions(B, r) +
                                    degree(X_0 + X_1, r) + 1 + degree(Y_0 + Y_1, r) + 1 +
                                    number_of_additions(B, r) + 
                                    2 * (degree(E, r) + 1) +
                                    degree(E, r) + 1;

number_of_multiplications_xy_Karatsuba := number_of_multiplications(A, r) +
                                          number_of_multiplications(B, r) +
                                          number_of_additions(B, r);

number_of_arithmetic_ops_xy_Karatsuba := number_of_additions_xy_Karatsuba + number_of_multiplications_xy_Karatsuba;