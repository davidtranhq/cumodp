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


## roots_of_unity.mpl

## Computes a generator of the cylic group
## of the units of Z/pZ for a prime p
CyclicGroupGenerator := proc(p)
    local G, alpha, f, fact, i, G1;
    if not isprime(p) then error "arg#1 is not prime"; end if;
    G := GF(p, 1);
    G1 := G:-ConvertIn(1);
    fact := ifactors(p-1);
    alpha := G:-PrimitiveElement();
    for i from 1 to nops(fact[2]) do 
        f := fact[2][i][1];
        if (G:-`^`(alpha,f) = G1) then error "bad generator" end if;
    end do;
    return (G:-ConvertOut(alpha));
end proc;

## Computes a N-th primitive root of unity
## of Z/pZ for a prime p
PrimitiveRootOfUnity := proc(p, N)
     local q, G, alpha, Galpha, Gomega;
     if irem(p-1,N) <> 0 then error "arg#2 does not divide arg#1 - 1"; end if;
     q := iquo(p-1,N);
     G := GF(p, 1);
     alpha := CyclicGroupGenerator(p);
     Galpha := G:-ConvertIn(alpha);
     Gomega :=  G:-`^`(Galpha, q);
     return (G:-ConvertOut(Gomega));
end proc;

## Returns true if omega is a N-th primitive 
## root of unity modulo p
IsPrimitiveRootOfUnity := proc(p, omega, N)
      local G, fact, i, f, Gomega, G1;
      G := GF(p, 1);
      Gomega := G:-ConvertIn(omega);
      G1 := G:-ConvertIn(1);
      if (G:-`^`(Gomega,N) <> G1) then return(false); end if;
      fact := ifactors(N);
      for i from 1 to nops(fact[2]) do 
          f := fact[2][i][1];
          e := iquo(N, f);
          if (G:-`^`(Gomega,e) = G1) then return(false); end if;
      end do;
      return (true);
end proc;

## Returns the N-th roots of beta modulo p
## Assuming that X^N - beta splits into linear factors
## modulo p
NthRoots := proc(p, beta,  M)
       local X, LF, i, f;
       LF := Factors(X^M - beta) mod p;
       for f in LF[2] do 
           if degree(f[1], X) <> 1 then error "Invalid input" end if;
       end do;
       racines := map(eval, LF[2], [X=0]);
       return ([seq(racines[i][1], i=1..nops(racines))]);
end proc;

## Computes a N-th primitive root of unity omega
## of Z/pZ such that omega^(N/ M) is eta and N is a power of M
PrimitiveRootOfUnityAsRootOf := proc(p, eta, M, N)
       local n, racines, L, omega, Taches, tache, racines_de_omega, plus_de_taches;
       if (N = M) then return(eta) end if;
       if (irem(N, M) <> 0) then error "Invalid input" end if;
       racines := NthRoots(p, eta, M);
       L := iquo(N, M^2);
       Taches := [seq([rac, L], rac=racines)];
       while (nops(Taches) <> 0) do
             tache := Taches[1]; Taches := Taches[2..nops(Taches)];
             omega := tache[1]; n := tache[2];
             if (n = 1) then
                if IsPrimitiveRootOfUnity(p, omega, N) then
                   return(omega);
                end if;
             else
                racines_de_omega := NthRoots(p, omega, M);
                plus_de_taches := [seq([rac, iquo(n,M)], rac=racines_de_omega)];
                Taches := [op(plus_de_taches), op(Taches)];
             end if;
       end do;
       error "Invalid input";
end proc;

## Computes the maximum s such that 2^s divides p-1
## Thus 2^s is the maximum size of input vector 
## for a two-way FFT, like our FFT codes
MaximumExponentOftheSizeOfAWayFFT := proc(p)
    local s, a;
    if not isprime(p) then error "arg#1 is not prime"; end if;
    s := 0;
    a := p-1;
    while (irem(a,2)=0) do
        a := iquo(a,2);
        s := s + 1;
    end do;
    return (s);
end proc;

BigPrimeFieldWithSparseRadixRepresentation := proc(k)
    local Results, i, r, p, s;
    Results := [];
    for i from 0 to 62 do
        r := 2^63 + 2^i;
        p := r^k + 1;
        b := 64 * k;
        s := i*k;
        if isprime(p) then 
           Results := [[k, ifactor(2^63) + ifactor(2^i), ifactor(2^s)], op(Results)];
        end if
    end do;
    for i from 0 to 63 do
        r := 2^64 - 2^i;
        p := r^k + 1;
        s := i*k;
        if isprime(p) then 
           Results := [[k, ifactor(2^64) - ifactor(2^i), ifactor(2^s)], op(Results)];
        end if
    end do;
   for i from 0 to 61 do
        r := 2^62 + 2^i;
        p := r^k + 1;
        s := i*k;
        if isprime(p) then 
           Results := [[k, ifactor(2^62) + ifactor(2^i), ifactor(2^s)], op(Results)];
        end if
    end do;
    for i from 0 to 62 do
        r := 2^63 - 2^i;
        p := r^k + 1;
        s := i*k;
        if isprime(p) then 
           Results := [[k, ifactor(2^63) - ifactor(2^i), ifactor(2^s)], op(Results)];
        end if
    end do;

   for i from 0 to 60 do
        r := 2^61 + 2^i;
        p := r^k + 1;
        s := i*k;
        if isprime(p) then 
           Results := [[k, ifactor(2^61) + ifactor(2^i), ifactor(2^s)], op(Results)];
        end if
    end do;
    for i from 0 to 61 do
        r := 2^62 - 2^i;
        p := r^k + 1;
        s := i*k;
        if isprime(p) then 
           Results := [[k, ifactor(2^62) - ifactor(2^i), ifactor(2^s)], op(Results)];
        end if
    end do;

    return (Results); 
end proc;

################
# Simple tests #
################
r := 2^63 + 2^34;
p := r^8 + 1;
s := MaximumExponentOftheSizeOfAWayFFT(p);
if (s <> 272) then error "In MaximumExponentOftheSizeOfAWayFFT" end ;

for i from 1 to 10 do
   omega := PrimitiveRootOfUnity(p,2^i);
   if not IsPrimitiveRootOfUnity(p, omega, 2^i) then 
      error "Bad primitive root of unity";
   end if; 
   if IsPrimitiveRootOfUnity(p, omega^2, 2^i) then 
      error "Bad primitive root of unity";
   end if; 
end do;

NthRoots(p, r, 16);

# for i from 1 to 10 do
#    N := 16^i;
#    omega := PrimitiveRootOfUnityAsRootOf(p, r, 16, N);
#    if not IsPrimitiveRootOfUnity(p, omega, N) then
#       error "Bad primitive root of unity";
#    end;
#    if IsPrimitiveRootOfUnity(p, omega^2, N) then
#       error "Bad primitive root of unity";
#    end;
#    if IsPrimitiveRootOfUnity(p, omega+1, N) then
#       error "Bad primitive root of unity";
#    end;
# end do;

# BigPrimeFieldWithSparseRadixRepresentation(2);
# BigPrimeFieldWithSparseRadixRepresentation(4);
# BigPrimeFieldWithSparseRadixRepresentation(8);
# BigPrimeFieldWithSparseRadixRepresentation(16);
# BigPrimeFieldWithSparseRadixRepresentation(32);
# BigPrimeFieldWithSparseRadixRepresentation(64);
# BigPrimeFieldWithSparseRadixRepresentation(128);

# writeFile :=fopen("roots_of_unity",WRITE);
# fprintf("%d",root);
# print(root);