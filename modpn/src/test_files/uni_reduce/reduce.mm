with(RegularChains):
with(ChainTools):
read "prime.txt";
read "poly.txt";
read "rc.txt";
RC := Chain(rc, Empty(R), R):
outputMaple := NormalForm(f,RC,R):
read "outputDF.txt";
read "outputBULL.txt";
check1 := Expand(outputDF - outputBULL) mod p;
check2 := Expand(outputMaple - outputBULL) mod p;
check3 := Expand(outputDF - outputMaple) mod p;





