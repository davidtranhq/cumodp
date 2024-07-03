with(RegularChains):
read "prime.txt";
read "poly1.txt";
read "poly2.txt";
resultant3 := Resultant(F1,F2,c) mod p:
read "resultant1.txt";
read "resultant2.txt";
check1 := Expand(resultant1 - resultant2) mod p;
check2 := Expand(resultant3 - resultant1) mod p;
check3 := Expand(resultant3 - resultant2) mod p;



