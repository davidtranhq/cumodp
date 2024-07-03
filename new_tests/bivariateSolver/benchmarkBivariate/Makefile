
maple_benchmark: maple_triangular maple_grobner
	cat results_maple_* > results_maple.txt
	less results_maple.txt

maple_triangular:
	maple benchmarkMaple_triangular.mpl > /dev/null

maple_grobner:
	maple benchmarkMaple_grobner.mpl > /dev/null

sage_benchmark: sage_triangular sage_grobner sage_variety
	cat result_sage_* > results_sage.txt
	less results_sage.txt

sage_triangular:
	sage benchmarkSage_triangular.sage

sage_grobner:
	sage benchmarkSage_grobner.sage

sage_variety:
	sage benchmarkSage_variety.sage

clean:
	@ rm -f *.sage.py
	@ rm -f results_*