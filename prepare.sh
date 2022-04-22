## This script prepares matrix for SVD. No parallel code here.
gfortran -o arrayconvert read.f
cd Res_P
mpiifort -O3 -axCORE-AVX2,AVX -ipo -o svd_p_low svd_p_low.f ../pdlawrite5.f ../pdlaread.f -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_thread -lmpifort -liomp5
mpiifort -O3 -axCORE-AVX2,AVX -ipo -o mdinvs_plow6000 mdinvs_plow6000.f ../pdlawrite5.f ../pdlaread.f -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_thread -lmpifort -liomp5
mpiifort -O3 -axCORE-AVX2,AVX -ipo -o rconst_p_low6000 rconst_p_low6000.f ../pdlawrite5.f ../pdlaread.f -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_thread -lmpifort -liomp5
mpiifort -O3 -axCORE-AVX2,AVX -ipo -o resopr_p_low_full6000 resopr_p_low_full6000.f ../pdlawrite5.f ../pdlaread.f -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_thread -lmpifort -liomp5
cd ../Res_S
mpiifort -O3 -axCORE-AVX2,AVX -ipo -o svd_sonly svd_sonly.f ../pdlawrite5.f ../pdlaread.f -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_thread -lmpifort -liomp5
mpiifort -O3 -axCORE-AVX2,AVX -ipo -o mdinvs_sonly6000 mdinvs_sonly6000.f ../pdlawrite5.f ../pdlaread.f -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_thread -lmpifort -liomp5
mpiifort -O3 -axCORE-AVX2,AVX -ipo -o rconst_sonly_low6000 rconst_sonly_low6000.f ../pdlawrite5.f ../pdlaread.f -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_thread -lmpifort -liomp5
mpiifort -O3 -axCORE-AVX2,AVX -ipo -o resopr_sonly_low_full6000 resopr_sonly_low_full6000.f ../pdlawrite5.f ../pdlaread.f -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_thread -lmpifort -liomp5
cd ..
# python run.py
