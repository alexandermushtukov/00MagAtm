objects =   ./obj/c.o \
            ./obj/RT.o


objects_math =	./obj/mf_algebra_polinom.o \
                ./obj/mf_bessel.o \
		./obj/mf_c.o \
		./obj/mf_equations.o \
		./obj/mf_ermit_pol.o \
		./obj/mf_erm_link.o \
		./obj/mf_fact.o \
		./obj/mf_fi_function.o \
		./obj/mf_gamma-function.o \
		./obj/mf_gener_fun.o \
		./obj/mf_g.o \
		./obj/mf_integrals.o \
		./obj/mf_integral_2d.o \
		./obj/mf_int_gauss_100.o \
		./obj/mf_int_gauss.o \
		./obj/mf_intGauss.o \
		./obj/mf_intGauss_p.o \
		./obj/mf_geom_2d.o \
		./obj/mf_geom_3d.o \
		./obj/mf_int_G.o \
		./obj/mf_intSimps.o \
		./obj/mf_intSimp_p.o \
		./obj/mf_interpol.o \
		./obj/mf_ksi.o \
		./obj/mf_laguer_exp.o \
		./obj/mf_laguer.o \
		./obj/mf_laguer_gen.o \
		./obj/mf_laguer_n.o \
		./obj/mf_matrix_matrix.o \
		./obj/mf_numbers.o \
		./obj/mf_vector_sk_mult.o \
		./obj/mf_coord_sys.o \
		./obj/mf_for_mass.o \
		./obj/mf_random.o \
		./obj/mf_statistics.o \
		./obj/mf_intExp.o

 #\
		#./obj/mf_bessel_VF.o \
		#./obj/mf_expmag_VF.o \
		#./obj/mf_tcrUL_VF.o


objects_phys =	./obj/PhFun_eMaxvell.o \
		./obj/PhFun_blackbody.o \
		./obj/PhFun_Compton.o \
		./obj/PhFun_ComptonBF_H.o \
		./obj/PhFun_ComptonBF_cons_laws.o \
		./obj/PhFun_ComptonBF_CrossSection.o \
		./obj/PhFun_MonteCarloCS_database.o \
		./obj/PhFun_NormalWaves.o \
		./obj/PhFun_Landau_Level_width.o \
		./obj/PhFun_Cyclotron_absorb.o \
		./obj/PhFun_magnetic_absorb.o \
		./obj/PhFun_PlasmaDispFun.o \
		./obj/PhFun_ComptonBF_Stokes.o \
		./obj/PhFun_Gaunt_factors.o \
		./obj/PhFun_MagneticDipole.o \
		./obj/PhFun_Polarization.o \
		./obj/PhFun_dielec_tensor.o \
		./obj/PhFun_SpR.o \
        ./obj/PhFun_GR_app.o #\


objects_astro =	./obj/af_ref_frame.o \
		./obj/af_mag_accretion.o \
		./obj/af_orbit_bin.o \
		./obj/af_RT_layer_an.o \
		./obj/af_NS_atm_structure.o


c : $(objects_math) $(objects_phys) $(objects_astro) $(objects)
	gfortran -fopenmp -o c $(objects_math) $(objects_phys)  $(objects_astro) $(objects)

./obj/c.o : c.f90
	gfortran -fopenmp -c -o ./obj/c.o c.f90
./obj/RT.o : RT.f90
	gfortran -c -o ./obj/RT.o RT.f90
./obj/RT_NonPolarized.o : RT_NonPolarized.f90
	gfortran -c -o ./obj/RT_NonPolarized.o RT_NonPolarized.f90
./obj/RT_NonPolarized2.o : RT_NonPolarized2.f90
	gfortran -c -o ./obj/RT_NonPolarized2.o RT_NonPolarized2.f90
./obj/RT_NonPolarized3.o : RT_NonPolarized3.f90
	gfortran -c -o ./obj/RT_NonPolarized3.o RT_NonPolarized3.f90
./obj/RT_Polarized.o : RT_Polarized.f90
	gfortran -c -o ./obj/RT_Polarized.o RT_Polarized.f90
./obj/RT_Polarized_2.o : RT_Polarized_2.f90
	gfortran -c -o ./obj/RT_Polarized_2.o RT_Polarized_2.f90
./obj/CrossSection.o : CrossSection.f90
	gfortran -c -o ./obj/CrossSection.o CrossSection.f90
./obj/s_for_mass.o : s_for_mass.f90
	gfortran -c -o ./obj/s_for_mass.o s_for_mass.f90
./obj/RT_timing.o : RT_timing.f90
	gfortran -c -o ./obj/RT_timing.o RT_timing.f90
./obj/mf.o : mf.f90
	gfortran -c -o ./obj/mf.o mf.f90

./obj/RT_MonteCarlo.o : RT_MonteCarlo.f90
	gfortran  -fopenmp -c -o ./obj/RT_MonteCarlo.o RT_MonteCarlo.f90


./obj/RT_MonteCarlo_2Layers.o : RT_MonteCarlo_2Layers.f90
	gfortran -c -o ./obj/RT_MonteCarlo_2Layers.o RT_MonteCarlo_2Layers.f90
./obj/RT_MonteCarloCS_v0.o : RT_MonteCarloCS_v0.f90
	gfortran -c -o ./obj/RT_MonteCarloCS_v0.o RT_MonteCarloCS_v0.f90
./obj/RT_MonteCarloCS_v1.o : RT_MonteCarloCS_v1.f90
	gfortran -c -o ./obj/RT_MonteCarloCS_v1.o RT_MonteCarloCS_v1.f90
./obj/RT_MonteCarloCS_v2.o : RT_MonteCarloCS_v2.f90
	gfortran -c -o ./obj/RT_MonteCarloCS_v2.o RT_MonteCarloCS_v2.f90


./obj/RT_MonteCarloCS_ave.o : RT_MonteCarloCS_ave.f90
	gfortran -fopenmp -c -o ./obj/RT_MonteCarloCS_ave.o RT_MonteCarloCS_ave.f90
./obj/RT_NS_acc_channel_DB.o : RT_NS_acc_channel_DB.f90
	gfortran -fopenmp -c -o ./obj/RT_NS_acc_channel_DB.o RT_NS_acc_channel_DB.f90

./obj/PhFun_MonteCarloCS_database.o : ../PhFun/PhFun_MonteCarloCS_database.f90
	gfortran -fopenmp -c -o ./obj/PhFun_MonteCarloCS_database.o ../PhFun/PhFun_MonteCarloCS_database.f90
./obj/RT_MonteCarloCS_v2add.o : RT_MonteCarloCS_v2add.f90
	gfortran -c -o ./obj/RT_MonteCarloCS_v2add.o RT_MonteCarloCS_v2add.f90
./obj/RT_MonteCarloCS_database.o : RT_MonteCarloCS_database.f90
	gfortran -fopenmp -c -o ./obj/RT_MonteCarloCS_database.o RT_MonteCarloCS_database.f90
./obj/RT_NS_illumination.o : RT_NS_illumination.f90
	gfortran -c -o ./obj/RT_NS_illumination.o RT_NS_illumination.f90
./obj/RT_NS_illumination_v2_DB.o : RT_NS_illumination_v2_DB.f90
	gfortran -c -o ./obj/RT_NS_illumination_v2_DB.o RT_NS_illumination_v2_DB.f90
./obj/new_AC.o : new_AC.f90
	gfortran -c -o ./obj/new_AC.o new_AC.f90
./obj/RT_make_DB_seedBBph.o : RT_make_DB_seedBBph.f90
	gfortran -c -o ./obj/RT_make_DB_seedBBph.o RT_make_DB_seedBBph.f90

./obj/RT_MonteCarlo_MultiLayers.o : RT_MonteCarlo_MultiLayers.f90
	gfortran -g -fcheck=all -Wall -Wextra -O0 -fopenmp -c -o ./obj/RT_MonteCarlo_MultiLayers.o RT_MonteCarlo_MultiLayers.f90


./obj/af_NS_atm_structure.o : ../AstroFun/af_NS_atm_structure.f90
	gfortran -c -o ./obj/af_NS_atm_structure.o ../AstroFun/af_NS_atm_structure.f90
./obj/af_ref_frame.o : ../AstroFun/af_ref_frame.f90
	gfortran -c -o ./obj/af_ref_frame.o ../AstroFun/af_ref_frame.f90
./obj/af_orbit_bin.o : ../AstroFun/af_orbit_bin.f90
	gfortran -c -o ./obj/af_orbit_bin.o ../AstroFun/af_orbit_bin.f90
./obj/af_mag_accretion.o : ../AstroFun/af_mag_accretion.f90
	gfortran -c -o ./obj/af_mag_accretion.o ../AstroFun/af_mag_accretion.f90
./obj/af_RT_layer_an.o : ../AstroFun/af_RT_layer_an.f90
	gfortran -c -o ./obj/af_RT_layer_an.o ../AstroFun/af_RT_layer_an.f90


./obj/PhFun_eMaxvell.o : ../PhFun/PhFun_eMaxvell.f90
	gfortran -c -o ./obj/PhFun_eMaxvell.o ../PhFun/PhFun_eMaxvell.f90
./obj/PhFun_Compton.o : ../PhFun/PhFun_Compton.f90
	gfortran -c -o ./obj/PhFun_Compton.o ../PhFun/PhFun_Compton.f90
./obj/PhFun_ComptonBF_H.o : ../PhFun/PhFun_ComptonBF_H.f90
	gfortran -fopenmp -c -o ./obj/PhFun_ComptonBF_H.o ../PhFun/PhFun_ComptonBF_H.f90
./obj/PhFun_ComptonBF_cons_laws.o : ../PhFun/PhFun_ComptonBF_cons_laws.f90
	gfortran -c -o ./obj/PhFun_ComptonBF_cons_laws.o ../PhFun/PhFun_ComptonBF_cons_laws.f90
./obj/PhFun_ComptonBF_CrossSection.o : ../PhFun/PhFun_ComptonBF_CrossSection.f90
	gfortran -fopenmp -c -o ./obj/PhFun_ComptonBF_CrossSection.o ../PhFun/PhFun_ComptonBF_CrossSection.f90
./obj/PhFun_NormalWaves.o : ../PhFun/PhFun_NormalWaves.f90
	gfortran  -fopenmp -c -o ./obj/PhFun_NormalWaves.o ../PhFun/PhFun_NormalWaves.f90
./obj/PhFun_Landau_Level_width.o : ../PhFun/PhFun_Landau_Level_width.f90
	gfortran -c -o ./obj/PhFun_Landau_Level_width.o ../PhFun/PhFun_Landau_Level_width.f90
./obj/PhFun_Cyclotron_absorb.o : ../PhFun/PhFun_Cyclotron_absorb.f90
	gfortran -c -o ./obj/PhFun_Cyclotron_absorb.o ../PhFun/PhFun_Cyclotron_absorb.f90
./obj/PhFun_magnetic_absorb.o : ../PhFun/PhFun_magnetic_absorb.f90
	gfortran -c -o ./obj/PhFun_magnetic_absorb.o ../PhFun/PhFun_magnetic_absorb.f90
./obj/PhFun_ComptonBF_Stokes.o : ../PhFun/PhFun_ComptonBF_Stokes.f90
	gfortran -fopenmp -c -o ./obj/PhFun_ComptonBF_Stokes.o ../PhFun/PhFun_ComptonBF_Stokes.f90
./obj/PhFun_Gaunt_factors.o : ../PhFun/PhFun_Gaunt_factors.f90
	gfortran -c -o ./obj/PhFun_Gaunt_factors.o ../PhFun/PhFun_Gaunt_factors.f90
./obj/PhFun_PlasmaDispFun.o : ../PhFun/PhFun_PlasmaDispFun.f
	gfortran -c -o ./obj/PhFun_PlasmaDispFun.o ../PhFun/PhFun_PlasmaDispFun.f
./obj/PhFun_MagneticDipole.o : ../PhFun/PhFun_MagneticDipole.f90
	gfortran -c -o ./obj/PhFun_MagneticDipole.o ../PhFun/PhFun_MagneticDipole.f90
./obj/PhFun_GR_app.o : ../PhFun/PhFun_GR_app.f90
	gfortran -fopenmp -c -o ./obj/PhFun_GR_app.o ../PhFun/PhFun_GR_app.f90
./obj/PhFun_Polarization.o : ../PhFun/PhFun_Polarization.f90
	gfortran -c -o ./obj/PhFun_Polarization.o ../PhFun/PhFun_Polarization.f90
./obj/PhFun_dielec_tensor.o : ../PhFun/PhFun_dielec_tensor.f90
	gfortran -c -o ./obj/PhFun_dielec_tensor.o ../PhFun/PhFun_dielec_tensor.f90
./obj/PhFun_blackbody.o : ../PhFun/PhFun_blackbody.f90
	gfortran -c -o ./obj/PhFun_blackbody.o ../PhFun/PhFun_blackbody.f90
./obj/PhFun_SpR.o : ../PhFun/PhFun_SpR.f90
	gfortran -c -o ./obj/PhFun_SpR.o ../PhFun/PhFun_SpR.f90

#./obj/PhFun_Compton_Harding.o : ../PhFun/PhFun_Compton_Harding.f90
#	gfortran -c -o ./obj/PhFun_Compton_Harding.o ../PhFun/PhFun_Compton_Harding.f90

./obj/PhFun_Compton_Harding.o : ../PhFun/PhFun_Compton_Harding.f90
	gfortran `pkg-config --cflags fgsl` -c -o ./obj/PhFun_Compton_Harding.o ../PhFun/PhFun_Compton_Harding.f90 integral `pkg-config --libs fgsl`


./obj/mf_algebra_polinom.o : ../mf/mf_algebra_polinom.f90
	gfortran -c -o ./obj/mf_algebra_polinom.o ../mf/mf_algebra_polinom.f90
./obj/mf_geom_2d.o : ../mf/mf_geom_2d.f90
	gfortran -c -o ./obj/mf_geom_2d.o ../mf/mf_geom_2d.f90
./obj/mf_geom_3d.o : ../mf/mf_geom_3d.f90
	gfortran -c -o ./obj/mf_geom_3d.o ../mf/mf_geom_3d.f90
./obj/mf_tcrUL_VF.o : ../mf/mf_tcrUL_VF.f
	gfortran -c -o ./obj/mf_tcrUL_VF.o ../mf/mf_tcrUL_VF.f
./obj/mf_expmag_VF.o : ../mf/mf_expmag_VF.f
	gfortran -c -o ./obj/mf_expmag_VF.o ../mf/mf_expmag_VF.f
./obj/mf_bessel_VF.o : ../mf/mf_bessel_VF.f
	gfortran -c -o ./obj/mf_bessel_VF.o ../mf/mf_bessel_VF.f
./obj/mf_bessel.o : ../mf/mf_bessel.f90
	gfortran -c -o ./obj/mf_bessel.o ../mf/mf_bessel.f90
./obj/mf_c.o : ../mf/mf_c.f90
	gfortran -c -o ./obj/mf_c.o ../mf/mf_c.f90
./obj/mf_ermit_pol.o : ../mf/mf_ermit_pol.f90
	gfortran -c -o ./obj/mf_ermit_pol.o ../mf/mf_ermit_pol.f90
./obj/mf_erm_link.o : ../mf/mf_erm_link.f90
	gfortran -c -o ./obj/mf_erm_link.o ../mf/mf_erm_link.f90
./obj/mf_fact.o : ../mf/mf_fact.f90
	gfortran -c -o ./obj/mf_fact.o ../mf/mf_fact.f90
./obj/mf_fi_function.o : ../mf/mf_fi_function.f90
	gfortran -c -o ./obj/mf_fi_function.o ../mf/mf_fi_function.f90
./obj/mf_gamma-function.o : ../mf/mf_gamma-function.f90
	gfortran -c -o ./obj/mf_gamma-function.o ../mf/mf_gamma-function.f90
./obj/mf_gener_fun.o : ../mf/mf_gener_fun.f90
	gfortran -c -o ./obj/mf_gener_fun.o ../mf/mf_gener_fun.f90
./obj/mf_g.o : ../mf/mf_g.f90
	gfortran -c -o ./obj/mf_g.o ../mf/mf_g.f90
./obj/mf_integrals.o : ../mf/mf_integrals.f90
	gfortran -c -o ./obj/mf_integrals.o ../mf/mf_integrals.f90
./obj/mf_integral_2d.o : ../mf/mf_integral_2d.f90
	gfortran -c -o ./obj/mf_integral_2d.o ../mf/mf_integral_2d.f90
./obj/mf_int_gauss_100.o : ../mf/mf_int_gauss_100.f90
	gfortran -c -o ./obj/mf_int_gauss_100.o ../mf/mf_int_gauss_100.f90
./obj/mf_int_gauss.o : ../mf/mf_int_gauss.f90
	gfortran -c -o ./obj/mf_int_gauss.o ../mf/mf_int_gauss.f90
./obj/mf_intGauss.o : ../mf/mf_intGauss.FOR
	gfortran -c -o ./obj/mf_intGauss.o ../mf/mf_intGauss.FOR
./obj/mf_intGauss_p.o : ../mf/mf_intGauss_p.FOR
	gfortran -c -o ./obj/mf_intGauss_p.o ../mf/mf_intGauss_p.FOR
./obj/mf_int_G.o : ../mf/mf_int_G.FOR
	gfortran -c -o ./obj/mf_int_G.o ../mf/mf_int_G.FOR
./obj/mf_intSimps.o : ../mf/mf_intSimps.f90
	gfortran -c -o ./obj/mf_intSimps.o ../mf/mf_intSimps.f90
./obj/mf_intSimp_p.o : ../mf/mf_intSimp_p.FOR
	gfortran -c -o ./obj/mf_intSimp_p.o ../mf/mf_intSimp_p.FOR
./obj/mf_ksi.o : ../mf/mf_ksi.f90
	gfortran -c -o ./obj/mf_ksi.o ../mf/mf_ksi.f90
./obj/mf_laguer_exp.o : ../mf/mf_laguer_exp.f90
	gfortran -c -o ./obj/mf_laguer_exp.o ../mf/mf_laguer_exp.f90
./obj/mf_laguer.o : ../mf/mf_laguer.f90
	gfortran -c -o ./obj/mf_laguer.o ../mf/mf_laguer.f90
./obj/mf_laguer_gen.o : ../mf/mf_laguer_gen.f90
	gfortran -c -o ./obj/mf_laguer_gen.o ../mf/mf_laguer_gen.f90
./obj/mf_laguer_n.o : ../mf/mf_laguer_n.f90
	gfortran -c -o ./obj/mf_laguer_n.o ../mf/mf_laguer_n.f90
./obj/mf_matrix_matrix.o : ../mf/mf_matrix_matrix.f90
	gfortran -c -o ./obj/mf_matrix_matrix.o ../mf/mf_matrix_matrix.f90
./obj/mf_num_mart.o : ../mf/mf_num_mart.f90
	gfortran -c -o ./obj/mf_num_mart.o ../mf/mf_num_mart.f90
./obj/mf_numbers.o : ../mf/mf_numbers.f90
	gfortran -c -o ./obj/mf_numbers.o ../mf/mf_numbers.f90
./obj/mf_vector_sk_mult.o : ../mf/mf_vector_sk_mult.f90
	gfortran -c -o ./obj/mf_vector_sk_mult.o ../mf/mf_vector_sk_mult.f90
./obj/mf_coord_sys.o : ../mf/mf_coord_sys.f90
	gfortran -c -o ./obj/mf_coord_sys.o ../mf/mf_coord_sys.f90
./obj/mf_for_mass.o : ../mf/mf_for_mass.f90
	gfortran -c -o ./obj/mf_for_mass.o ../mf/mf_for_mass.f90
./obj/mf_random.o : ../mf/mf_random.f90
	gfortran -c -o ./obj/mf_random.o ../mf/mf_random.f90
./obj/mf_statistics.o : ../mf/mf_statistics.f90
	gfortran -c -o ./obj/mf_statistics.o ../mf/mf_statistics.f90
./obj/mf_equations.o : ../mf/mf_equations.f90
	gfortran -c -o ./obj/mf_equations.o ../mf/mf_equations.f90
./obj/mf_intExp.o : ../mf/mf_intExp.f90
	gfortran -c -o ./obj/mf_intExp.o ../mf/mf_intExp.f90
./obj/mf_interpol.o : ../mf/mf_interpol.f90
	gfortran -c -o ./obj/mf_interpol.o ../mf/mf_interpol.f90

