#FC = mpif90 -fc=/opt/pgi/osx86-64/2013/bin/pgfortran -C
#FC = mpif90 -f90=/srv/pgi/linux86-64/9.0/bin/pgfortran -fast -O4 -fastsse -Mstandard -Mfprelaxed -Mipa=fast,inline
FC = mpif90 -O3 -fdefault-double-8 -fdefault-real-8 -fno-align-commons
FO=libmartin.a

MODS = 
OBJS = mar_para.o mar_data.o calc_mod.o model.o  mat_funcs.o pinput.o files.o calc.o spline.o sys_resp.o for_funcs.o ez_tx.o inline_f_fht.o hr_f_fht.o ex_f_fht.o mded_f_fht.o 
                  
$(FO): ${MODS} ${OBJS}
	ar urv libmartin.a ${MODS} ${OBJS}
	cp *.a ../obj/
	cp *.mod ../

mar_para.o: mar_para.f90; $(FC) -c mar_para.f90
mar_data.o: mar_data.f90; $(FC) -c mar_data.f90
calc_mod.o: calc_mod.f90; $(FC) -c calc_mod.f90
model.o: model.f90; $(FC) -c model.f90
mat_funcs.o: mat_funcs.f90; $(FC) -c mat_funcs.f90
pinput.o: pinput.f90; $(FC) -c pinput.f90
files.o: files.f90; $(FC) -c files.f90
calc.o: calc.f90; $(FC) -c calc.f90
spline.o: spline.f90; $(FC) -c spline.f90
sys_resp.o: sys_resp.f90; $(FC) -c sys_resp.f90
for_funcs.o: for_funcs.f90; $(FC) -c for_funcs.f90
ez_tx.o: ez_tx.f90; $(FC) -c ez_tx.f90
inline_f_fht.o: inline_f_fht.f90; $(FC) -c inline_f_fht.f90
hr_f_fht.o: hr_f_fht.f90; $(FC) -c hr_f_fht.f90
ex_f_fht.o: ex_f_fht.f90; $(FC) -c ex_f_fht.f90
mded_f_fht.o: mded_f_fht.f90; $(FC) -c mded_f_fht.f90
clean:	
	/bin/rm -f *.o *.mod *.l core *.a
