MAKEFLAGS += --no-builtin-rules --no-builtin-variables
FC := gfortran
RM := F:/MSYS2/usr/bin/rm.exe
OBJ_Suf := .o
FPP_Suf := .f90
Mod_SRS := basic_data$(FPP_Suf) solver_data$(FPP_Suf) mat_eqn_solve$(FPP_Suf) \
get_matrix$(FPP_Suf) integration$(FPP_Suf) element$(FPP_Suf) read$(FPP_Suf) \
solver_kernel$(FPP_Suf) output$(FPP_Suf)
Mod_OBJ := $(Mod_SRS:%$(FPP_Suf)=%$(OBJ_Suf))
Mod_MOD := $(Mod_SRS:%$(FPP_Suf)=%.mod)
HeatStress: HeatStress$(FPP_Suf) $(Mod_OBJ)
	$(FC) $^ -o $@
$(Mod_OBJ): %$(OBJ_Suf): %$(FPP_Suf)
	$(FC) -c $<
$(Mod_MOD): %.mod: %$(FPP_Suf)
	$(FC) -c $<
.PHONY: clean
clean:
	$(RM) *.mod *$(OBJ_Suf)