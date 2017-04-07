#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _Ca_HVA_reg();
extern void _Ca_LVAst_reg();
extern void _Ih_reg();
extern void _ca_h_reg();
extern void _ca_r_reg();
extern void _cad_reg();
extern void _distr_reg();
extern void _kadist_reg();
extern void _kaprox_reg();
extern void _kca_reg();
extern void _km_reg();
extern void _kv_reg();
extern void _multiclamp_reg();
extern void _na_reg();
extern void _zoidsyn_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," Ca_HVA.mod");
fprintf(stderr," Ca_LVAst.mod");
fprintf(stderr," Ih.mod");
fprintf(stderr," ca_h.mod");
fprintf(stderr," ca_r.mod");
fprintf(stderr," cad.mod");
fprintf(stderr," distr.mod");
fprintf(stderr," kadist.mod");
fprintf(stderr," kaprox.mod");
fprintf(stderr," kca.mod");
fprintf(stderr," km.mod");
fprintf(stderr," kv.mod");
fprintf(stderr," multiclamp.mod");
fprintf(stderr," na.mod");
fprintf(stderr," zoidsyn.mod");
fprintf(stderr, "\n");
    }
_Ca_HVA_reg();
_Ca_LVAst_reg();
_Ih_reg();
_ca_h_reg();
_ca_r_reg();
_cad_reg();
_distr_reg();
_kadist_reg();
_kaprox_reg();
_kca_reg();
_km_reg();
_kv_reg();
_multiclamp_reg();
_na_reg();
_zoidsyn_reg();
}
