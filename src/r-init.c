
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "r-bcv-svd-gabriel.h"
#include "r-bcv-svd-wold.h"
#include "r-svd-impute.h"

static R_CallMethodDef callMethods[] = {
    { "R_svd_impute",      (DL_FUNC) &R_svd_impute,      4 },
    { "R_bcv_svd_gabriel", (DL_FUNC) &R_bcv_svd_gabriel, 6 },
    { "R_bcv_svd_wold",    (DL_FUNC) &R_bcv_svd_gabriel, 6 },
    { NULL,                NULL,                         0 }
};

void
R_init_bcv (DllInfo *info)
{
    R_registerRoutines (info, NULL, callMethods, NULL, NULL);
}

void
R_unload_bcv (DllInfo *info)
{
    
}
