
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "r-cv-svd-gabriel.h"
#include "r-cv-svd-wold.h"
#include "r-impute-svd.h"

static R_CallMethodDef callMethods[] = {
    { "R_impute_svd",     (DL_FUNC) &R_impute_svd,     4 },
    { "R_cv_svd_gabriel", (DL_FUNC) &R_cv_svd_gabriel, 6 },
    { "R_cv_svd_wold",    (DL_FUNC) &R_cv_svd_wold,    6 },
    { NULL,               NULL,                        0 }
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
