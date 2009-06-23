
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "driver.h"

static R_CallMethodDef callMethods[] = {
    { "driver_svd", (DL_FUNC) &driver_svd, 6 },
    { NULL        , NULL                 , 0 }
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
