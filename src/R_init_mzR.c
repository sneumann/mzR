#include <R_ext/Rdynload.h>
#include "rnetCDF.h"

static const R_CMethodDef cMethods[] = {
    {"NetCDFStrError", (DL_FUNC) &NetCDFStrError, 3},
    {"NetCDFOpen", (DL_FUNC) &NetCDFOpen, 3},
    {"NetCDFClose", (DL_FUNC) &NetCDFClose, 2},
    {"NetCDFVarID", (DL_FUNC) &NetCDFVarID, 4},
    {"NetCDFVarLen", (DL_FUNC) &NetCDFVarLen, 4},
    {"NetCDFVarDouble", (DL_FUNC) &NetCDFVarDouble, 4},
    {"NetCDFVarInt", (DL_FUNC) &NetCDFVarInt, 4},
    {"NetCDFMSPoints", (DL_FUNC) &NetCDFMSPoints, 7},
    {NULL, NULL, 0}
};

void R_init_mzR(DllInfo * info)
{
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
