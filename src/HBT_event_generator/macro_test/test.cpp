#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_version.h>

//#if GSL_MAJOR_VERSION < 2
//	#include "fit_v1.h"
//#else
//	#include "fit_v2.h"
#include "fit_v1.h"
#include "fit_v2.h"

int main (void)
{
	fit_driver();

	return 0;
}