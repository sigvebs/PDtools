#ifndef FORCES
#define FORCES

#include "force.h"

// Bond Based Peridynamics
#include "PdForces/pd_bondforce.h"
#include "PdForces/pd_dampenedbondforce.h"
#include "PdForces/pd_bondforcegaussian.h"
#include "PdForces/viscousdamper.h"
#include "PdForces/pd_pmb.h"
#include "PdForces/pd_pmb_linear_integrator.h"

// Ordinary State Based Peridynamics
#include "PdForces/pd_osp.h"
#include "PdForces/pd_lps.h"
#include "PdForces/lps_mc.h"
#include "PdForces/pd_lps_adrmc.h"
#include "PdForces/pd_lpsdampenedcontact.h"

// Non-Ordinary State Based Peridynamics
#include "PdForces/pd_nopd.h"

// Other
#include "PdForces/viscousdamper.h"
#include "PdForces/contactforce.h"
#include "DemForces/demforce.h"

#endif // FORCES

