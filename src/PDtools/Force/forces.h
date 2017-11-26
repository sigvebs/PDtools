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
#include "PdForces/LPS/pd_lps.h"
#include "PdForces/LPS/pd_lps2.h"
#include "PdForces/LPS/lps_mc.h"
#include "PdForces/LPS/pd_lps_adrmc.h"
#include "PdForces/LPS/pd_lpsdampenedcontact.h"
#include "PdForces/LPS/pd_lps_crit_strain.h"
#include "PdForces/LPS/pd_lps_adr_strain.h"
#include "PdForces/LPS/pd_lps_k.h"

#include "PdForces/LPS_porosity/pd_lps_p.h"
#include "PdForces/LPS_porosity/lps_p_mc.h"
#include "PdForces/LPS_porosity/pd_lps_p_adrmc.h"
#include "PdForces/LPS_porosity/pd_lpsdampenedcontact_p.h"

// Non-Ordinary State Based Peridynamics
#include "PdForces/pd_nopd.h"
#include "PdForces/LPSS/pd_lpss.h"
#include "PdForces/LPSS/pd_lpss_g.h"
#include "PdForces/LPSS/pd_lpss_opt.h"

// Other
#include "PdForces/viscousdamper.h"
#include "PdForces/contactforce.h"
#include "DemForces/demforce.h"

#endif // FORCES

