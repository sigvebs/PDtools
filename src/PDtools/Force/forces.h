#ifndef FORCES
#define FORCES

#include "force.h"

// Bond Based Peridynamics
#include "PDtools/Force/PdForces/pd_bondforce.h"
#include "PDtools/Force/PdForces/pd_bondforcegaussian.h"
#include "PDtools/Force/PdForces/pd_dampenedbondforce.h"
#include "PDtools/Force/PdForces/pd_pmb.h"
#include "PDtools/Force/PdForces/pd_pmb_linear_integrator.h"
#include "PDtools/Force/PdForces/viscousdamper.h"

// Ordinary State Based Peridynamics
#include "PDtools/Force/PdForces/LPS/lps_mc.h"
#include "PDtools/Force/PdForces/LPS/pd_lps.h"
#include "PDtools/Force/PdForces/LPS/pd_lps2.h"
#include "PDtools/Force/PdForces/LPS/pd_lps_adr_strain.h"
#include "PDtools/Force/PdForces/LPS/pd_lps_adrmc.h"
#include "PDtools/Force/PdForces/LPS/pd_lps_crit_strain.h"
#include "PDtools/Force/PdForces/LPS/pd_lps_k.h"
#include "PDtools/Force/PdForces/LPS/pd_lpsdampenedcontact.h"
#include "PDtools/Force/PdForces/pd_osp.h"

#include "PDtools/Force/PdForces/LPS_porosity/lps_p_mc.h"
#include "PDtools/Force/PdForces/LPS_porosity/pd_lps_p.h"
#include "PDtools/Force/PdForces/LPS_porosity/pd_lps_p_adrmc.h"
#include "PDtools/Force/PdForces/LPS_porosity/pd_lpsdampenedcontact_p.h"

// Non-Ordinary State Based Peridynamics
#include "PDtools/Force/PdForces/LPSS/pd_lpss.h"
#include "PDtools/Force/PdForces/LPSS/pd_lpss_g.h"
#include "PDtools/Force/PdForces/LPSS/pd_lpss_opt.h"
#include "PDtools/Force/PdForces/pd_nopd.h"

// Other
#include "PDtools/Force/DemForces/demforce.h"
#include "PDtools/Force/PdForces/contactforce.h"
#include "PDtools/Force/PdForces/viscousdamper.h"

#endif // FORCES
