#include <cctk.h>

extern "C" int FCKleinGordon_startup(void) {
  const char *banner = "FCKleinGordon_startup: Evolutions of a Klein gordon field over an "
                       "arbitrary background using a first order flux conservative formulation";
  CCTK_RegisterBanner(banner);

  return 0;
}