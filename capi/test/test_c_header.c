/* Compiled as C: miradac.h must be a valid C header. */
#include <string.h>
#include "miradac.h"

int main(void) {
    if (mdac_abi_version() != MDAC_ABI_VERSION) return 1;
    if (strcmp(mdac_version(), MDAC_EXPECTED_VERSION) != 0) return 2;
    return 0;
}
