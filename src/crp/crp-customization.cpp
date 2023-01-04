#include "crp.h"

// TODO: add actual customization routines

namespace crp
{
void CRPAlgorithm::customize()
{
    this->reverse = g->reversed();
    this->params.customizer(g, overlay.get());
}
} // namespace crp
