#include "Actions/action.h"
#include "Math/lattice.h"

//! Links an action object to a GluonField object
/*!
  \param lattice The gluonfield to link to the action instance.
*/
void Action::initAction(GluonField* lattice){
    m_lat = lattice;
}
