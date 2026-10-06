/** \file G4INCLXXInterfaceMessenger.hh
 * \brief Messenger class for the Geant4 INCL++ interface.
 *
 * \date 26th April 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLXXInterfaceMessenger_hh
#define G4INCLXXInterfaceMessenger_hh

#include "G4INCLXXInterfaceStore.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4String.hh"

class G4INCLXXInterfaceStore;

class G4INCLXXInterfaceMessenger : public G4UImessenger
{

  public:
    G4INCLXXInterfaceMessenger (G4INCLXXInterfaceStore *anInterfaceStore);
    ~G4INCLXXInterfaceMessenger ();
    void SetNewValue (G4UIcommand *command, G4String newValues);
    //    Identifies the command which has been invoked by the user, extracts the
    //    parameters associated with that command (held in newValues, and uses
    //    these values with the appropriate member function of G4INCLXXInterfaceStore.
    //
  private:
    static const G4String theUIDirectory;
    G4INCLXXInterfaceStore *theINCLXXInterfaceStore;
    G4UIdirectory *theINCLXXDirectory;
    G4UIcmdWithAString *accurateNucleusCmd;
    G4UIcmdWithAnInteger *maxClusterMassCmd;
    G4UIcmdWithADoubleAndUnit *cascadeMinEnergyPerNucleonCmd;
    G4UIcmdWithAString *inclPhysicsCmd;
    G4UIcommand *useAblaCmd;
};

#endif

