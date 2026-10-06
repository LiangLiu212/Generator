/** \file G4INCLCascadeAction.hh
 * \brief Class containing default actions to be performed at intermediate cascade steps
 *
 * \date 22nd October 2013
 * \author Davide Mancusi
 */

#ifndef G4INCLCASCADEACTION_HH
#define G4INCLCASCADEACTION_HH 1

#include "G4INCLIAvatar.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLIAvatar.hh"
#include "G4INCLConfig.hh"

namespace G4INCL {

  class CascadeAction {
    // class INCL must be a friend because it needs to call private methods
    friend class INCL;

    public:
    CascadeAction();
    virtual ~CascadeAction();

    virtual void beforeRunUserAction(Config const *) {}
    virtual void beforeCascadeUserAction(IPropagationModel *) {}
    virtual void beforePropagationUserAction(IPropagationModel *) {}
    virtual void beforeAvatarUserAction(IAvatar *, Nucleus *) {}
    virtual void afterAvatarUserAction(IAvatar *, Nucleus *, FinalState *) {}
    virtual void afterPropagationUserAction(IPropagationModel *, IAvatar *) {}
    virtual void afterCascadeUserAction(Nucleus *) {}
    virtual void afterRunUserAction() {}

    private:
    // These four methods should be private because the user must not be
    // allowed to override them
    void beforeRunAction(Config const *config);
    void beforeCascadeAction(IPropagationModel *);
    void beforePropagationAction(IPropagationModel *pm);
    void beforeAvatarAction(IAvatar *a, Nucleus *n);
    void afterAvatarAction(IAvatar *a, Nucleus *n, FinalState *fs);
    void afterPropagationAction(IPropagationModel *pm, IAvatar *avatar);
    void afterCascadeAction(Nucleus *);
    void afterRunAction();

    void beforeRunDefaultAction(Config const *config);
    void beforeCascadeDefaultAction(IPropagationModel *pm);
    void beforePropagationDefaultAction(IPropagationModel *pm);
    void beforeAvatarDefaultAction(IAvatar *a, Nucleus *n);
    void afterAvatarDefaultAction(IAvatar *a, Nucleus *n, FinalState *fs);
    void afterPropagationDefaultAction(IPropagationModel *pm, IAvatar *avatar);
    void afterCascadeDefaultAction(Nucleus *);
    void afterRunDefaultAction();

    private: // data members
    long stepCounter;
  };

}
#endif // G4INCLCASCADEACTION_HH
