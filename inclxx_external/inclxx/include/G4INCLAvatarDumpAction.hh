/** \file G4INCLAvatarDumpAction.hh
 * \brief Alternative CascadeAction for dumping avatars
 *
 * \date 15th October 2014
 * \author Davide Mancusi
 */

#ifndef G4INCLAVATARDUMPACTION_HH
#define G4INCLAVATARDUMPACTION_HH 1

#include "G4INCLCascadeAction.hh"
#include <fstream>

namespace G4INCL {

  class AvatarDumpAction : public CascadeAction {

    public:
      AvatarDumpAction();
      virtual ~AvatarDumpAction();

      virtual void beforeCascadeUserAction(IPropagationModel *);
      virtual void afterAvatarUserAction(IAvatar *avatar, Nucleus *nucleus, FinalState *);
      virtual void afterCascadeUserAction(Nucleus *);

    private:
      std::ofstream *oFile;
      int eventCounter;

  };

}
#endif // G4INCLAVATARDUMPACTION_HH
