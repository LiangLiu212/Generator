/******************************************************************************
 *             ABLA++ de-excitation model Copyright (C) 2018-2025             *
 *           J.L. Rodriguez-Sánchez, J.-C. David, and A. Kelic-Heil           *
 *                                                                            *
 * This software is distributed under the terms of the GNU General Public     *
 * License (GPL) version 3, which can be found in the LICENSE file.           *
 *                                                                            *
 * This software is provided "as is", without warranty of any kind, express   *
 * or implied, including but not limited to the warranties of merchantability *
 * and fitness for a particular purpose. See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with this software. If not, see <http://www.gnu.org/licenses/>.            *
 ******************************************************************************/

#pragma once

#include "G4INCLConfig.hh"
#include "G4INCLGeant4Compat.hh"

/**
 * An interface to data used by ABLA. This interface allows
 * us to abstract the actual source of data. Currently the data is
 * read from datafiles by using class G4AblaDataFile.  @see
 * G4AblaDataFile
 */

class G4AblaVirtualData
{
  protected:
    /**
     * Constructor, destructor
     */
    G4AblaVirtualData(G4INCL::Config*);

    virtual ~G4AblaVirtualData() = default;

  public:
    /**
     * Set the value of Alpha.
     */
    G4bool setAlpha(G4int N, G4int Z, G4double value);

    /**
     * Set the value of Ecnz.
     */
    G4bool setEcnz(G4int N, G4int Z, G4double value);

    /**
     * Set the value of Vgsld.
     */
    G4bool setVgsld(G4int N, G4int Z, G4double value);

    /**
     * Set the value of RMS.
     */
    G4bool setRms(G4int N, G4int Z, G4double value);

    /**
     * Set the value of experimental masses.
     */
    G4bool setMexp(G4int N, G4int Z, G4double value);

    /**
     * Set the value of experimental masses ID.
     */
    G4bool setMexpID(G4int N, G4int Z, G4int value);

    /**
     * Set the value of beta2 deformation.
     */
    G4bool setBeta2(G4int N, G4int Z, G4double value);

    /**
     * Set the value of beta4 deformation.
     */
    G4bool setBeta4(G4int N, G4int Z, G4double value);

    /**
     * Set the value of experimental Sn.
     */
    G4bool setSnexp(G4int N, G4int Z, G4double value);

    /**
     * Set the value of experimental Sp.
     */
    G4bool setSpexp(G4int N, G4int Z, G4double value);

    /**
     * Set the fission barrier from mmFRLDM.
     */
    G4bool setmmFb(G4int N, G4int Z, G4double value);

    /**
     * Get the value of Alpha.
     */
    G4double getAlpha(G4int N, G4int Z);

    /**
     * Get the value of Ecnz.
     */
    G4double getEcnz(G4int N, G4int Z);

    /**
     * Get the value of Vgsld.
     */
    G4double getVgsld(G4int N, G4int Z);

    /**
     * Get the value of RMS.
     */
    G4double getRms(G4int N, G4int Z);

    /**
     * Get the value of experimental masses.
     */
    G4double getMexp(G4int N, G4int Z);

    /**
     * Get the value of experimental masses ID.
     */
    G4int getMexpID(G4int N, G4int Z);

    /**
     * Get the value of beta2 deformation.
     */
    G4double getBeta2(G4int N, G4int Z);

    /**
     * Get the value of beta4 deformation.
     */
    G4double getBeta4(G4int N, G4int Z);

    /**
     * Get the value of experimental Sn.
     */
    G4double getSnexp(G4int N, G4int Z);

    /**
     * Get the value of experimental Sp.
     */
    G4double getSpexp(G4int N, G4int Z);

    /**
     * Get the fission barrier from mmFRLDM.
     */
    G4double getmmFb(G4int N, G4int Z);

    virtual G4bool readData() = 0;

  private:
    static const G4int sRows = 180;
    static const G4int sCols = 122;

    static const G4int betaRows = sRows + sCols;
    static const G4int betaCols = 137;

    G4double alpha[sRows][sCols];
    G4double ecnz[sRows][sCols];
    G4double vgsld[sRows][sCols];
    G4double rms[sRows][sCols];
    G4double mexp[sRows][sCols];
    G4int mexpid[sRows][sCols];
    G4double beta2[betaRows][betaCols];
    G4double beta4[betaRows][betaCols];
    G4double snexp[sRows][sCols];
    G4double spexp[sRows][sCols];
    G4double mmfb[sRows][sCols];
};
