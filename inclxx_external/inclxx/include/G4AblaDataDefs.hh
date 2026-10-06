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

// Data structures needed by ABLA evaporation code.

#pragma once

#include "G4INCLGeant4Compat.hh"

#include <cmath>
#include <vector>

constexpr G4int nrows = 180;
constexpr G4int zcols = 122;

constexpr G4int lpcols = 13;
constexpr G4int lprows = 154;

constexpr G4int nrowsbeta = 251;
constexpr G4int zcolsbeta = 137;

constexpr G4int indexpart = 300;

constexpr G4double fmp = 938.27231, fmn = 939.56563, fml = 1115.683;

constexpr G4double PI = 3.14159265358979323846;
constexpr G4double AMU = 931.49410242; // MeV/c²
constexpr G4double C = 29.9792458;     // cm/ns

class G4Mexp
{

  public:
    G4Mexp(){};

    virtual ~G4Mexp() = default;

    G4double massexp[lprows][lpcols] = { { 0. } };
    G4double bind[lprows][lpcols] = { { 0. } };
    G4int mexpiop[lprows][lpcols] = { { 0 } };
};

class G4Ec2sub
{
  public:
    G4Ec2sub(){};

    virtual ~G4Ec2sub() = default;

    G4double ecnz[nrows][zcols] = { { 0. } };
};

class G4Ald
{
  public:
    G4Ald()
        : av(0.0)
        , as(0.0)
        , ak(0.0)
        , optafan(0.0){};

    virtual ~G4Ald() = default;

    G4double av, as, ak, optafan = 0.;
};

/**
 * Shell corrections and deformations.
 */

class G4Ecld
{

  public:
    G4Ecld(){};

    virtual ~G4Ecld() = default;

    /**
     * Ground state shell correction frldm for a spherical ground state.
     */
    G4double ecgnz[nrows][zcols] = { { 0. } };

    /**
     * Shell correction for the saddle point.
     */
    G4double ecfnz[nrows][zcols] = { { 0. } };

    /**
     * Difference between deformed ground state and ldm value.
     */
    G4double vgsld[nrows][zcols] = { { 0. } };

    /**
     * Alpha ground state deformation (this is not beta2!)
     * beta2 = std::sqrt(5/(4pi)) * alpha
     */
    G4double alpha[nrows][zcols] = { { 0. } };

    /**
     * RMS function for lcp emission barriers
     */
    G4double rms[nrows][zcols] = { { 0. } };

    /**
     * Beta2 deformations
     */
    G4double beta2[nrowsbeta][zcolsbeta] = { { 0. } };

    /**
     * Beta4 deformations
     */
    G4double beta4[nrowsbeta][zcolsbeta] = { { 0. } };
};

class G4Fiss
{
    /**
     * Options and parameters for fission channel.
     */

  public:
    G4Fiss()
        : bet(0.0)
        , bethyp(0.0)
        , ifis(0.0)
        , ucr(0.0)
        , dcr(0.0)
        , optshp(0)
        , optxfis(0)
        , optct(0)
        , optcol(0)
        , at(0)
        , zt(0){};

    virtual ~G4Fiss() = default;

    G4double bet, bethyp, ifis, ucr, dcr;
    G4int optshp, optxfis, optct, optcol, at, zt;
};

/**
 * Fission and emission particle barriers.
 */

class G4Fb
{

  public:
    G4Fb()
        : h2barfact(0.0)
        , h3barfact(0.0)
        , he3barfact(0.0)
        , he4barfact(0.0)
        , he6barfact(0.0){};

    virtual ~G4Fb() = default;

    float h2barfact, h3barfact, he3barfact, he4barfact, he6barfact;
    G4double efa[nrows][zcols] = { { 0. } };
    G4double mmfb[nrows][zcols] = { { 0. } };
};

/**
 * Options
 */

class G4Opt
{

  public:
    G4Opt()
        : optemd(0)
        , optcha(0)
        , optshpimf(0)
        , optimfallowed(0)
        , nblan0(0){};

    virtual ~G4Opt() = default;

    G4int optemd, optcha, optshpimf, optimfallowed, nblan0;
};

class G4AblaOutput
{
  public:
    G4AblaOutput() { clear(); };

    virtual ~G4AblaOutput() = default;

    /**
     * Clear and initialize all variables and arrays.
     */
    void clear()
    {
        ntrack = 0;
        kfis = 0;
        estfis = 0;
        izfis = 0;
        iafis = 0;
        esci = 0;
        iasci.clear();
        izsci.clear();
        fissmode = -1;
        itypcasc.clear();
        avv.clear();
        zvv.clear();
        svv.clear();
        jvv.clear();
        enerj.clear();
        pxlab.clear();
        pylab.clear();
        pzlab.clear();
    }

    /**
     * Fission 1/0=Y/N.
     */
    G4int kfis;

    /**
     * Excit energy at fis.
     */
    G4double estfis;

    /**
     * Z of fiss nucleus.
     */
    G4int izfis;

    /**
     * A of fiss nucleus.
     */
    G4int iafis;

    /**
     * Excit energy at fis.
     */
    G4double esci;

    /**
     * Z of fiss nucleus.
     */
    std::vector<G4int> izsci;

    /**
     * A of fiss nucleus.
     */
    std::vector<G4int> iasci;

    /**
     * Number of particles.
     */
    G4int ntrack;

    /**
     * Fission mode.
     */
    G4int fissmode;

    /**
     * emitted in cascade (0) or evaporation (1).
     */
    std::vector<G4int> itypcasc;

    /**
     * A (-1 for pions).
     */
    std::vector<G4int> avv;

    /**
     * Z number.
     */
    std::vector<G4int> zvv;

    /**
     * S (-1 for lambda_0).
     */
    std::vector<G4int> svv;

    /**
     * J angular momemtum.
     */
    std::vector<G4double> jvv;

    /**
     * Kinetic energy.
     */
    std::vector<G4double> enerj;

    /**
     * Momentum.
     */
    std::vector<G4double> plab;
    std::vector<G4double> pxlab;
    std::vector<G4double> pylab;
    std::vector<G4double> pzlab;

    /**
     * Theta angle.
     */
    std::vector<G4double> tetlab;

    /**
     * Phi angle.
     */
    std::vector<G4double> philab;
};
