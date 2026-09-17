/*
 * =====================================================================================
 *
 *       Filename:  twobody_currents_sf.cxx
 *
 *    Description:  Native C++ port of N. Rocco's Fortran module dirac_matrices_intf
 *                  (ACHILLES, src/Achilles/fortran/currents_intf.f90 @ e02d266):
 *                  one-body current and the pion-in-flight + seagull + pion-pole +
 *                  Delta two-body currents (direct minus exchange) between a struck
 *                  nucleon and a spectator that stays in the Fermi sea.
 *
 *                  Reference: Lovato, Rocco, Steinberg, arXiv:2312.12545
 *
 *                  Differences with respect to the Fortran, none of which changes
 *                  a number when q is along z and C4V = C5V = 0:
 *                   - the sums over the Delta Lorentz indices are done
 *                     analytically (the pion vertex is proportional to the unit
 *                     matrix, the gamma-N-Delta vertex has two tensor structures),
 *                     operators sharing an isospin coefficient are added before
 *                     the spin matrix elements are taken, and products with Dirac
 *                     matrices use their sparsity: ~10 times faster;
 *                   - the C4V/C5V terms multiply the unit matrix, see
 *                     ModelParams::c4c5_broadcast;
 *                   - the pion-pole denominator uses the full q^2 instead of
 *                     q0^2 - qz^2;
 *                   - the in-medium Delta potential table (rho_0p5.dat) is not
 *                     read: the Fortran looks it up but never uses the value.
 *
 *         Author:  Liang Liu (L. Liu), liangliu@fnal.gov
 *		    Fermi National Accelerator Laboratory
 *  Collaboration:  GENIE
 *
 *  \cpright  Copyright (c) 2003-2026, The GENIE Collaboration
 *            For the full text of the license visit http://copyright.genie-mc.org
 *            or see $GENIE/LICENSE
 *
 * =====================================================================================
 */
#include <complex>
#include <array>
#include <cmath>

#include "Physics/HadronTensors/twobody_currents_sf.h"
#include "Physics/HadronTensors/TensorUtil.h"

namespace genie {
  namespace twobody_currents_sf {

    namespace {

      typedef std::complex<double> cplx;
      using genie::TensorUtil::Matrix2cd;
      using genie::TensorUtil::Matrix4cd;
      using genie::TensorUtil::Vector2cd;
      using genie::TensorUtil::Vector4cd;

      const cplx kCZero{0.0, 0.0};
      const cplx kCOne {1.0, 0.0};
      const cplx kCI   {0.0, 1.0};
      const double kMetric[4] = {1.0, -1.0, -1.0, -1.0};

      //---------------------------------------------------------------------------
      // Constant Dirac algebra (Fortran: dirac_matrices_in)
      //---------------------------------------------------------------------------
      struct DiracAlgebra {
        Matrix2cd sig[3];
        Matrix2cd id2;
        Matrix4cd id4;
        Matrix4cd ones4;            // every element = 1, see c4c5_broadcast
        Matrix4cd gam[5];           // gamma^0..gamma^3, gamma^5
        Matrix4cd gam5gam[4];       // gamma^5 gamma^mu
        Matrix4cd gamgam5[4];       // gamma^mu gamma^5
        Matrix4cd sigmunu[4][4];    // sigma^{mu nu} = i/2 [gamma^mu, gamma^nu]

        DiracAlgebra()
        {
          id2.setZero();
          id2(0,0) = kCOne;  id2(1,1) = kCOne;

          for(int i = 0; i < 3; ++i) sig[i].setZero();
          sig[0](0,1) =  kCOne;  sig[0](1,0) =  kCOne;
          sig[1](0,1) = -kCI;    sig[1](1,0) =  kCI;
          sig[2](0,0) =  kCOne;  sig[2](1,1) = -kCOne;

          id4.setZero();
          id4.block<2,2>(0,0) = id2;
          id4.block<2,2>(2,2) = id2;

          for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j) ones4(i,j) = kCOne;

          for(int k = 0; k < 5; ++k) gam[k].setZero();
          gam[0].block<2,2>(0,0) =  id2;
          gam[0].block<2,2>(2,2) = -id2;
          for(int mu = 1; mu <= 3; ++mu)
          {
            gam[mu].block<2,2>(0,2) =  sig[mu-1];
            gam[mu].block<2,2>(2,0) = -sig[mu-1];
          }
          gam[4].block<2,2>(0,2) = id2;
          gam[4].block<2,2>(2,0) = id2;

          for(int mu = 0; mu < 4; ++mu)
          {
            gam5gam[mu] = gam[4] * gam[mu];
            gamgam5[mu] = gam[mu] * gam[4];
            for(int nu = 0; nu < 4; ++nu)
            {
              sigmunu[mu][nu] = (kCI * 0.5) * (gam[mu] * gam[nu] - gam[nu] * gam[mu]);
            }
          }
        }
      };

      const DiracAlgebra & DA()
      {
        static const DiracAlgebra d;
        return d;
      }

      //---------------------------------------------------------------------------
      // Small 4-vector helpers, index 0..3 = (t,x,y,z)
      //---------------------------------------------------------------------------
      typedef std::array<double, 4> V4;

      double Dot(const V4 & a, const V4 & b)
      {
        return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
      }

      V4 Add(const V4 & a, const V4 & b)
      {
        return V4{{a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3]}};
      }

      V4 Sub(const V4 & a, const V4 & b)
      {
        return V4{{a[0]-b[0], a[1]-b[1], a[2]-b[2], a[3]-b[3]}};
      }

      Matrix4cd Slash(const V4 & a)
      {
        Matrix4cd m;
        for(int i = 0; i < 4; ++i) m += (kMetric[i] * a[i]) * DA().gam[i];
        return m;
      }

      // ubar * M * u   (Fortran sum(a*b): no implicit conjugation)
      cplx Sandwich(const Vector4cd & ubar, const Matrix4cd & m, const Vector4cd & u)
      {
        Vector4cd tmp = m * u;
        cplx val = kCZero;
        for(int a = 0; a < 4; ++a) val += ubar(a) * tmp(a);
        return val;
      }

      //---------------------------------------------------------------------------
      // Kinematics, spinors and isospinors of one phase-space point
      // (Fortran: current_init + define_spinors)
      //---------------------------------------------------------------------------
      struct Kinematics {
        V4 p1, p2, pp1, pp2, q;
        V4 k1, k2, k1e, k2e;
        Matrix4cd q_sl;
        Matrix4cd Pi_k1, Pi_k2, Pi_k1e, Pi_k2e;
        std::array<Vector4cd, 2> up1, up2;
        std::array<Vector4cd, 2> ubarpp1, ubarpp2;
        Vector2cd t1, t1p;
        int  ax;
        bool nc;
      };

      void BuildSpinors(const V4 & p, double xmn,
          std::array<Vector4cd, 2> & u, std::array<Vector4cd, 2> & ubar)
      {
        const DiracAlgebra & d = DA();
        Matrix2cd sigp = Matrix2cd::Zero();
        for(int i = 0; i < 3; ++i) sigp += d.sig[i] * p[i+1];

        const double c     = std::sqrt(p[0] + xmn);
        const double denom = p[0] + xmn;

        Vector2cd chi[2];
        chi[0] << kCOne,  kCZero;
        chi[1] << kCZero, kCOne;

        for(int s = 0; s < 2; ++s)
        {
          u[s].setZero();
          ubar[s].setZero();

          Vector2cd lower = sigp * chi[s];
          genie::TensorUtil::RowVector2cd row = chi[s].transpose() * sigp;
          for(int a = 0; a < 2; ++a)
          {
            u[s](a)      = c * chi[s](a);
            u[s](a+2)    = c * lower(a) / denom;
            ubar[s](a)   = c * chi[s](a);
            ubar[s](a+2) = -c * row(a) / denom;
          }
        }
      }

      Vector2cd Isospinor(int pdg)
      {
        Vector2cd t;
        if(pdg == 2212) t << kCOne,  kCZero;
        else            t << kCZero, kCOne;
        return t;
      }

      void InitKinematics(const ModelParams & par,
          const double p1_in[4], const double pp1_in[4],
          const double p2_in[4], const double q_in[4],
          int pdg_in, int pdg_out, bool has_axial, bool nc,
          Kinematics & k)
      {
        const DiracAlgebra & d = DA();

        k.nc = nc;
        k.ax = has_axial ? 1 : 0;
        k.t1  = Isospinor(pdg_in);
        k.t1p = Isospinor(pdg_out);

        for(int i = 0; i < 4; ++i)
        {
          k.p1[i]  = p1_in[i];
          k.pp1[i] = pp1_in[i];
          k.p2[i]  = p2_in[i];
          k.pp2[i] = p2_in[i];
          k.q[i]   = q_in[i];
        }

        // de Forest: struck nucleon on shell, energy transfer shifted accordingly
        const double w = k.q[0];
        k.q[0]  = w + k.p1[0];
        k.p1[0] = std::sqrt(k.p1[1]*k.p1[1] + k.p1[2]*k.p1[2] + k.p1[3]*k.p1[3]
            + par.xmn*par.xmn);
        k.q[0]  = k.q[0] - k.p1[0];

        k.k1  = Sub(k.pp1, k.p1);
        k.k2  = Sub(k.q,   k.k1);
        k.k1e = Sub(k.pp2, k.p1);
        k.k2e = Sub(k.pp1, k.p2);

        k.q_sl = Slash(k.q);

        const double xmpi2 = par.xmpi * par.xmpi;
        k.Pi_k1  = (d.gam[4] * Slash(k.k1))  / (Dot(k.k1,  k.k1)  - xmpi2);
        k.Pi_k2  = (d.gam[4] * Slash(k.k2))  / (Dot(k.k2,  k.k2)  - xmpi2);
        k.Pi_k1e = (d.gam[4] * Slash(k.k1e)) / (Dot(k.k1e, k.k1e) - xmpi2);
        k.Pi_k2e = (d.gam[4] * Slash(k.k2e)) / (Dot(k.k2e, k.k2e) - xmpi2);

        std::array<Vector4cd, 2> unused_u, unused_ubar;
        BuildSpinors(k.p1,  par.xmn, k.up1,    unused_ubar);
        BuildSpinors(k.p2,  par.xmn, k.up2,    unused_ubar);
        BuildSpinors(k.pp1, par.xmn, unused_u, k.ubarpp1);
        BuildSpinors(k.pp2, par.xmn, unused_u, k.ubarpp2);
      }

      //---------------------------------------------------------------------------
      // Isospin matrix elements (Fortran: me, iden, Ivz, Ivplus, IDelta*)
      //---------------------------------------------------------------------------
      cplx Me(int i, const Vector2cd & it, const Vector2cd & itp)
      {
        Vector2cd tmp = DA().sig[i] * it;
        return itp(0) * tmp(0) + itp(1) * tmp(1);
      }

      cplx Iden(const Vector2cd & it, const Vector2cd & itp)
      {
        return itp(0) * it(0) + itp(1) * it(1);
      }

      cplx Ivz(const Vector2cd & it1, const Vector2cd & it2,
          const Vector2cd & itp1, const Vector2cd & itp2)
      {
        return kCI * (Me(0,it1,itp1) * Me(1,it2,itp2) - Me(1,it1,itp1) * Me(0,it2,itp2));
      }

      cplx Ivplus(const Vector2cd & it1, const Vector2cd & it2,
          const Vector2cd & itp1, const Vector2cd & itp2)
      {
        return kCI * ( (Me(1,it1,itp1) * Me(2,it2,itp2) - Me(2,it1,itp1) * Me(1,it2,itp2))
            + kCI * (Me(2,it1,itp1) * Me(0,it2,itp2) - Me(0,it1,itp1) * Me(2,it2,itp2)) );
      }

      // Isospin coefficients of the four Delta diagrams and of the pion currents
      struct IsospinCoeffs { cplx a, b, c, d, pi; };

      IsospinCoeffs Isospin(bool charged, const Vector2cd & it1, const Vector2cd & it2,
          const Vector2cd & itp1, const Vector2cd & itp2)
      {
        IsospinCoeffs r;
        cplx iv, c2, c1;
        if(!charged)
        {
          iv = Ivz(it1, it2, itp1, itp2);
          c2 = Me(2, it2, itp2) * Iden(it1, itp1);
          c1 = Me(2, it1, itp1) * Iden(it2, itp2);
        }
        else
        {
          iv = Ivplus(it1, it2, itp1, itp2);
          c2 = (Me(0,it2,itp2) + kCI * Me(1,it2,itp2)) * Iden(it1, itp1);
          c1 = (Me(0,it1,itp1) + kCI * Me(1,it1,itp1)) * Iden(it2, itp2);
        }
        r.a  = 2.*c2/3. - iv/3.;
        r.b  = 2.*c2/3. + iv/3.;
        r.c  = 2.*c1/3. + iv/3.;
        r.d  = 2.*c1/3. - iv/3.;
        r.pi = iv;
        return r;
      }

      //---------------------------------------------------------------------------
      // One-body current operator (Fortran: det_J1)
      //---------------------------------------------------------------------------
      void DetJ1(const ModelParams & par, const FormFactors & ff, const Kinematics & k,
          Matrix4cd J_1b[4])
      {
        const DiracAlgebra & d = DA();
        for(int mu = 0; mu < 4; ++mu)
        {
          J_1b[mu].setZero();
          for(int nu = 0; nu < 4; ++nu)
          {
            J_1b[mu] += (kCI * ff.f2 * kMetric[nu] * k.q[nu] / 2.0 / par.xmn) * d.sigmunu[mu][nu];
          }
          J_1b[mu] += ff.f1 * d.gam[mu];
          J_1b[mu] += ff.fa * d.gamgam5[mu];
          J_1b[mu] += (ff.fap * k.q[mu] / par.xmn) * d.gam[4];
        }
      }

      //---------------------------------------------------------------------------
      // Delta decay width (Fortran: delta_se)
      //---------------------------------------------------------------------------
      double DeltaWidth(const ModelParams & par, double pd2)
      {
        const double thr = par.xmpi + par.xmn;
        if(pd2 < thr*thr) return 0.0;
        const double dif   = par.xmn - par.xmpi;
        const double kpi   = std::sqrt(0.25 / pd2 * (pd2 - thr*thr) * (pd2 - dif*dif));
        const double eknuc = std::sqrt(kpi*kpi + par.xmn*par.xmn);
        return std::pow(4.0 * par.fpind, 2) / 12.0 / M_PI / (par.xmpi*par.xmpi)
          * kpi*kpi*kpi / std::sqrt(pd2) * (par.xmn + eknuc);
      }

      //---------------------------------------------------------------------------
      // The Dirac matrices (and gamma^mu gamma^5, gamma^5 gamma^mu) have a single
      // non-zero element per row and column, so that multiplying by them is a
      // permutation of columns (rows) with a phase
      //---------------------------------------------------------------------------
      struct SparseGamma {
        int  col_src[4];   // (M g)(i,c) = M(i, col_src[c]) * col_val[c]
        cplx col_val[4];
        int  row_src[4];   // (g M)(r,j) = row_val[r] * M(row_src[r], j)
        cplx row_val[4];

        SparseGamma() {}
        explicit SparseGamma(const Matrix4cd & g)
        {
          for(int c = 0; c < 4; ++c)
            for(int r = 0; r < 4; ++r)
              if(g(r,c) != kCZero) { col_src[c] = r;  col_val[c] = g(r,c); }
          for(int r = 0; r < 4; ++r)
            for(int c = 0; c < 4; ++c)
              if(g(r,c) != kCZero) { row_src[r] = c;  row_val[r] = g(r,c); }
        }
      };

      struct SparseAlgebra {
        SparseGamma gam[5];
        SparseGamma gamgam5[4];   // gamma^mu gamma^5
        SparseGamma gam5gam[4];   // gamma^5 gamma^mu
        SparseAlgebra()
        {
          for(int i = 0; i < 5; ++i) gam[i] = SparseGamma(DA().gam[i]);
          for(int i = 0; i < 4; ++i)
          {
            gamgam5[i] = SparseGamma(DA().gamgam5[i]);
            gam5gam[i] = SparseGamma(DA().gam5gam[i]);
          }
        }
      };

      const SparseAlgebra & SA()
      {
        static const SparseAlgebra s;
        return s;
      }

      // M * gamma
      Matrix4cd MulR(const Matrix4cd & m, const SparseGamma & g)
      {
        Matrix4cd out;
        for(int i = 0; i < 4; ++i)
          for(int c = 0; c < 4; ++c) out(i,c) = m(i, g.col_src[c]) * g.col_val[c];
        return out;
      }

      // gamma * M
      Matrix4cd MulL(const SparseGamma & g, const Matrix4cd & m)
      {
        Matrix4cd out;
        for(int r = 0; r < 4; ++r)
          for(int j = 0; j < 4; ++j) out(r,j) = g.row_val[r] * m(g.row_src[r], j);
        return out;
      }

      // out += s * m
      void AddScaled(Matrix4cd & out, const cplx & s, const Matrix4cd & m)
      {
        for(int i = 0; i < 4; ++i)
          for(int j = 0; j < 4; ++j) out(i,j) += s * m(i,j);
      }

      // (Pslash + M) * numerator of the Rarita-Schwinger propagator, contracted with
      // the pion momentum kpi on its first (left = true) or second Lorentz index,
      // divided by P^2 - (M - i Gamma/2)^2. Returns one matrix per free index.
      //   left : sum_i kpi_i [ g^{ij} - g^i g^j/3 - 2 P^i P^j/(3M^2) - (g^i P^j - g^j P^i)/(3M) ]
      //   right: the same with the pion momentum contracted on j
      void ContractedRS(const ModelParams & par, const V4 & P, const V4 & kpi, bool left,
          Matrix4cd out[4])
      {
        const SparseAlgebra & s = SA();
        const double M  = par.xmd;
        const double P2 = Dot(P, P);
        const double kP = Dot(kpi, P);
        const cplx   xmd_c = cplx(M, -0.5 * DeltaWidth(par, P2));
        const cplx   inv_D = 1.0 / (P2 - xmd_c * xmd_c);

        const Matrix4cd proj = Slash(P) + M * DA().id4;
        Matrix4cd PG[4];                       // proj * gamma^j
        for(int j = 0; j < 4; ++j) PG[j] = MulR(proj, s.gam[j]);
        Matrix4cd PK;                          // proj * kslash
        for(int j = 0; j < 4; ++j) AddScaled(PK, kMetric[j] * kpi[j], PG[j]);

        for(int j = 0; j < 4; ++j)
        {
          Matrix4cd C;
          AddScaled(C, kpi[j] - 2.0 * kP * P[j] / 3.0 / (M*M), proj);
          if(left)
          {
            // - proj kslash g^j / 3 - (P^j proj kslash - kP proj g^j) / (3M)
            AddScaled(C, -1.0 / 3.0, MulR(PK, s.gam[j]));
            AddScaled(C, -P[j] / (3.0 * M), PK);
            AddScaled(C,  kP   / (3.0 * M), PG[j]);
          }
          else
          {
            // - proj g^j kslash / 3 - (kP proj g^j - P^j proj kslash) / (3M),
            // with g^j kslash = 2 k^j - kslash g^j
            AddScaled(C, -2.0 * kpi[j] / 3.0, proj);
            AddScaled(C,  1.0 / 3.0, MulR(PK, s.gam[j]));
            AddScaled(C, -kP   / (3.0 * M), PG[j]);
            AddScaled(C,  P[j] / (3.0 * M), PK);
          }
          out[j] = inv_D * C;
        }
      }

      //---------------------------------------------------------------------------
      // Delta current operators of diagrams a, b, c, d (Fortran: det_JaJb_JcJd)
      //   i_fl = 1 direct, 2 exchange
      //
      // With the gamma-N-Delta vertex (Delta index a, current index mu)
      //   V(a,mu) = C3 (g^{a mu} qslash - q^a gamma^mu)
      //           + C4/m (g^{a mu} q.P - q^a P^mu) U + C5/m (g^{a mu} q.p - q^a p^mu) U
      //   A(a,mu) = C5A m g^{a mu}
      // the sums over the Delta index collapse to
      //   J_a^mu = RS[mu] B1 - Rq B2(mu)          (N -> Delta: V gamma5 + A)
      //   J_b^mu = B1' RS[mu] - B2'(mu) Rq        (Delta -> N: gamma5 V + A)
      // where Rq = sum_a g_aa q^a RS[a]. U is the unit matrix (or the all-ones
      // matrix, see ModelParams::c4c5_broadcast).
      //---------------------------------------------------------------------------
      struct DeltaOps { Matrix4cd Ja[4], Jb[4], Jc[4], Jd[4]; };

      void DeltaDiagram(const ModelParams & par, const FormFactors & ff, const Kinematics & k,
          const V4 & Pdel, const V4 & pnuc, const V4 & kpi, bool n_to_delta, double norm,
          Matrix4cd J[4])
      {
        const DiracAlgebra & d = DA();
        const SparseAlgebra & s = SA();

        Matrix4cd RS[4];
        ContractedRS(par, Pdel, kpi, n_to_delta, RS);
        Matrix4cd Rq;
        for(int a = 0; a < 4; ++a) AddScaled(Rq, kMetric[a] * k.q[a], RS[a]);

        const double c4 = ff.cv4 / par.xmn;
        const double c5 = ff.cv5 / par.xmn;
        const double s1 = c4 * Dot(k.q, Pdel) + c5 * Dot(k.q, pnuc);
        const bool   has_c45 = (c4 != 0.0 || c5 != 0.0);

        // U gamma5 (N -> Delta) or gamma5 U (Delta -> N)
        Matrix4cd Ug5;
        if(has_c45)
        {
          if(!par.c4c5_broadcast) Ug5 = d.gam[4];
          else Ug5 = n_to_delta ? d.ones4 * d.gam[4] : d.gam[4] * d.ones4;
        }

        // B1 = C3 qslash gamma5 + s1 U gamma5 + C5A m   (or gamma5 qslash ...)
        Matrix4cd B1 = n_to_delta ? MulR(k.q_sl, s.gam[4]) : MulL(s.gam[4], k.q_sl);
        B1 *= cplx(ff.cv3);
        if(has_c45) AddScaled(B1, s1, Ug5);
        if(ff.ca5 != 0.0) AddScaled(B1, ff.ca5 * par.xmn, d.id4);

        // Rq gamma5 (or gamma5 Rq), needed for the C4/C5 part of B2 when U = 1
        Matrix4cd Rq5;
        if(has_c45 && !par.c4c5_broadcast)
          Rq5 = n_to_delta ? MulR(Rq, s.gam[4]) : MulL(s.gam[4], Rq);

        for(int mu = 0; mu < 4; ++mu)
        {
          // B2(mu) = C3 gamma^mu gamma5 + (C4/m P^mu + C5/m p^mu) U gamma5
          const double s2 = c4 * Pdel[mu] + c5 * pnuc[mu];

          J[mu] = n_to_delta ? RS[mu] * B1 : B1 * RS[mu];
          AddScaled(J[mu], -ff.cv3,
              n_to_delta ? MulR(Rq, s.gamgam5[mu]) : MulL(s.gam5gam[mu], Rq));
          if(has_c45)
          {
            if(!par.c4c5_broadcast) AddScaled(J[mu], -s2, Rq5);
            else AddScaled(J[mu], -s2, n_to_delta ? Rq * Ug5 : Ug5 * Rq);
          }
          J[mu] *= cplx(norm);
        }
      }

      void DetJaJbJcJd(const ModelParams & par, const FormFactors & ff, const Kinematics & k,
          int i_fl, DeltaOps & ops)
      {
        V4 pa, pb, pc, pd, k_1, k_2, ppn1, ppn2;
        pa = Add(k.p1, k.q);
        pc = Add(k.p2, k.q);
        if(i_fl == 1)
        {
          pb = Sub(k.pp1, k.q);  pd = Sub(k.pp2, k.q);
          k_1 = k.k1;  k_2 = k.k2;  ppn1 = k.pp1;  ppn2 = k.pp2;
        }
        else
        {
          pb = Sub(k.pp2, k.q);  pd = Sub(k.pp1, k.q);
          k_1 = k.k1e;  k_2 = k.k2e;  ppn1 = k.pp2;  ppn2 = k.pp1;
        }

        const double lpi2   = par.lpi   * par.lpi;
        const double lpind2 = par.lpind * par.lpind;
        const double xmpi2  = par.xmpi  * par.xmpi;
        const double fpik1   = (lpi2 - xmpi2) / (lpi2 - Dot(k_1, k_1));
        const double fpik2   = (lpi2 - xmpi2) / (lpi2 - Dot(k_2, k_2));
        const double fpindk1 = lpind2 / (lpind2 - Dot(k_1, k_1));
        const double fpindk2 = lpind2 / (lpind2 - Dot(k_2, k_2));

        const double norm = std::sqrt(par.fpinn2) * par.fstar / xmpi2 / par.xmn;
        const double norm_ab = fpik2 * fpindk2 * norm;
        const double norm_cd = fpik1 * fpindk1 * norm;

        // In the direct term a and b multiply the spectator spin trace of a
        // single pion vertex, Tr[(p2slash + m) gamma5 k2slash] = 0
        if(i_fl != 1)
        {
          DeltaDiagram(par, ff, k, pa, k.p1, k_2, true,  norm_ab, ops.Ja);
          DeltaDiagram(par, ff, k, pb, ppn1, k_2, false, norm_ab, ops.Jb);
        }
        DeltaDiagram(par, ff, k, pc, k.p2, k_1, true,  norm_cd, ops.Jc);
        DeltaDiagram(par, ff, k, pd, ppn2, k_1, false, norm_cd, ops.Jd);
      }

      //---------------------------------------------------------------------------
      // Pion-in-flight, seagull and pion-pole operators (Fortran: det_Jpi).
      // The operators acting on the struck-nucleon line (1) and on the
      // spectator line (2) are summed: they share their isospin coefficient.
      //---------------------------------------------------------------------------
      struct PionOps { Matrix4cd J1[4], J2[4]; };

      void DetJpi(const ModelParams & par, const FormFactors & ff, const Kinematics & k,
          int i_fl, PionOps & ops)
      {
        const DiracAlgebra & d = DA();
        const V4 & k_1 = (i_fl == 1) ? k.k1 : k.k1e;
        const V4 & k_2 = (i_fl == 1) ? k.k2 : k.k2e;
        const Matrix4cd & Pi_1 = (i_fl == 1) ? k.Pi_k1 : k.Pi_k1e;

        const double lpi2  = par.lpi  * par.lpi;
        const double xmpi2 = par.xmpi * par.xmpi;
        const double fpik1 = (lpi2 - xmpi2) / (lpi2 - Dot(k_1, k_1));
        const double fpik2 = (lpi2 - xmpi2) / (lpi2 - Dot(k_2, k_2));
        const double c     = par.fpinn2 / xmpi2;

        for(int mu = 0; mu < 4; ++mu)
        {
          // pion in flight + seagull on line 1, seagull on line 2
          ops.J1[mu].setZero();
          AddScaled(ops.J1[mu], ff.fpiem * (k_1[mu] - k_2[mu]) * fpik1 * fpik2 * c, Pi_1);
          AddScaled(ops.J1[mu], -ff.fpiem * fpik2 * fpik2 * c, d.gam5gam[mu]);
          ops.J2[mu].setZero();
          AddScaled(ops.J2[mu],  ff.fpiem * fpik1 * fpik1 * c, d.gam5gam[mu]);
        }

        if(k.ax == 0) return;

        // axial seagull and pion-pole pieces
        const double xmrho2 = par.xmrho * par.xmrho;
        const double frho1 = 1.0 / (1.0 - Dot(k_1, k_1) / xmrho2);
        const double frho2 = 1.0 / (1.0 - Dot(k_2, k_2) / xmrho2);
        const double pole  = Dot(k.q, k.q) - xmpi2;
        for(int mu = 0; mu < 4; ++mu)
        {
          AddScaled(ops.J1[mu], -frho1 / par.ga * fpik2 * fpik2 * c, d.gam[mu]);
          AddScaled(ops.J2[mu],  frho2 / par.ga * fpik1 * fpik1 * c, d.gam[mu]);
          AddScaled(ops.J1[mu],  frho1 / par.ga * k.q[mu] / pole * c * fpik2 * fpik2, k.q_sl);
          AddScaled(ops.J2[mu], -frho2 / par.ga * k.q[mu] / pole * c * fpik1 * fpik1, k.q_sl);
        }
      }

      // Isospin-independent two-body operators of one phase-space point
      struct TwoBodyOps {
        DeltaOps del_dir, del_exc;
        PionOps  pi_dir,  pi_exc;
      };

      //---------------------------------------------------------------------------
      // Two-body matrix elements for a spectator of isospin t2, spectator spin
      // summed, direct minus exchange
      // (Fortran: twobody_del_curr_matrix_el + twobody_pi_curr_matrix_el)
      //---------------------------------------------------------------------------
      void TwoBodyMatrixElements(const Kinematics & k, const TwoBodyOps & ops,
          const Vector2cd & t2, SpinCurrent & J)
      {
        const bool charged = (k.ax == 1 && !k.nc);

        for(int i1 = 0; i1 < 2; ++i1)
          for(int f1 = 0; f1 < 2; ++f1) J[f1][i1].fill(kCZero);

        // ---- direct: <pp1 pp2| j |p1 p2>, pp2 = p2 -------------------------------
        {
          const IsospinCoeffs iso = Isospin(charged, k.t1, t2, k.t1p, t2);

          cplx J_1[2][2];   // pion vertex of line 1
          for(int i1 = 0; i1 < 2; ++i1)
            for(int f1 = 0; f1 < 2; ++f1)
              J_1[f1][i1] = Sandwich(k.ubarpp1[f1], k.Pi_k1, k.up1[i1]);
          // Operators on line 1 (diagrams a, b, pion in flight and the seagull of
          // line 1) come with the spectator spin trace of the pion vertex of
          // line 2, Tr[(p2slash + m) gamma5 k2slash] = 0: nothing to add.

          for(int mu = 0; mu < 4; ++mu)
          {
            // trace of the operator on line 2 (times the pion vertex of line 1)
            Matrix4cd O2;
            AddScaled(O2, iso.c, ops.del_dir.Jc[mu]);
            AddScaled(O2, iso.d, ops.del_dir.Jd[mu]);
            AddScaled(O2, iso.pi, ops.pi_dir.J2[mu]);
            cplx T2 = kCZero;
            for(int i2 = 0; i2 < 2; ++i2) T2 += Sandwich(k.ubarpp2[i2], O2, k.up2[i2]);

            for(int i1 = 0; i1 < 2; ++i1)
            {
              for(int f1 = 0; f1 < 2; ++f1)
              {
                J[f1][i1][mu] += 0.5 * T2 * J_1[f1][i1];
              }
            }
          }
        }

        // ---- exchange: <pp2 pp1| j |p1 p2> -----------------------------------------
        {
          const IsospinCoeffs iso = Isospin(charged, k.t1, t2, t2, k.t1p);

          cplx Je_1[2][2], Je_2[2][2];   // [f][i]: <pp2 f|Pi|p1 i>, <pp1 f|Pi|p2 i>
          for(int i = 0; i < 2; ++i)
          {
            for(int f = 0; f < 2; ++f)
            {
              Je_1[f][i] = Sandwich(k.ubarpp2[f], k.Pi_k1e, k.up1[i]);
              Je_2[f][i] = Sandwich(k.ubarpp1[f], k.Pi_k2e, k.up2[i]);
            }
          }

          for(int mu = 0; mu < 4; ++mu)
          {
            Matrix4cd O1;   // between pp2 and p1
            AddScaled(O1, iso.a, ops.del_exc.Ja[mu]);
            AddScaled(O1, iso.b, ops.del_exc.Jb[mu]);
            AddScaled(O1, iso.pi, ops.pi_exc.J1[mu]);
            Matrix4cd O2;   // between pp1 and p2
            AddScaled(O2, iso.c, ops.del_exc.Jc[mu]);
            AddScaled(O2, iso.d, ops.del_exc.Jd[mu]);
            AddScaled(O2, iso.pi, ops.pi_exc.J2[mu]);

            cplx X[2][2], Y[2][2];
            for(int i = 0; i < 2; ++i)
            {
              for(int f = 0; f < 2; ++f)
              {
                X[f][i] = Sandwich(k.ubarpp2[f], O1, k.up1[i]);
                Y[f][i] = Sandwich(k.ubarpp1[f], O2, k.up2[i]);
              }
            }

            for(int i1 = 0; i1 < 2; ++i1)
            {
              for(int f1 = 0; f1 < 2; ++f1)
              {
                cplx sum = kCZero;
                for(int i2 = 0; i2 < 2; ++i2)
                  sum += X[i2][i1] * Je_2[f1][i2] + Y[f1][i2] * Je_1[i2][i1];
                J[f1][i1][mu] -= 0.5 * sum;
              }
            }
          }
        }

        // Spectator spinors are normalised to 2 E
        const double inv_2E2 = 1.0 / (2.0 * k.p2[0]);
        for(int i1 = 0; i1 < 2; ++i1)
          for(int f1 = 0; f1 < 2; ++f1)
            for(int mu = 0; mu < 4; ++mu) J[f1][i1][mu] *= inv_2E2;
      }

      void OneBodyMatrixElements(const ModelParams & par, const FormFactors & ff,
          const Kinematics & k, SpinCurrent & J1b)
      {
        Matrix4cd J_1b[4];
        DetJ1(par, ff, k, J_1b);
        for(int i1 = 0; i1 < 2; ++i1)
          for(int f1 = 0; f1 < 2; ++f1)
            for(int mu = 0; mu < 4; ++mu)
              J1b[f1][i1][mu] = Sandwich(k.ubarpp1[f1], J_1b[mu], k.up1[i1]);
      }

      void BuildTwoBodyOps(const ModelParams & par, const FormFactors & ff,
          const Kinematics & k, TwoBodyOps & ops)
      {
        DetJaJbJcJd(par, ff, k, 1, ops.del_dir);
        DetJpi(par, ff, k, 1, ops.pi_dir);
        DetJaJbJcJd(par, ff, k, 2, ops.del_exc);
        DetJpi(par, ff, k, 2, ops.pi_exc);
      }

    }  // anonymous namespace


    //---------------------------------------------------------------------------
    void ComputeCurrents(const ModelParams & par, const FormFactors & ff,
        const double p1[4], const double pp1[4],
        const double p2[4], const double q[4],
        int pdg_in, int pdg_out, int pdg_spect,
        bool has_axial, bool nc,
        SpinCurrent & J1b, SpinCurrent & J2b)
    {
      Kinematics k;
      InitKinematics(par, p1, pp1, p2, q, pdg_in, pdg_out, has_axial, nc, k);
      OneBodyMatrixElements(par, ff, k, J1b);

      TwoBodyOps ops;
      BuildTwoBodyOps(par, ff, k, ops);
      TwoBodyMatrixElements(k, ops, Isospinor(pdg_spect), J2b);
    }

    //---------------------------------------------------------------------------
    void ComputeCurrents(const ModelParams & par, const FormFactors & ff,
        const double p1[4], const double pp1[4],
        const double p2[4], const double q[4],
        int pdg_in, int pdg_out, bool has_axial, bool nc,
        SpinCurrent & J1b, SpinCurrent & J2b_pspect, SpinCurrent & J2b_nspect)
    {
      Kinematics k;
      InitKinematics(par, p1, pp1, p2, q, pdg_in, pdg_out, has_axial, nc, k);
      OneBodyMatrixElements(par, ff, k, J1b);

      // The operators do not depend on the isospin of the spectator
      TwoBodyOps ops;
      BuildTwoBodyOps(par, ff, k, ops);
      TwoBodyMatrixElements(k, ops, Isospinor(2212), J2b_pspect);
      TwoBodyMatrixElements(k, ops, Isospinor(2112), J2b_nspect);
    }

    //---------------------------------------------------------------------------
    void InterferenceTensor(const SpinCurrent & J1b, const SpinCurrent & J2b,
        std::complex<double> R[4][4])
    {
      for(int mu = 0; mu < 4; ++mu)
      {
        for(int nu = 0; nu < 4; ++nu)
        {
          R[mu][nu] = kCZero;
          for(int i1 = 0; i1 < 2; ++i1)
          {
            for(int f1 = 0; f1 < 2; ++f1)
            {
              R[mu][nu] += J2b[f1][i1][mu] * std::conj(J1b[f1][i1][nu])
                + std::conj(J2b[f1][i1][mu]) * J1b[f1][i1][nu];
            }
          }
        }
      }
    }

    //---------------------------------------------------------------------------
    void InterferenceHadronTensor(const SpinCurrent & J1b, const SpinCurrent & J2b,
        std::complex<double> A[4][4])
    {
      for(int mu = 0; mu < 4; ++mu)
      {
        for(int nu = 0; nu < 4; ++nu)
        {
          A[mu][nu] = kCZero;
          for(int i1 = 0; i1 < 2; ++i1)
          {
            for(int f1 = 0; f1 < 2; ++f1)
            {
              A[mu][nu] += std::conj(J1b[f1][i1][mu]) * J2b[f1][i1][nu]
                + std::conj(J2b[f1][i1][mu]) * J1b[f1][i1][nu];
            }
          }
        }
      }
    }

    //---------------------------------------------------------------------------
    void SquaredHadronTensor(const SpinCurrent & J, std::complex<double> A[4][4])
    {
      for(int mu = 0; mu < 4; ++mu)
      {
        for(int nu = 0; nu < 4; ++nu)
        {
          A[mu][nu] = kCZero;
          for(int i1 = 0; i1 < 2; ++i1)
            for(int f1 = 0; f1 < 2; ++f1)
              A[mu][nu] += std::conj(J[f1][i1][mu]) * J[f1][i1][nu];
        }
      }
    }

  }  // namespace twobody_currents_sf
}  // namespace genie
