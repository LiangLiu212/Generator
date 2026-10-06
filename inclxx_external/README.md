# GENIE + INCL++ v6.34 for SBND (`genie v3_06_02_sbn4 -q e26:incl634:prof`)

This branch (`feature/genie-incl634-sbnd`) is `feature/for_Anna` plus this
directory. It carries INCL++ v6.34 (the GENIE-interface fork) as a **binary
external**: libraries, headers and run-time data tables, **no INCL++ source
code**. INCL++ is not a separate UPS product; the build installs it inside the
genie product, in `${GENIE_FQ_DIR}/inclxx`.

| Path | What it is |
|---|---|
| `inclxx/lib/` | the 9 INCL++ / de-excitation shared libraries |
| `inclxx/include/` | INCL++ headers needed to compile the GENIE interface |
| `inclxx/bin/` | `inclxx-config` (used by GENIE's `Make.include`), `thisinclxx.sh` |
| `inclxx/share/` | data tables read at run time (INCL, ABLA07, ABLA++, GEMINI++) |
| `inclxx/BUILD_INFO`, `inclxx/LICENSE.pdf` | provenance of the binaries, INCL++ rules of use |
| `sbnd/` | UPS build recipe (ssibuildshims) and build / test drivers for SBND |
| `build_inclxx_external.sh` | how `inclxx/` was produced (needs the INCL++ source) |

The binaries are built with the `e26:prof` toolchain of sbndcode v10_14_02_05
(SL7, gcc 12.1.0, ROOT v6_28_12, boost v1_82_0). They only fit a GENIE built
with the same toolchain, which is what the recipe below does.

INCL++ is distributed under its own rules (`inclxx/LICENSE.pdf`): it may not be
passed on to other people without the authorization of its authors.

## Build the genie UPS product in the SBND environment

Work on an SBND gpvm (or any machine with cvmfs and apptainer). The scripts
enter the SL7 container themselves.

1. Get the branch.

   ```bash
   git clone -b feature/genie-incl634-sbnd <url of this repository> Generator
   cd Generator/inclxx_external/sbnd
   ```

2. Choose a UPS products area to build into. Any writable directory with a
   `.upsfiles/` directory works, for example an `mrb` `localProducts_*`
   directory, or a new one:

   ```bash
   PROD=/path/to/my_products
   mkdir -p $PROD && cp -r /cvmfs/larsoft.opensciencegrid.org/products/.upsfiles $PROD/
   ```

3. Build (about 25 minutes on one core).

   ```bash
   ./build_genie_sbn4_incl634.sh $PROD > build.log 2>&1
   ```

   This runs `bootstrap.sh $PROD` and `build_genie.sh $PROD e26:incl634 prof tar`.
   It builds the **committed** state of the branch checked out in this clone
   (it makes a `git archive` of it); uncommitted changes are not built. It
   first deletes any previous `$PROD/genie/v3_06_02_sbn4`.

   The log must end with `BUILD DONE rc=0`. The results are

   - `$PROD/genie/v3_06_02_sbn4/` and `$PROD/genie/v3_06_02_sbn4.version/` — the declared product,
   - `$PROD/genie-3.06.02.sbn4-sl7-x86_64-e26-incl634-prof.tar.bz2` — the binary
     tarball (product directory, table file and version file),
   - `$PROD/genie-3.06.02.sbn4-source.tar.bz2` — GENIE / Reweight source and the recipe.

4. Test.

   ```bash
   ./test_genie_sbn4_incl634.sh $PROD    # unpack the tarball elsewhere, run gevgen from there
   ./test_lar_sbn4_incl634.sh   $PROD    # GENIEGen job with sbndcode
   ```

   Both must end with `RESULT: PASS`. The first one unpacks the binary tarball
   into an empty products area and checks with `strace` that INCL++ is read
   only from inside the product, which is what matters for a cvmfs
   installation. Both need the INCL26 cross-section splines (`genie_xsec
   v3_06_00 -q INCL2609a00000:k250:e1000`); set `XSEC_PROD` to the products area
   that holds them (default: Liang's area on the SBND gpvms).

5. Install somewhere else (for example on cvmfs): unpack the binary tarball in
   the top directory of the products area. Nothing else is needed; the product
   has no absolute paths to the build area.

   ```bash
   tar -xjf genie-3.06.02.sbn4-sl7-x86_64-e26-incl634-prof.tar.bz2 -C /path/to/products
   ```

Options (environment variables): `GENIE_SOURCE_URL` and `GENIE_SOURCE_REF`
build another repository / branch instead of this clone, for example
`GENIE_SOURCE_URL=https://github.com/<user>/Generator.git
GENIE_SOURCE_REF=feature/genie-incl634-sbnd`; `SBNDCODE_VERSION` selects the
sbndcode release (default `v10_14_02_05`).

The product flavor is the one UPS reports on the build host (`Linux64bit+5.14-2.17`
on the AL9 gpvms with the SL7 container), as for the other SBN genie builds.

## Use it with sbndcode

```bash
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup sbndcode v10_14_02_05 -q prof:e26
export PRODUCTS=$PROD:$PRODUCTS          # not needed once the product is on cvmfs
unsetup genie
setup genie v3_06_02_sbn4 -q e26:incl634:prof
unsetup genie_xsec
setup genie_xsec v3_06_00 -q INCL2609a00000:k250:e1000     # or INCL2607a00000:k250:e1000
```

`setup genie` also defines the INCL++ environment: `INCLXX_DIR`
(`${GENIE_FQ_DIR}/inclxx`), `INCLXX_DATA_DIR` (`${INCLXX_DIR}/share`), and it
adds `${INCLXX_DIR}/lib` to `LD_LIBRARY_PATH`. sbndcode does not have to be
rebuilt: nugen loads the GENIE libraries by their unversioned names.

`gevgen` and the other GENIE applications work directly. **`lar` jobs need the
INCL++ libraries in `LD_PRELOAD`**: `lar` loads Geant4, whose own built-in
INCL++ uses the same symbol names, and without the preload the job crashes
when GENIE initialises INCL++.

```bash
INCL_PRELOAD=$(ls $INCLXX_DIR/lib/lib{INCL_Utils,INCL_IO,INCL_Physics,DeExcitation,ABLA07,ABLAXX,FERMI_BREAKUP,GEMINIXX,SMM}.so | tr '\n' ':')
LD_PRELOAD=${INCL_PRELOAD%:} lar -c my_genie.fcl -n 10
```

This is safe as long as the Geant4 physics list does not use the Geant4
INCLXX models.

## How the recipe handles INCL++

- `sbnd/bootstrap.sh` archives the GENIE branch; `inclxx_external/inclxx`
  comes with it. No INCL++ download, no UPS `inclxx` product.
- `sbnd/build_genie.sh` accepts the extra qualifier `incl634` (only as
  `e26:incl634` `prof`). For it, it moves `inclxx_external/inclxx` to
  `${GENIE_FQ_DIR}/inclxx` and configures GENIE with `--enable-incl
  --with-incl-inc=${GENIE_FQ_DIR}/inclxx/include
  --with-incl-lib=${GENIE_FQ_DIR}/inclxx/lib` plus the boost paths. The other
  configure options are those of the SBN genie builds, except that this recipe
  enables Pythia6 only (`--enable-pythia6 --disable-pythia8`), like
  `v3_06_02_incl`.
- `sbnd/ups/genie.table` is the table of `genie v3_06_02_sbn3` plus the
  `e26:incl634:prof` entry, which sets the INCL++ environment relative to
  `${GENIE_FQ_DIR}` and adds boost.

For a build without UPS, pass the same configure options as `build_genie.sh`
with the INCL++ paths pointing at `inclxx_external/inclxx/include` and
`inclxx_external/inclxx/lib`, and source `inclxx/bin/thisinclxx.sh` for the
run-time environment.

## Rebuilding the INCL++ external

Only needed to change the INCL++ version or its build options, and only
possible with access to the INCL++ source (the GENIE-interface fork
`LiangLiu212/inclxx`, not public).

```bash
cd inclxx_external
INCLXX_REPO=/path/to/inclxx INCLXX_REF=<commit> ./build_inclxx_external.sh > build_inclxx_external.log 2>&1
git status --short inclxx      # review, then commit
```

The script builds INCL++ with the `e26:prof` toolchain and replaces `inclxx/`.
It keeps the libraries, the headers, `inclxx-config`, `thisinclxx.sh` and the
data tables, and removes everything else that `make install` produces (the
standalone `INCLCascade` program and the GEMINI++ sources and documentation
that INCL++ installs under `share/`). `inclxx/BUILD_INFO` records the commit,
compiler and cmake options. Rebuild genie afterwards.
