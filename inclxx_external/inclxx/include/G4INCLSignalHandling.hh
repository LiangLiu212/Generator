#ifndef G4INCLSignalHandling_hh
#define G4INCLSignalHandling_hh

#ifdef INCL_SIGNAL_HANDLING

bool INTTERMsignalled = false, HUPsignalled = false;

#if defined(__linux__) && defined(__GNUC__)
#include <csignal>
#include <cstdlib>
#include <features.h>
#include <fenv.h>
#include <iostream>
#include <unistd.h>
#ifdef INCL_ROOT_USE
#include "TSystem.h"
#endif

struct sigaction fpeaction, hupaction, inttermaction, oldaction;

void FPESignalHandler(int /*sig*/, siginfo_t *siginfo, void * /*context*/) {
  const char *message;
  switch (siginfo->si_code) {
  case FPE_INTDIV:
    message = "ERROR - Integer divide by zero.          \n";
    break;
  case FPE_INTOVF:
    message = "ERROR - Integer overflow.                \n";
    break;
  case FPE_FLTDIV:
    message = "ERROR - Floating point divide by zero.   \n";
    break;
  case FPE_FLTOVF:
    message = "ERROR - Floating point overflow.         \n";
    break;
  case FPE_FLTUND:
    message = "ERROR - Floating point underflow.        \n";
    break;
  case FPE_FLTRES:
    message = "ERROR - Floating point inexact result.   \n";
    break;
  case FPE_FLTINV:
    message = "ERROR - Floating point invalid operation.\n";
    break;
  case FPE_FLTSUB:
    message = "ERROR - Subscript out of range.          \n";
    break;
  default:
    message = "ERROR - Unknown floating-point exception.\n";
    break;
  }
  size_t size = 42;
  size_t sizew = (size_t)::write(STDERR_FILENO, message,
                                 size); // the answer to the ultimate question
                                        // on Life, the Universe and Everything
  if (sizew != size) {
    INCL_ERROR("FPESignalHandler");
  }
  (void)::raise(SIGFPE);
}

void INTTERMSignalHandler(int /*sig*/) {
  if (INTTERMsignalled) {
    const char *message =
        "Caught SIGINT while waiting for an event to finish, exiting.\n";
    size_t size = 61;
    size_t sizew = (size_t)::write(STDERR_FILENO, message, size);
    if (sizew != size) {
      INCL_ERROR("INTTERMSignalHandler");
    }
    ::_exit(1);
  } else {
    INTTERMsignalled = true;
    const char *message = "Caught SIGINT, will try to quit cleanly at the end "
                          "of the current event.\n";
    size_t size = 74;
    size_t sizew = (size_t)::write(STDERR_FILENO, message, size);
    if (sizew != size) {
      INCL_ERROR("INTTERMSignalHandler");
    }
  }
}

void HUPSignalHandler(int /*sig*/) {
  HUPsignalled = true;
  const char *message = "Caught SIGHUP, will flush log and output files at the "
                        "end of the current event.\n";
  size_t size = 81;
  size_t sizew = (size_t)::write(STDERR_FILENO, message, size);
  if (sizew != size) {
    INCL_ERROR("HUPSignalHandler");
  }
}

void enableSignalHandling() {
#ifdef INCL_ROOT_USE
  // We don't want ROOT's signal handlers
  for (int sig = 0; sig < kMAXSIGNALS; ++sig) {
    gSystem->ResetSignal((ESignals)sig);
  }
#endif

  (void)feenableexcept(FE_DIVBYZERO);
  (void)feenableexcept(FE_INVALID);
  //(void) feenableexcept( FE_OVERFLOW );
  //(void) feenableexcept( FE_UNDERFLOW );

  // Set the handler for SIGFPE
  sigset_t *def_set;
  def_set = &fpeaction.sa_mask;
  (void)sigfillset(def_set);
  (void)sigdelset(def_set, SIGFPE);
  fpeaction.sa_sigaction = FPESignalHandler;
  fpeaction.sa_flags = SA_SIGINFO | SA_RESETHAND;
  (void)sigaction(SIGFPE, &fpeaction, &oldaction);

  // Set the handler for SIGHUP
  sigset_t *def_set_hup;
  def_set_hup = &hupaction.sa_mask;
  (void)sigfillset(def_set_hup);
  (void)sigdelset(def_set_hup, SIGHUP);
  hupaction.sa_handler = HUPSignalHandler;
  hupaction.sa_flags = 0;
  (void)sigaction(SIGHUP, &hupaction, &oldaction);

  // Set the handler for SIGINT and SIGTERM
  sigset_t *def_set_intterm;
  def_set_intterm = &inttermaction.sa_mask;
  (void)sigfillset(def_set_intterm);
  (void)sigdelset(def_set_intterm, SIGINT);
  (void)sigdelset(def_set_intterm, SIGTERM);
  inttermaction.sa_handler = INTTERMSignalHandler;
  inttermaction.sa_flags = 0;
  (void)sigaction(SIGINT, &inttermaction, &oldaction);
  (void)sigaction(SIGTERM, &inttermaction, &oldaction);
}
#else
void enableSignalHandling() {
  std::cerr << "Signal handling is not supported by your platform."
            << std::endl;
}
#endif

#endif

#endif
