#ifndef INCLTREECONVERTER_HH
#define INCLTREECONVERTER_HH

#include <iostream>
#include <string>

const std::string version = "0.4";

void printWelcomeMessage() {
  static std::string welcomeMessage =
      "INCLTreeConverter v" + version +
#ifdef HAS_THBOOK
      "\n\nConverts INCL4.6-style `101\' ntuples or `h101' trees in ROOT files to INCL++-style `et\' trees.\n";
#else
      "\n\nConverts INCL4.6-style `h101\' trees in ROOT files to INCL++-style `et\' trees.\n";
#endif

  std::cerr << welcomeMessage;
}

void printUsageMessage() {
#ifdef HAS_THBOOK
  std::cerr << "Usage: INCLTreeConverter [-h] -o <output_file> <root/hbk_file>:<text_file> [<root/hbk_file>:<text_file> [<root/hbk_file>:<text_file> ... ]]\n\n";
#else
  std::cerr << "Usage: INCLTreeConverter [-h] -o <output_file> <root_file>:<text_file> [<root_file>:<text_file> [<root_file>:<text_file> ... ]]\n\n";
#endif
}

#endif // INCLTREECONVERTER_HH
