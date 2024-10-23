#include <iostream>

#include "zuker.hpp"

int main() {
  try {
    const std::string seq =
        "ACGUUCGUAGCUCAAUAGGUUAGAGCAUUACCAUGACAUGGUAGAGGUUAGUGGUUCAAGUCCACUCGAA"
        "CGUA";

    Zuker zuker(seq);
    zuker.dp();
    std::cout << "Secondary structure for " << zuker.getSeq() << " is:\n"
              << zuker.traceback() << std::endl;

    const std::string answer =
        "(((((((..((((.........)))).(((((.......))))).....(((((.......)))))))))"
        "))).";
    std::cout << "Sensitivity: " << zuker.calcSensitivity(answer) << std::endl;
    std::cout << "PPV: " << zuker.calcPPV(answer) << std::endl;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
}
