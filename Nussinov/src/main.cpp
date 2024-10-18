#include <fstream>
#include <iostream>

#include "nussinov.hpp"

int main() {
  try {
    std::string seq, answer;

    std::ifstream ifs("data/seq.txt");

    if (ifs.is_open()) {
      std::string line;
      int lineCount = 0;

      while (std::getline(ifs, line)) {
        if (lineCount == 0 && line.substr(0, 9) == "sequence:") {
          seq = line.substr(10);
        } else if (lineCount == 1 && line.substr(0, 10) == "structure:") {
          answer = line.substr(11);
        }
        lineCount++;
      }
      ifs.close();
    } else {
      throw std::runtime_error("Cannot open file");
    }

    Nussinov nussinov(seq);
    nussinov.dp();
    std::cout << "Secondary structure for " << nussinov.getSeq() << " is:\n"
              << nussinov.traceback() << std::endl;

    std::cout << "Sensitivity: " << nussinov.calcSensitivity(answer)
              << std::endl;
    std::cout << "PPV: " << nussinov.calcPPV(answer) << std::endl;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
}
