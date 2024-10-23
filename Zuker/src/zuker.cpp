#include "zuker.hpp"

/*
 * Constructor
 * @param seq: RNA sequence
 */
Zuker::Zuker(const std::string& seq) : seq_(seq) {
  const int n = seq_.size();

  if (seq.find_first_not_of("AUGC") != std::string::npos) {
    throw std::invalid_argument(
        "Zuker: Invalid nucleotide: sequence must be consist of AUGC, but "
        "got " +
        seq);
  }

  W_.resize(n, std::vector<double>(n, 0));
  V_.resize(n, std::vector<double>(n, kInf));
  M_.resize(n, std::vector<double>(n, kInf));
}

/*
 * Get bracket pair from secondary structure
 * @param s: secondary structure
 * @return: map of bracket pair
 */
std::unordered_map<int, int> Zuker::getBracketPair(const std::string& s) const {
  std::stack<int> stk;
  std::unordered_map<int, int> pair;
  for (size_t i = 0; i < s.length(); ++i) {
    if (s.at(i) == '(') {
      stk.push(i);
    } else if (s.at(i) == ')') {
      if (stk.empty()) {
        throw std::invalid_argument("Invalid secondary structure");
      }
      pair.emplace(stk.top(), i);
      stk.pop();
    }
  }
  return pair;
}

/*
 * Calculate True Positive(TP), False Positive(FP), False Negative(FN)
 * @param pred: predicted secondary structure
 * @return: tuple of TP, FP, FN
 */
std::tuple<int, int, int> Zuker::countTPFPFN(const std::string& answer) const {
  std::string pred = traceback();
  int TP = 0, FP = 0, FN = 0;

  if (pred.length() != answer.length()) {
    throw std::invalid_argument(
        "Length of answer and prediction must be the same");
  }

  std::stack<int> stk;
  std::unordered_map<int, int> pairPred = getBracketPair(pred);
  std::unordered_map<int, int> pairAnswer = getBracketPair(answer);

  for (const auto& [i, j] : pairPred) {
    if (pairAnswer.find(i) != pairAnswer.end() && pairAnswer.at(i) == j) {
      ++TP;
    } else {
      ++FP;
    }
  }

  for (const auto& [i, j] : pairAnswer) {
    if (pairPred.find(i) == pairPred.end() || pairPred.at(i) != j) {
      ++FN;
    }
  }

  return std::make_tuple(TP, FP, FN);
}

/*
 * Free energy parameters for hairpin, internal, bulge loop
 * @param i: loop type (0: hairpin, 1: internal, 2: bulge)
 * @param j: loop size
 */
double Zuker::getLoopData(int i, int j) const {
  assert(j < 31);

  const double loopData[3][31] = {
      // Hairpin
      {kInf, kInf, 4.4, 4.6, 4.7, 4.4, 5.0, 4.5, 5.4, 5.5, 5.6,
       5.7,  5.8,  5.9, 5.9, 6.0, 6.1, 6.1, 6.2, 6.2, 6.3, 6.3,
       6.4,  6.4,  6.5, 6.5, 6.5, 6.6, 6.6, 6.7, kInf},
      // Internal
      {kInf, 1.0, 1.0, 1.1, 2.0, 2.0, 2.1, 2.3, 2.4, 2.5, 2.6,
       2.7,  2.8, 2.9, 2.9, 3.0, 3.1, 3.1, 3.2, 3.3, 3.3, 3.4,
       3.4,  3.5, 3.5, 3.5, 3.6, 3.6, 3.7, 3.7, kInf},
      // Bulge
      {1.0, 1.0, 1.0, 1.1, 2.0, 2.0, 2.1, 2.3, 2.4, 2.5, 2.6,
       2.7, 2.8, 2.9, 2.9, 3.0, 3.1, 3.1, 3.2, 3.3, 3.3, 3.4,
       3.4, 3.5, 3.5, 3.5, 3.6, 3.6, 3.7, 3.7, kInf}};

  return loopData[i][j];
}

/*
 * Free energy calculation for hairpin loop
 * @param i: start index of hairpin loop
 * @param j: end index of hairpin loop
 */
double Zuker::F1(int i, int j) const {
  if (j - i - 1 > 30) {
    return kInf;
  }
  return getLoopData(0, (j - i - 1) - 1);
}

/*
 * Classify loop type; Helper function for F2
 * @param i: start index of loop
 * @param j: end index of loop
 * @param h: start index of internal loop
 * @param l: end index of internal loop
 */
LoopType Zuker::classifyLoop(int i, int j, int h, int l) const {
  if (h == i + 1 && l == j - 1) {
    return LoopType::STACKING;
  } else if (i + 1 < h && l == j - 1) {
    return LoopType::BULGE;
  } else if (h == i + 1 && l < j - 1) {
    return LoopType::BULGE;
  } else if (i + 1 < h && l < j - 1) {
    return LoopType::INTERNAL;
  }
  throw std::invalid_argument("classifyLoop: Unexpected loop type");
}

/*
 * Classify base pair based on the indices i and j
 * @param i: index i
 * @param j: index j
 */
BasePair Zuker::classifyBases(int i, int j) const {
  if ((seq_.at(i) == 'G' && seq_.at(j) == 'C') ||
      (seq_.at(i) == 'C' && seq_.at(j) == 'G')) {
    return BasePair::GC;
  } else if ((seq_.at(i) == 'A' && seq_.at(j) == 'U') ||
             (seq_.at(i) == 'U' && seq_.at(j) == 'A')) {
    return BasePair::AU;
  } else if ((seq_.at(i) == 'G' && seq_.at(j) == 'U') ||
             (seq_.at(i) == 'U' && seq_.at(j) == 'G')) {
    return BasePair::GU;
  } else {
    return BasePair::UNMATCH;
  }
}

/*
 * Free energy parameters for stacking loop
 * @param a: stacking loop type
 * @param b: internal loop type
 */
double Zuker::getStackLoopData(BasePair a, BasePair b) const {
  const double stackLoopData[4][4] = {
      // GC, AU, GU, UNMATCH
      {-3.0, -2.0, -2.0, 0},  // GC
      {-2.0, -0.5, -0.5, 0},  // AU
      {-2.0, -0.5, -0.5, 0},  // GU
      {0, 0, 0, 0}            // UNMATCH
  };

  return stackLoopData[a][b];
}

/*
 * Free energy calculation for stacking loop
 * @param i: start index of stacking loop
 * @param j: end index of stacking loop
 */
double Zuker::calcStackingEnergy(int i, int j) const {
  BasePair ij = classifyBases(i, j);
  BasePair hl = classifyBases(i + 1, j - 1);

  return getStackLoopData(ij, hl);
}

/*
 * Free energy calculation for stacking, internal, bulge loop;
 * Determine the type of loop and calculate the free energy
 */
double Zuker::F2(int i, int j, int h, int l) const {
  switch (classifyLoop(i, j, h, l)) {
    case LoopType::STACKING:
      return calcStackingEnergy(i, j);
    case LoopType::INTERNAL:
      return getLoopData(1, (h - i - 1) + (j - l - 1) - 1);
    case LoopType::BULGE:
      return getLoopData(2, (h - i - 1) + (j - l - 1) - 1);
    default:
      throw std::invalid_argument("F2: Unexpected loop type");
  }
}

/*
 * Update W table
 * @param i: start index
 * @param j: end index
 */
void Zuker::updateW(int i, int j) {
  const int seqLen = seq_.size();

  W_.at(i).at(j) = V_.at(i).at(j);

  if (i + 1 < seqLen) {
    utils::chmin(W_.at(i).at(j), W_.at(i + 1).at(j));
  }

  if (j > 0) {
    utils::chmin(W_.at(i).at(j), W_.at(i).at(j - 1));
  }

  for (int k = i; k < j; ++k) {
    utils::chmin(W_.at(i).at(j), W_.at(i).at(k) + W_.at(k + 1).at(j));
  }
}

void Zuker::updateV(int i, int j) {
  if (classifyBases(i, j) == BasePair::UNMATCH) {
    V_.at(i).at(j) = kInf;
    return;
  }

  V_.at(i).at(j) = F1(i, j);

  for (int h = i + 1; h < j - 1; ++h) {
    for (int l = j - 1; l > h; --l) {
      if ((h - i - 1) + (j - l - 1) > 30) break;
      utils::chmin(V_.at(i).at(j), F2(i, j, h, l) + V_.at(h).at(l));
    }
  }

  for (int k = i + 1; k < j - 1; ++k) {
    utils::chmin(V_.at(i).at(j),
                 M_.at(i + 1).at(k) + M_.at(k + 1).at(j - 1) + a);
  }
}

void Zuker::updateM(int i, int j) {
  const int seqLen = seq_.size();

  M_.at(i).at(j) = V_.at(i).at(j) + b;

  if (i + 1 < seqLen) {
    utils::chmin(M_.at(i).at(j), M_.at(i + 1).at(j) + c);
  }

  if (j > 0) {
    utils::chmin(M_.at(i).at(j), M_.at(i).at(j - 1) + c);
  }

  for (int k = i; k < j; ++k) {
    utils::chmin(M_.at(i).at(j), M_.at(i).at(k) + M_.at(k + 1).at(j));
  }
}

void Zuker::dp() {
  const int seqLen = seq_.size();

  /* Initialization */
  for (int i = 0; i < seqLen; ++i) {
    W_.at(i).at(i) = 0;
    V_.at(i).at(i) = kInf;
    M_.at(i).at(i) = kInf;

    if (i + 1 < seqLen) {
      W_.at(i).at(i + 1) = 0;
      V_.at(i).at(i + 1) = kInf;
      M_.at(i).at(i + 1) = kInf;
    }
    if (i + 2 < seqLen) {
      W_.at(i).at(i + 2) = 0;
      V_.at(i).at(i + 2) = kInf;
      M_.at(i).at(i + 2) = kInf;
    }
    if (i + 3 < seqLen) {
      W_.at(i).at(i + 3) = 0;
      V_.at(i).at(i + 3) = kInf;
      M_.at(i).at(i + 3) = kInf;
    }
  }

  /* Recursion */
  for (int i = seqLen - 5; i >= 0; --i) {
    int idx = i;
    for (int j = seqLen - 1; idx >= 0; --j, --idx) {
      updateV(idx, j);
      updateW(idx, j);
      updateM(idx, j);
    }
  }
}

std::string Zuker::traceback() const {
  std::stack<std::tuple<int, int, char>> stack;
  std::string result(seq_.size(), '.');

  stack.push({0, static_cast<int>(seq_.size()) - 1, 'W'});

  while (!stack.empty()) {
    auto [i, j, table] = stack.top();
    stack.pop();

    if (i >= j) continue;

    if (table == 'W') {
      if (std::abs(W_.at(i).at(j) - V_.at(i).at(j)) < kEps) {
        stack.push({i, j, 'V'});
      } else if (std::abs(W_.at(i).at(j) - W_.at(i + 1).at(j)) < kEps) {
        stack.push({i + 1, j, 'W'});
      } else if (std::abs(W_.at(i).at(j) - W_.at(i).at(j - 1)) < kEps) {
        stack.push({i, j - 1, 'W'});
      } else {
        bool found = false;
        for (int k = i; k < j; ++k) {
          if (std::abs(W_.at(i).at(j) - W_.at(i).at(k) - W_.at(k + 1).at(j)) <
              kEps) {
            stack.push({i, k, 'W'});
            stack.push({k + 1, j, 'W'});
            found = true;
            break;
          }
        }
        if (!found) {
          throw std::runtime_error("traceback: W table is not correct");
        }
      }
    } else if (table == 'V') {
      if (classifyBases(i, j) == BasePair::UNMATCH) {
        continue;
      }

      bool found = false;
      result[i] = '(';
      result[j] = ')';

      if (std::abs(V_.at(i).at(j) - F1(i, j)) < kEps) {
        continue;
      }

      for (int h = i + 1; h < j - 1; ++h) {
        for (int l = j - 1; l > h; --l) {
          if ((h - i - 1) + (j - l - 1) > 30) break;

          if (std::abs(V_.at(i).at(j) - (F2(i, j, h, l) + V_[h][l])) < kEps) {
            stack.push({h, l, 'V'});
            found = true;
            break;
          }
        }
      }

      if (!found) {
        for (int k = i + 1; k < j - 1; ++k) {
          if (std::abs(V_.at(i).at(j) -
                       (M_.at(i + 1).at(k) + M_.at(k + 1).at(j) + a))) {
            stack.push({i + 1, k, 'M'});
            stack.push({k + 1, j - 1, 'M'});
            found = true;
            break;
          }
        }
        if (!found) {
          throw std::runtime_error("traceback: V table is not correct");
        }
      }
    } else if (table == 'M') {
      if (std::abs(M_.at(i).at(j) - (V_.at(i).at(j) + b)) < kEps) {
        stack.push({i, j, 'V'});
      } else if (std::abs(M_.at(i).at(j) - (M_.at(i + 1).at(j) + c)) < kEps) {
        stack.push({i + 1, j, 'M'});
      } else if (std::abs(M_.at(i).at(j) - (M_.at(i).at(j - 1) + c)) < kEps) {
        stack.push({i, j - 1, 'M'});
      } else {
        bool found = false;
        for (int k = i; k < j; ++k) {
          if (std::abs(M_.at(i).at(j) - M_.at(i).at(k) - M_.at(k + 1).at(j)) <
              kEps) {
            stack.push({i, k, 'M'});
            stack.push({k + 1, j, 'M'});
            found = true;
            break;
          }
        }
        if (!found) {
          throw std::runtime_error("traceback: M table is not correct");
        }
      }
    }
  }

  return result;
}

double Zuker::calcSensitivity(const std::string& answer) const {
  const auto [TP, _, FN] = countTPFPFN(answer);
  return static_cast<double>(TP) / (TP + FN);
}

double Zuker::calcPPV(const std::string& answer) const {
  const auto [TP, FP, _] = countTPFPFN(answer);
  return static_cast<double>(TP) / (TP + FP);
}

const std::string& Zuker::getSeq() const { return seq_; }
