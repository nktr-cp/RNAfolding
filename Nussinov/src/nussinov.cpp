#include "nussinov.hpp"

#ifdef DEBUG
#include <iomanip>
#include <iostream>
#endif  // DEBUG

Nussinov::Nussinov(std::string seq) : seq_(seq) {
  if (seq.find_first_not_of("AUGC") != std::string::npos) {
    throw std::invalid_argument(
        "Invalid nucleotide: sequence must be consist of AUGC, but got " + seq);
  }
  dp_.resize(seq_.size(), std::vector<int64_t>(seq_.size(), -kInf));
}

Nussinov::~Nussinov() {}

/*
 * Get bracket pair from secondary structure
 * @param s: secondary structure
 * @return: map of bracket pair
 */
std::unordered_map<int, int> Nussinov::getBracketPair(
    const std::string& s) const {
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
std::tuple<int, int, int> Nussinov::countTPFPFN(
    const std::string& answer) const {
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
 * delta function for Nussinov algorithm
 * @param i: index of nucleotide in sequence
 * @param j: index of nucleotide in sequence
 * @return: 1 if i-j is a valid base pair, 0 otherwise
 */
int64_t Nussinov::delta(int i, int j) const {
  if (j - i - 1 < 3) {
    return 0;
  }
  switch (seq_.at(i)) {
    case 'A':
      return seq_.at(j) == 'U' ? 1 : 0;
    case 'U':
      return (seq_.at(j) == 'A' || seq_.at(j) == 'G') ? 1 : 0;
    case 'G':
      return (seq_.at(j) == 'C' || seq_.at(j) == 'U') ? 1 : 0;
    case 'C':
      return seq_.at(j) == 'G' ? 1 : 0;
    default:
      throw std::invalid_argument("Invalid nucleotide");
  }
}

/*
 * Dynamic programming function for Nussinov algorithm
 * Time complexity: O(n^3)
 */
void Nussinov::dp() {
  const int seqLen = seq_.length();
  /* Initialization */
  for (int i = 0; i < seqLen; ++i) {
    dp_.at(i).at(i) = 0;
    if (i > 0) {
      dp_.at(i).at(i - 1) = 0;
    }
  }

  /* Recursion */
  for (int i = seqLen - 2; i >= 0; --i) {
    int idx = i;
    for (int j = seqLen - 1; idx >= 0; --j, --idx) {
      if (idx + 1 < seqLen) {
        dp_.at(idx).at(j) = std::max(dp_.at(idx).at(j), dp_.at(idx + 1).at(j));
      }

      if (j > 0) {
        dp_.at(idx).at(j) = std::max(dp_.at(idx).at(j), dp_.at(idx).at(j - 1));
      }

      if (idx + 1 < seqLen && j > 0) {
        dp_.at(idx).at(j) = std::max(dp_.at(idx).at(j),
                                     dp_.at(idx + 1).at(j - 1) + delta(idx, j));
      }

      for (int k = idx; k < j; ++k) {
        dp_.at(idx).at(j) = std::max(dp_.at(idx).at(j),
                                     dp_.at(idx).at(k) + dp_.at(k + 1).at(j));
      }
    }
  }

#ifdef DEBUG
  std::cout << "===> Nussinov: DP table <===" << std::endl;
  int maxWidth = 0;
  for (const auto& row : dp_) {
    for (const auto& col : row) {
      std::string colStr = std::to_string(col);
      if (col == -kInf) {
        colStr = "-Inf";
      }
      maxWidth = std::max(maxWidth, static_cast<int>(colStr.length()));
    }
  }

  for (const auto& row : dp_) {
    for (const auto& col : row) {
      std::cout << std::setw(maxWidth + 1) << std::right
                << (col == -kInf ? "-Inf" : std::to_string(col));
    }
    std::cout << std::endl;
  }
#endif  // DEBUG
}

/*
 * Traceback function for Nussinov algorithm
 * @return: secondary structure of the sequence
 */
std::string Nussinov::traceback() const {
  if (dp_.back().back() == -kInf) {
    throw std::runtime_error("DP table is not filled");
  }

  const int seqLen = seq_.length();
  std::stack<std::pair<int, int>> stk;

  /* Initialization */
  stk.emplace(0, seqLen - 1);

  /* Recursion */
  std::vector<std::pair<int, int>> basePair;
  while (!stk.empty()) {
    const auto [i, j] = stk.top();
    stk.pop();

    if (i >= j) continue;

    if (i + 1 < seqLen && dp_.at(i + 1).at(j) == dp_.at(i).at(j)) {
      stk.emplace(i + 1, j);
    } else if (j > 0 && dp_.at(i).at(j - 1) == dp_.at(i).at(j)) {
      stk.emplace(i, j - 1);
    } else if (i + 1 < seqLen && j > 0 &&
               dp_.at(i + 1).at(j - 1) + delta(i, j) == dp_.at(i).at(j)) {
      stk.emplace(i + 1, j - 1);
      basePair.emplace_back(i, j);
    } else {
      for (int k = i; k < j; ++k) {
        if (dp_.at(i).at(k) + dp_.at(k + 1).at(j) == dp_.at(i).at(j)) {
          stk.emplace(k + 1, j);
          stk.emplace(i, k);
          break;
        }
      }
    }
  }

  std::string result(seqLen, '.');
  for (const auto& bp : basePair) {
    result.at(bp.first) = '(';
    result.at(bp.second) = ')';
  }

  return result;
}

/*
 * Calculate Sensitivity
 * @param answer: answer secondary structure
 * @return: Sensitivity
 */
double Nussinov::calcSensitivity(const std::string& answer) const {
  const auto [TP, _, FN] = countTPFPFN(answer);
  return static_cast<double>(TP) / (TP + FN);
}

double Nussinov::calcPPV(const std::string& answer) const {
  const auto [TP, FP, _] = countTPFPFN(answer);
  return static_cast<double>(TP) / (TP + FP);
}

const std::string& Nussinov::getSeq() const { return seq_; }
