#ifndef ZUKKER_HPP_
#define ZUKKER_HPP_

#include <cassert>
#include <stack>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

enum LoopType {
  HAIRPIN,
  STACKING,
  INTERNAL,
  BULGE,
  MULTILOOP,
};

enum BasePair {
  GC,
  AU,
  GU,
  UNMATCH,
};

namespace utils {
/*
 * utils::chmin: Update a value if the new value is smaller
 * @param a: value to be updated
 * @param b: new value
 */
template <typename T>
static bool chmin(T& a, const T& b) {
  if (b < a) {
    a = b;
    return true;
  }
  return false;
}
}  // namespace utils

class Zuker {
 private:
  static constexpr int64_t kInf = 1e18;
  // static constexpr double kInitialValue = 42;
  static constexpr double kEps = 1e-10;

  static constexpr double a = 5.0;
  static constexpr double b = -1.0;
  static constexpr double c = 0.1;

  std::string seq_;

  std::vector<std::vector<double>> W_;
  std::vector<std::vector<double>> V_;
  std::vector<std::vector<double>> M_;

  LoopType classifyLoop(int i, int j, int h, int l) const;

  double getLoopData(int i, int j) const;
  BasePair classifyBases(int i, int j) const;
  double getStackLoopData(BasePair a, BasePair b) const;
  double calcStackingEnergy(int i, int j) const;

  double F1(int i, int j) const;
  double F2(int i, int j, int h, int l) const;

  void updateW(int i, int j);
  void updateV(int i, int j);
  void updateM(int i, int j);

  std::unordered_map<int, int> getBracketPair(const std::string& s) const;
  // TP, FP, FN
  std::tuple<int, int, int> countTPFPFN(const std::string& pred) const;

 public:
  Zuker(const std::string& seq);

  void dp();
  std::string traceback() const;

  const std::string& getSeq() const;

  double calcSensitivity(const std::string& answer) const;
  double calcPPV(const std::string& answer) const;
};

#endif  // ZUKKER_HPP_