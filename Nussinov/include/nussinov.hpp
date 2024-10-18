#ifndef NUSSINOV_HPP_
#define NUSSINOV_HPP_

#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

class Nussinov {
 private:
  static constexpr int64_t kInf = 1e9;
  std::string seq_;
  std::vector<std::vector<int64_t>> dp_;

  int64_t delta(int i, int j) const;

  std::unordered_map<int, int> getBracketPair(const std::string& s) const;
  // TP, FP, FN
  std::tuple<int, int, int> countTPFPFN(const std::string& pred) const;

 public:
  Nussinov(std::string seq);
  ~Nussinov();

  void dp();
  std::string traceback() const;

  double calcSensitivity(const std::string& answer) const;
  double calcPPV(const std::string& answer) const;

  const std::string& getSeq() const;
};

#endif  // NUSSINOV_HPP_