#include "flag.h"
struct LeftPart
{
  Flag firstTerm;
  Flag secondTerm;
  Flag thirdTerm;
  Flag fourthTerm;
  Flag fifthTerm;

  static LeftPart newLeftPart()
  {
    LeftPart lp;
    lp.firstTerm.index = 1;
    lp.secondTerm.index = 2;
    lp.thirdTerm.index = 3;
    lp.fourthTerm.index = 4;
    lp.fifthTerm.index = 5;

    return lp;
  }

  static void disable(LeftPart &lp, const std::vector<int> &indices)
  {
    if (indices.empty())
      return;

    for (int index : indices)
    {
      switch (index)
      {
      case 1:
        lp.firstTerm.shouldCalculate = false;
        break;
      case 2:
        lp.secondTerm.shouldCalculate = false;
        break;
      case 3:
        lp.thirdTerm.shouldCalculate = false;
        break;
      case 4:
        lp.fourthTerm.shouldCalculate = false;
        break;
      case 5:
        lp.fifthTerm.shouldCalculate = false;
        break;
      default:
        std::cout << "Invalid index provided: " << index << ". Exiting program..." << std::endl;
        exit(1);
      }
    }
  }
};