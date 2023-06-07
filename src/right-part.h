#ifndef _RIGHT_PART_H_
#define _RIGHT_PART_H_

#include "flag.h"
// #include <iostream>

struct RightPart
{
  Flag firstTerm;
  Flag secondTerm;
  Flag thirdTerm;
  Flag fourthTerm;
  Flag fifthTerm;

  static RightPart newRightPart()
  {
    RightPart rp;
    rp.firstTerm.index = 1;
    rp.secondTerm.index = 2;
    rp.thirdTerm.index = 3;
    rp.fourthTerm.index = 4;
    rp.fifthTerm.index = 5;

    return rp;
  }

  static void disable(RightPart &rp, const std::vector<int> &indices)
  {
    if (indices.empty())
      return;

    for (int index : indices)
    {
      switch (index)
      {
      case 1:
        rp.firstTerm.shouldCalculate = false;
        break;
      case 2:
        rp.secondTerm.shouldCalculate = false;
        break;
      case 3:
        rp.thirdTerm.shouldCalculate = false;
        break;
      case 4:
        rp.fourthTerm.shouldCalculate = false;
        break;
      case 5:
        rp.fifthTerm.shouldCalculate = false;
        break;
      default:
        std::cout << "Invalid index provided: " << index << ". Exiting program..." << std::endl;
        exit(1);
      }
    }
  }
};

#endif