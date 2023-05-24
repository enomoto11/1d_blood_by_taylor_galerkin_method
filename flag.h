#include <cstdio>

struct Flag
{
  int index;
  bool shouldCalculate = true;
};

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
        cout << "Invalid index provided: " << index << ". Exiting program..." << endl;
        exit(1);
      }
    }
  }
};

class CheckArgs
{
private:
  double b_area_ele0;
  double b_area_ele1;
  double b_flowQuantity_ele0;
  double b_flowQuantity_ele1;
  int i;
  int j;
  int index;

public:
  static CheckArgs newCheckArgs(
      double b_area_ele0,
      double b_area_ele1,
      double b_flowQuantity_ele0,
      double b_flowQuantity_ele1,
      int i,
      int j,
      int index)
  {
    CheckArgs args;
    args.b_area_ele0 = b_area_ele0;
    args.b_area_ele1 = b_area_ele1;
    args.b_flowQuantity_ele0 = b_flowQuantity_ele0;
    args.b_flowQuantity_ele1 = b_flowQuantity_ele1;
    args.i = i;
    args.j = j;
    args.index = index;

    return args;
  }

  double getB_area_ele0() const { return b_area_ele0; }
  double getB_area_ele1() const { return b_area_ele1; }
  double getB_flowQuantity_ele0() const { return b_flowQuantity_ele0; }
  double getB_flowQuantity_ele1() const { return b_flowQuantity_ele1; }
  int getI() const { return i; }
  int getJ() const { return j; }
  int getIndex() const { return index; }
};

class ExceptionManager
{
public:
  static void check(const CheckArgs &args)
  {
    cout << "i = " << args.getI() << endl;
    cout << "j = " << args.getJ() << endl;
    cout << args.getIndex() << "th term is nan Exit..." << endl;
    exit(1);
  }
};
