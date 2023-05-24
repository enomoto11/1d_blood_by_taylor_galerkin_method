struct LeftPartFlag
{
  bool shouldCalculateFirstTerm = true;
  bool shouldCalculateSecondTerm = true;
  bool shouldCalculateThirdTerm = true;
  bool shouldCalculateFourthTerm = true;
  bool shouldCalculateFifthTerm = true;
};

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
};

struct CheckArgs
{
  double b_area_ele0;
  double b_area_ele1;
  double b_flowQuantity_ele0;
  double b_flowQuantity_ele1;
  int i;
  int j;
  int index;

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
};

struct ExceptionManager
{
  static void check(CheckArgs args)
  {
    cout << "i = " << args.i << endl;
    cout << "j = " << args.j << endl;
    cout << args.index << "th term is nan Exit..." << endl;
    exit(1);
  }
};