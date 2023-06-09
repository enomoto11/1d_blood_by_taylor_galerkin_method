#ifndef _EXCEPTION_ARG_H_
#define _EXCEPTION_ARG_H_

class ExceptionArgs
{
private:
  double b_area_ele0;
  double b_area_ele1;
  double b_flowQuantity_ele0;
  double b_flowQuantity_ele1;
  double A;
  double Q;
  int i;
  int j;
  int index;

public:
  static ExceptionArgs newCheckArgs(
      double b_area_ele0,
      double b_area_ele1,
      double b_flowQuantity_ele0,
      double b_flowQuantity_ele1,
      double A,
      double Q,
      int i,
      int j,
      int index)
  {
    ExceptionArgs args;
    args.b_area_ele0 = b_area_ele0;
    args.b_area_ele1 = b_area_ele1;
    args.b_flowQuantity_ele0 = b_flowQuantity_ele0;
    args.b_flowQuantity_ele1 = b_flowQuantity_ele1;
    args.A = A;
    args.Q = Q;
    args.i = i;
    args.j = j;
    args.index = index;

    return args;
  }

  double getB_area_ele0() const { return b_area_ele0; }
  double getB_area_ele1() const { return b_area_ele1; }
  double getB_flowQuantity_ele0() const { return b_flowQuantity_ele0; }
  double getB_flowQuantity_ele1() const { return b_flowQuantity_ele1; }
  double getA() const { return A; }
  double getQ() const { return Q; }
  int getI() const { return i; }
  int getJ() const { return j; }
  int getIndex() const { return index; }
};

#endif