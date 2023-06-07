#ifndef _EXCEPTION_MANAGER_H_
#define _EXCEPTION_MANAGER_H_

class ExceptionManager
{
public:
  static void check(const ExceptionArgs &args)
  {
    std::cout << "i = " << args.getI() << std::endl;
    std::cout << "j = " << args.getJ() << std::endl;
    std::cout << "b_area_ele0 = " << args.getB_area_ele0() << std::endl;
    std::cout << "b_area_ele1 = " << args.getB_area_ele1() << std::endl;
    std::cout << "b_flowQuantity_ele0 = " << args.getB_flowQuantity_ele0() << std::endl;
    std::cout << "b_flowQuantity_ele1 = " << args.getB_flowQuantity_ele1() << std::endl;
    std::cout << "A = " << args.getA() << std::endl;
    std::cout << "Q = " << args.getQ() << std::endl;
    std::cout << args.getIndex() << "th term is nan Exit..." << std::endl;
    exit(1);
  }
};

#endif