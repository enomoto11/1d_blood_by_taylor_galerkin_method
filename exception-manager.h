class ExceptionManager
{
public:
  static void check(const ExceptionArgs &args)
  {
    std::cout << "i = " << args.getI() << std::endl;
    std::cout << "j = " << args.getJ() << std::endl;
    std::cout << args.getIndex() << "th term is nan Exit..." << std::endl;
    exit(1);
  }
};
