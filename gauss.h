class Gauss
{
public:
	Gauss() {}
	Gauss(const int num)
	{
		switch (num)
		{
		case 1:
			point[0] = -0.577350296189626;
			point[1] = 0.577350296189626;
			weight[0] = 1e0;
			weight[1] = 1e0;
			break;
		case 2:
			point[0] = -0.774596669241483;
			point[1] = 0e0;
			point[2] = 0.774596669241483;
			weight[0] = 0.555555555555555;
			weight[1] = 0.888888888888888;
			weight[2] = 0.555555555555555;
			break;
		case 3:
			point[0] = -0.861135311594053;
			point[1] = -0.339981043584856;
			point[2] = 0.339981043584856;
			point[3] = 0.861135311594053;
			weight[0] = 0.347854845137454;
			weight[1] = 0.652145154862546;
			weight[2] = 0.652145154862546;
			weight[3] = 0.347854845137454;
			break;
		default:
			printf("undefined order is set in the gauss integral\n");
		}
	}
	double point[4], weight[4];
};
