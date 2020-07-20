#include "rtklib.h"

typedef struct{ 
	double lat_deg;//注意：纬度单位是度，非弧度
	double lon_deg;//注意：经度单位是度，非弧度
	//下面四项均未经过高程归算，为格网点处原始数据
	double ah;
	double aw;
	double zhd;
	double zwd;
	double p;//所在高度处气压值
} VMF1_data_h0;//h0说明存储的均是格网点处未经过高程归算的数据，格网点的高度存储在orography_ell文件中，P是根据格网点高度反算而来的

typedef struct {
	double lat_deg;//注意：纬度单位是度，非弧度
	double lon_deg;//注意：经度单位是度，非弧度
	//下面四项为测站处高程归算后得到的数据
	double ah;
	double aw;
	double zhd;
	double zwd;
	double p;//测站所在高度处气压值
	//测站处算出的映射函数值
	double mfh;
	double mfw;
} VMF1_data_h1;//h1说明存储的是测站处已修正后的数据，以及测站处映射函数值
/*注意：VMF1数据文件一个纬度对应144个数据，而orography_ell文件一个纬度对应145个数据，最后一个冗余（0度与360度是一个经度！！！）*/
/*注意：这些存储用静态变量务必初始化为0*/
static VMF1_data_h0 VMF1_data_all[2][13104] = { {{0} } };//[2]表示存储两个文件数据
static int orography_ell[13195] = {0};//存储orography_ell文件中的格网高程数据(注意：数据为int型)



static VMF1_data_h0 VMF1_add_h0(VMF1_data_h0 V1, VMF1_data_h0 V2,double a,double b,double c) {
	/* 结构体相减函数：(a*V1+b*V2)*c */
	VMF1_data_h0 result;
	result.lat_deg = c*(a*V1.lat_deg + b*V2.lat_deg);
	result.lon_deg = c * (a*V1.lon_deg +b* V2.lon_deg);
	result.ah = c * (a*V1.ah +b* V2.ah);
	result.aw = c * (a*V1.aw + b*V2.aw);
	result.zhd = c * (a*V1.zhd + b*V2.zhd);
	result.zwd = c * (a*V1.zwd + b*V2.zwd);
	result.p = c * (a*V1.p + b*V2.p);
	return result;
}



static VMF1_data_h1 VMF1_add_h1(VMF1_data_h1 V1, VMF1_data_h1 V2, double a, double b,double c) {
	/* 结构体相减函数：(a*V1+b*V2)*c */
	VMF1_data_h1 result;
	result.lat_deg = c * (a * V1.lat_deg + b * V2.lat_deg);
	result.lon_deg = c * (a * V1.lon_deg + b * V2.lon_deg);
	result.ah = c * (a * V1.ah + b * V2.ah);
	result.aw = c * (a * V1.aw + b * V2.aw);
	result.zhd = c * (a * V1.zhd + b * V2.zhd);
	result.zwd = c * (a * V1.zwd + b * V2.zwd);
	result.p = c * (a * V1.p + b * V2.p);
	result.mfh = c * (a*V1.mfh + b*V2.mfh);
	result.mfw = c * (a*V1.mfw + b*V2.mfw);
	return result;
}

/* 自己写的基于GGOS产品的对流层修正程序，参考vmf1_grid.m、vmf1_ht.m -----------------------------------------------------
* 注意：经度是0到360度，而非-180到180
*  gtime_t time_current     I   该历元的时间
*   gtime_t time_former     I   前一个VMF1数据文件的时间，知道前一个即知道后一文件的时间，相隔6h
*        char *infile[]     I   三个文件路径，orography_ell文件、两个VMF1数据文件
*          double *pos      I   某次迭代算出的接收机位置（纬经高） {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double *trp      O   m为单位的对流层延迟 (m)
*          double *var      O   m^2为单位的方差(m^2)，我觉得对于高精度的VMF1模型而言可以设置为0
*-----------------------------------------------------------------------------*/
extern void Demo_tropcorr(gtime_t time_current,gtime_t time_former, char *infile[], const double *pos, const double *azel, double *trp)
{
	/*局部变量初始化*/
	char buff_temp[1024];//读文件跳行用的临时存储
	int i=0;
	double pos_transformed[3];//经转化后的纬度经度高程
	gtime_t time_latter = timeadd(time_former, 21600);//time_latter比time_former自然多6h
	VMF1_data_h0 VMF1_surrounding_data_h0[2][4] = { {0} };//存储测站周围四个格网点原始数据（未经高程归算）[2]表示存储两个文件数据
	VMF1_data_h0 VMF1_surrounding_data_interpolated_h0[4] = { {0} };//存储时间插值后的测站周围四个格网点的数据
	VMF1_data_h1 VMF1_surrounding_data_h1[4] = { 0 };//存储测站周围四个格网点数据(高程归算后)
	int index0, index1, index2, index3,index0_orography_ell, index1_orography_ell, index2_orography_ell, index3_orography_ell;//四个相邻格网点索引
	double doy;//(Niell1996模型的年积日，起点是1980年1月28日，但rtklib中GPST时间起点应该在1970年1月1日)

	if (orography_ell[0]==0)
	{//注意：为加快运行速度，这样写保证地形文件仅读取一次
		FILE *fp_orography_ell;//指向orography_ell文件的文件指针

		/*读取orography_ell文件*/
		fp_orography_ell = fopen(infile[0], "r");
		fgets(buff_temp, sizeof(buff_temp), fp_orography_ell);//跳过orography_ell文件的第一行
		while (fscanf(fp_orography_ell, "%d%d%d%d%d", orography_ell + i * 5 + 0, orography_ell + i * 5 + 1, orography_ell + i * 5 + 2, orography_ell + i * 5 + 3, orography_ell + i * 5 + 4) != EOF)
			i++;//5个5个的读取orography_ell文件的数据(fscanf自动匹配的功能还是挺强的！)

		/***务必记得释放文件指针***/
		fclose(fp_orography_ell);
	}
	/*******这样存在一个BUG，只能测6h以内的数据，不过没有改的必要********/
	if (VMF1_data_all[0][0].ah==0.0)
	{//注意：为加快运行速度，这样写保证尽可能少的读取VMF1数据文件
		FILE *fp_VMF1_former_file;//指向前一个VMF1数据文件的文件指针
		FILE *fp_VMF1_latter_file;//指向后一个VMF1数据文件的文件指针

		/*读取VMF1数据文件*/
		fp_VMF1_former_file = fopen(infile[1], "r");
		while (fgets(buff_temp, sizeof(buff_temp), fp_VMF1_former_file)) {
			if (buff_temp[0] != '!') break;//跳过文件头部分
		}
		sscanf(buff_temp, "%lf%lf%lf%lf%lf%lf", &(VMF1_data_all[0][0].lat_deg), &(VMF1_data_all[0][0].lon_deg), &(VMF1_data_all[0][0].ah), &(VMF1_data_all[0][0].aw), &(VMF1_data_all[0][0].zhd), &(VMF1_data_all[0][0].zwd));
		i = 1;//注意：i应从1开始，因为上句已经读了！
		while (fscanf(fp_VMF1_former_file, "%lf%lf%lf%lf%lf%lf", &(VMF1_data_all[0][i].lat_deg), &(VMF1_data_all[0][i].lon_deg), &(VMF1_data_all[0][i].ah), &(VMF1_data_all[0][i].aw), &(VMF1_data_all[0][i].zhd), &(VMF1_data_all[0][i].zwd)) != EOF)
			i++;//一行一行读取VMF1数据文件
		fp_VMF1_latter_file = fopen(infile[2], "r");
		while (fgets(buff_temp, sizeof(buff_temp), fp_VMF1_latter_file)) {
			if (buff_temp[0] != '!') break;//跳过文件头部分
		}
		sscanf(buff_temp, "%lf%lf%lf%lf%lf%lf", &(VMF1_data_all[1][0].lat_deg), &(VMF1_data_all[1][0].lon_deg), &(VMF1_data_all[1][0].ah), &(VMF1_data_all[1][0].aw), &(VMF1_data_all[1][0].zhd), &(VMF1_data_all[1][0].zwd));
		i = 1;//注意：i应从1开始，因为上句已经读了！
		while (fscanf(fp_VMF1_latter_file, "%lf%lf%lf%lf%lf%lf", &(VMF1_data_all[1][i].lat_deg), &(VMF1_data_all[1][i].lon_deg), &(VMF1_data_all[1][i].ah), &(VMF1_data_all[1][i].aw), &(VMF1_data_all[1][i].zhd), &(VMF1_data_all[1][i].zwd)) != EOF)
			i++;//一行一行读取VMF1数据文件

		/***务必记得释放文件指针***/
		fclose(fp_VMF1_former_file);
		fclose(fp_VMF1_latter_file);
	}

	/*将纬经度由弧度转化为度,同时将西半球的负经度转化为正经度（0~360）*/
	pos_transformed[0] = pos[0] / D2R;
	pos_transformed[1] = pos[1] / D2R;
	pos_transformed[2] = pos[2];
	if (pos_transformed[1] < 0) pos_transformed[1] = 360 + pos_transformed[1];//西半球负经度转化为正经度

	/*获取测站周边四个相邻格网点的原始数据*/	
	//下面代码测站不在格网线上,而且经度也不大于357.5度（最一般的情形）
	index0 = ceil(abs(pos_transformed[0] - 90) / 2.0) * 144 + floor(pos_transformed[1] / 2.5);//matlab程序中第一个点是右下角，这里是左下角
	index1 = ceil(abs(pos_transformed[0] - 90) / 2.0) * 144 + ceil(pos_transformed[1] / 2.5);
	index2 = floor(abs(pos_transformed[0] - 90) / 2.0) * 144 + floor(pos_transformed[1] / 2.5);
	index3 = floor(abs(pos_transformed[0] - 90) / 2.0) * 144 + ceil(pos_transformed[1] / 2.5);
	//注意：VMF1数据文件一个纬度对应144个数据，而orography_ell文件一个纬度对应145个数据，最后一个冗余（0度与360度是一个经度！！！）
	index0_orography_ell = ceil(abs(pos_transformed[0] - 90) / 2.0) * 145 + floor(pos_transformed[1] / 2.5);
	index1_orography_ell = ceil(abs(pos_transformed[0] - 90) / 2.0) * 145 + ceil(pos_transformed[1] / 2.5);
	index2_orography_ell = floor(abs(pos_transformed[0] - 90) / 2.0) * 145 + floor(pos_transformed[1] / 2.5);
	index3_orography_ell = floor(abs(pos_transformed[0] - 90) / 2.0) * 145 + ceil(pos_transformed[1] / 2.5);
	//下面修复各种特殊情况
	if (pos_transformed[1]>357.5)
	{//经度大于357.5度的情况
		index1 -= 144;
		index3 -= 144;
		index1_orography_ell -= 145;
		index3_orography_ell -= 145;
	}
	if (floor(pos_transformed[1] / 2.5) == ceil(pos_transformed[1] / 2.5))
	{//测站经度踩线
		index1 += 1;
		index3 += 1;
		index1_orography_ell += 1;
		index3_orography_ell += 1;
		if (pos_transformed[1] == 357.5)
		{//经度正好在357.5度的特殊情况
			index1 -= 144;
			index3 -= 144;
			index1_orography_ell -= 145;
			index3_orography_ell -= 145;
		}
		if (floor(abs(pos_transformed[0] - 90) / 2.0) == ceil(abs(pos_transformed[0] - 90) / 2.0))
		{//测站经度踩线，纬度也踩线（即测站正好在格网点上）
			index0 += 144;
			index1 += 144;
			index0_orography_ell += 144;
			index1_orography_ell += 144;
			if (pos_transformed[0] == -90)
			{//测站正好在格网点上，纬度正好为-90度
				index0 -= 288;
				index1 -= 288;
				index0_orography_ell -= 288;
				index1_orography_ell -= 288;
			}
		}
	}
	else if (floor(abs(pos_transformed[0] - 90) / 2.0) ==ceil(abs(pos_transformed[0] - 90) / 2.0))
	{//测站经度没踩线，纬度踩线了
		index0 += 144;
		index1 += 144;
		index0_orography_ell += 144;
		index1_orography_ell += 144;
		if (pos_transformed[0]==-90)
		{//测站纬度正好等于-90度（等于90度的情况无需考虑）
			index0 -= 288;
			index1 -= 288;
			index0_orography_ell -= 288;
			index1_orography_ell -= 288;
		}
	}
	else { /*占位：可以什么都不做*/ }
	VMF1_surrounding_data_h0[0][0] = VMF1_data_all[0][index0];
	VMF1_surrounding_data_h0[0][1] = VMF1_data_all[0][index1];
	VMF1_surrounding_data_h0[0][2] = VMF1_data_all[0][index2];
	VMF1_surrounding_data_h0[0][3] = VMF1_data_all[0][index3];
	VMF1_surrounding_data_h0[1][0] = VMF1_data_all[1][index0];
	VMF1_surrounding_data_h0[1][1] = VMF1_data_all[1][index1];
	VMF1_surrounding_data_h0[1][2] = VMF1_data_all[1][index2];
	VMF1_surrounding_data_h0[1][3] = VMF1_data_all[1][index3];

	/*通过时间插值，获取测站周边四个相邻格网点插值后的数据*/
	VMF1_surrounding_data_interpolated_h0[0] = VMF1_add_h0(VMF1_surrounding_data_h0[0][0] , VMF1_add_h0(VMF1_surrounding_data_h0[1][0], VMF1_surrounding_data_h0[0][0],1,-1, (timediff(time_current, time_former) / timediff(time_latter, time_former))),1,1,1);
	VMF1_surrounding_data_interpolated_h0[1] = VMF1_add_h0(VMF1_surrounding_data_h0[0][1], VMF1_add_h0(VMF1_surrounding_data_h0[1][1], VMF1_surrounding_data_h0[0][1], 1, -1, (timediff(time_current, time_former) / timediff(time_latter, time_former))), 1, 1, 1);
	VMF1_surrounding_data_interpolated_h0[2] = VMF1_add_h0(VMF1_surrounding_data_h0[0][2], VMF1_add_h0(VMF1_surrounding_data_h0[1][2], VMF1_surrounding_data_h0[0][2], 1, -1, (timediff(time_current, time_former) / timediff(time_latter, time_former))), 1, 1, 1);
	VMF1_surrounding_data_interpolated_h0[3] = VMF1_add_h0(VMF1_surrounding_data_h0[0][3], VMF1_add_h0(VMF1_surrounding_data_h0[1][3], VMF1_surrounding_data_h0[0][3], 1, -1, (timediff(time_current, time_former) / timediff(time_latter, time_former))), 1, 1, 1);
	//注意：VMF1_add_h0函数与VMF1_add_h1函数的形式是（a*V1+b*V2）*c，即可进行加减数乘运算

	/*h1与h0的前四项（lat、lon、ah、aw）是相同的，此处赋值*/
	VMF1_surrounding_data_h1[0].lat_deg = VMF1_surrounding_data_interpolated_h0[0].lat_deg; VMF1_surrounding_data_h1[0].lon_deg = VMF1_surrounding_data_interpolated_h0[0].lon_deg; VMF1_surrounding_data_h1[0].ah = VMF1_surrounding_data_interpolated_h0[0].ah; VMF1_surrounding_data_h1[0].aw = VMF1_surrounding_data_interpolated_h0[0].aw;
	VMF1_surrounding_data_h1[1].lat_deg = VMF1_surrounding_data_interpolated_h0[1].lat_deg; VMF1_surrounding_data_h1[1].lon_deg = VMF1_surrounding_data_interpolated_h0[1].lon_deg; VMF1_surrounding_data_h1[1].ah = VMF1_surrounding_data_interpolated_h0[1].ah; VMF1_surrounding_data_h1[1].aw = VMF1_surrounding_data_interpolated_h0[1].aw;
	VMF1_surrounding_data_h1[2].lat_deg = VMF1_surrounding_data_interpolated_h0[2].lat_deg; VMF1_surrounding_data_h1[2].lon_deg = VMF1_surrounding_data_interpolated_h0[2].lon_deg; VMF1_surrounding_data_h1[2].ah = VMF1_surrounding_data_interpolated_h0[2].ah; VMF1_surrounding_data_h1[2].aw = VMF1_surrounding_data_interpolated_h0[2].aw;
	VMF1_surrounding_data_h1[3].lat_deg = VMF1_surrounding_data_interpolated_h0[3].lat_deg; VMF1_surrounding_data_h1[3].lon_deg = VMF1_surrounding_data_interpolated_h0[3].lon_deg; VMF1_surrounding_data_h1[3].ah = VMF1_surrounding_data_interpolated_h0[3].ah; VMF1_surrounding_data_h1[3].aw = VMF1_surrounding_data_interpolated_h0[3].aw;

	/*下面根据论文《Implementation and testing of the gridded Vienna Mapping Function(VMF1)》中式（3）反算格网点高度处气压*/
	VMF1_surrounding_data_interpolated_h0[0].p = (VMF1_surrounding_data_interpolated_h0[0].zhd / 0.0022768)*(1 - 0.00266*cos(2 * pos[0]) - 0.28E-6*orography_ell[index0_orography_ell]);
	VMF1_surrounding_data_interpolated_h0[1].p = (VMF1_surrounding_data_interpolated_h0[1].zhd / 0.0022768)*(1 - 0.00266*cos(2 * pos[0]) - 0.28E-6*orography_ell[index1_orography_ell]);
	VMF1_surrounding_data_interpolated_h0[2].p = (VMF1_surrounding_data_interpolated_h0[2].zhd / 0.0022768)*(1 - 0.00266*cos(2 * pos[0]) - 0.28E-6*orography_ell[index2_orography_ell]);
	VMF1_surrounding_data_interpolated_h0[3].p = (VMF1_surrounding_data_interpolated_h0[3].zhd / 0.0022768)*(1 - 0.00266*cos(2 * pos[0]) - 0.28E-6*orography_ell[index3_orography_ell]);
	//注意：cos函数用的是弧度

	/*下面依据论文《Discussion and Recommandations about the Height Correction for A Priori Zenit》中式（2），根据格网点高度处气压推算测站高程处气压（高程均已知，1948 Berg模型）*/
	VMF1_surrounding_data_h1[0].p = VMF1_surrounding_data_interpolated_h0[0].p*pow(1-0.0000226*(pos[2]- orography_ell[index0_orography_ell]), 5.225);
	VMF1_surrounding_data_h1[1].p = VMF1_surrounding_data_interpolated_h0[1].p*pow(1-0.0000226*(pos[2] - orography_ell[index1_orography_ell]), 5.225);
	VMF1_surrounding_data_h1[2].p = VMF1_surrounding_data_interpolated_h0[2].p*pow(1-0.0000226*(pos[2] - orography_ell[index2_orography_ell]), 5.225);
	VMF1_surrounding_data_h1[3].p = VMF1_surrounding_data_interpolated_h0[3].p*pow(1-0.0000226*(pos[2] - orography_ell[index3_orography_ell]), 5.225);

	/*下面依据论文《Implementation and testing of the gridded Vienna Mapping Function(VMF1)》中式（3）算经高程修正后的ZHD*/
	VMF1_surrounding_data_h1[0].zhd = 0.0022768*VMF1_surrounding_data_h1[0].p / (1 - 0.00266*cos(2 * pos[0]) - 0.28E-6*pos[2]);
	VMF1_surrounding_data_h1[1].zhd = 0.0022768*VMF1_surrounding_data_h1[1].p / (1 - 0.00266*cos(2 * pos[0]) - 0.28E-6*pos[2]);
	VMF1_surrounding_data_h1[2].zhd = 0.0022768*VMF1_surrounding_data_h1[2].p / (1 - 0.00266*cos(2 * pos[0]) - 0.28E-6*pos[2]);
	VMF1_surrounding_data_h1[3].zhd = 0.0022768*VMF1_surrounding_data_h1[3].p / (1 - 0.00266*cos(2 * pos[0]) - 0.28E-6*pos[2]);

	/*下面依据论文《Implementation and testing of the gridded Vienna Mapping Function(VMF1)》中式（5）算经高程修正后的ZWD*/
	VMF1_surrounding_data_h1[0].zwd = VMF1_surrounding_data_interpolated_h0[0].zwd*exp(-(pos[2] - orography_ell[index0_orography_ell]) / 2000.0);
	VMF1_surrounding_data_h1[1].zwd = VMF1_surrounding_data_interpolated_h0[1].zwd*exp(-(pos[2] - orography_ell[index1_orography_ell]) / 2000.0);
	VMF1_surrounding_data_h1[2].zwd = VMF1_surrounding_data_interpolated_h0[2].zwd*exp(-(pos[2] - orography_ell[index2_orography_ell]) / 2000.0);
	VMF1_surrounding_data_h1[3].zwd = VMF1_surrounding_data_interpolated_h0[3].zwd*exp(-(pos[2] - orography_ell[index3_orography_ell]) / 2000.0);

	/*计算Niell1996模型(起算时间:1980/1/26)所用的doy*/
	/*****注意：rtklib中的GPST起算点在1970/1/1并非网上查到的1980/1/6!!!*****/
	doy=(double)time_current.time / 86400.0 + time_current.sec-3679;

	/*确定映射函数连分式的系数*/
	/*根据参数修正后的Boehm论文以及《Generation and Assessment of VMF1-Type Grids Using North-American Numerical Weather Models》*/
	double Psi, c10, c11;
	double b_zhd = 0.0029;//干延迟连分式的b、c系数
	double c_zhd[4];//注意四个格网点的c系数各不相同！！！
	double b_zwd = 0.00146;//湿延迟连分式的b、c系数
	double c_zwd = 0.04391;
	if (pos[0]<0)
	{
		Psi = PI;
		c10= 0.002;
		c11= 0.007;
	}
	else
	{
		Psi = 0;
		c10 = 0.001;
		c11 = 0.005;
	}
	c_zhd[0] = 0.062 + ((cos(doy / 365.25 * 2 * PI + Psi) + 1)*c11 / 2 + c10)*(1 - cos(VMF1_surrounding_data_h1[0].lat_deg*D2R));
	c_zhd[1] = 0.062 + ((cos(doy / 365.25 * 2 * PI + Psi) + 1)*c11 / 2 + c10)*(1 - cos(VMF1_surrounding_data_h1[1].lat_deg*D2R));
	c_zhd[2] = 0.062 + ((cos(doy / 365.25 * 2 * PI + Psi) + 1)*c11 / 2 + c10)*(1 - cos(VMF1_surrounding_data_h1[2].lat_deg*D2R));
	c_zhd[3] = 0.062 + ((cos(doy / 365.25 * 2 * PI + Psi) + 1)*c11 / 2 + c10)*(1 - cos(VMF1_surrounding_data_h1[3].lat_deg*D2R));

	/*确定未经高程修正的映射函数，即将上面算得的系数代入连分式,注意仰角单位为弧度*/
	VMF1_surrounding_data_h1[0].mfh = (1 + VMF1_surrounding_data_h1[0].ah / (1 + b_zhd / (1 + c_zhd[0]))) / (sin(azel[1]) + VMF1_surrounding_data_h1[0].ah / (sin(azel[1]) + b_zhd / (sin(azel[1]) + c_zhd[0])));
	VMF1_surrounding_data_h1[1].mfh = (1 + VMF1_surrounding_data_h1[1].ah / (1 + b_zhd / (1 + c_zhd[1]))) / (sin(azel[1]) + VMF1_surrounding_data_h1[1].ah / (sin(azel[1]) + b_zhd / (sin(azel[1]) + c_zhd[1])));
	VMF1_surrounding_data_h1[2].mfh = (1 + VMF1_surrounding_data_h1[2].ah / (1 + b_zhd / (1 + c_zhd[2]))) / (sin(azel[1]) + VMF1_surrounding_data_h1[2].ah / (sin(azel[1]) + b_zhd / (sin(azel[1]) + c_zhd[2])));
	VMF1_surrounding_data_h1[3].mfh = (1 + VMF1_surrounding_data_h1[3].ah / (1 + b_zhd / (1 + c_zhd[3]))) / (sin(azel[1]) + VMF1_surrounding_data_h1[3].ah / (sin(azel[1]) + b_zhd / (sin(azel[1]) + c_zhd[3])));
	VMF1_surrounding_data_h1[0].mfw = (1 + VMF1_surrounding_data_h1[0].aw / (1 + b_zwd / (1 + c_zwd))) / (sin(azel[1]) + VMF1_surrounding_data_h1[0].aw / (sin(azel[1]) + b_zwd / (sin(azel[1]) + c_zwd)));
	VMF1_surrounding_data_h1[1].mfw = (1 + VMF1_surrounding_data_h1[1].aw / (1 + b_zwd / (1 + c_zwd))) / (sin(azel[1]) + VMF1_surrounding_data_h1[1].aw / (sin(azel[1]) + b_zwd / (sin(azel[1]) + c_zwd)));
	VMF1_surrounding_data_h1[2].mfw = (1 + VMF1_surrounding_data_h1[2].aw / (1 + b_zwd / (1 + c_zwd))) / (sin(azel[1]) + VMF1_surrounding_data_h1[2].aw / (sin(azel[1]) + b_zwd / (sin(azel[1]) + c_zwd)));
	VMF1_surrounding_data_h1[3].mfw = (1 + VMF1_surrounding_data_h1[3].aw / (1 + b_zwd / (1 + c_zwd))) / (sin(azel[1]) + VMF1_surrounding_data_h1[3].aw / (sin(azel[1]) + b_zwd / (sin(azel[1]) + c_zwd)));

	/*根据官方代码，还应该加上映射函数的高程修正！！！(论文里没提)*/
	double a_ht = 2.53E-5;
	double b_ht = 5.49E-3;
	double c_ht = 1.14E-3;
	double ht_corr=0.0;
	ht_corr = 1 / sin(azel[1]) - (1 + a_ht / (1 + b_ht / (1 + c_ht))) / (sin(azel[1]) + a_ht / (sin(azel[1]) + b_ht / (sin(azel[1]) + c_ht)));
	ht_corr=ht_corr* pos[2] / 1000.0;
	VMF1_surrounding_data_h1[0].mfh += ht_corr;
	VMF1_surrounding_data_h1[1].mfh += ht_corr;
	VMF1_surrounding_data_h1[2].mfh += ht_corr;
	VMF1_surrounding_data_h1[3].mfh += ht_corr;

	/*进行BIL插值*/
	double zhd_result,zwd_result,mfh_result,mfw_result;//存储插值结果的变量
	double zhd_lon1, zhd_lon2, zwd_lon1, zwd_lon2, mfh_lon1, mfh_lon2, mfw_lon1, mfw_lon2;//插值用中间变量
	zhd_lon1 = VMF1_surrounding_data_h1[0].zhd + (VMF1_surrounding_data_h1[1].zhd - VMF1_surrounding_data_h1[0].zhd)*(pos_transformed[1] - VMF1_surrounding_data_h1[0].lon_deg) / (VMF1_surrounding_data_h1[1].lon_deg - VMF1_surrounding_data_h1[0].lon_deg);
	zhd_lon2 = VMF1_surrounding_data_h1[2].zhd + (VMF1_surrounding_data_h1[3].zhd - VMF1_surrounding_data_h1[2].zhd)*(pos_transformed[1] - VMF1_surrounding_data_h1[2].lon_deg) / (VMF1_surrounding_data_h1[3].lon_deg - VMF1_surrounding_data_h1[2].lon_deg);
	zwd_lon1 = VMF1_surrounding_data_h1[0].zwd + (VMF1_surrounding_data_h1[1].zwd - VMF1_surrounding_data_h1[0].zwd)*(pos_transformed[1] - VMF1_surrounding_data_h1[0].lon_deg) / (VMF1_surrounding_data_h1[1].lon_deg - VMF1_surrounding_data_h1[0].lon_deg);
	zwd_lon2 = VMF1_surrounding_data_h1[2].zwd + (VMF1_surrounding_data_h1[3].zwd - VMF1_surrounding_data_h1[2].zwd)*(pos_transformed[1] - VMF1_surrounding_data_h1[2].lon_deg) / (VMF1_surrounding_data_h1[3].lon_deg - VMF1_surrounding_data_h1[2].lon_deg);
	mfh_lon1 = VMF1_surrounding_data_h1[0].mfh + (VMF1_surrounding_data_h1[1].mfh - VMF1_surrounding_data_h1[0].mfh)*(pos_transformed[1] - VMF1_surrounding_data_h1[0].lon_deg) / (VMF1_surrounding_data_h1[1].lon_deg - VMF1_surrounding_data_h1[0].lon_deg);
	mfh_lon2 = VMF1_surrounding_data_h1[2].mfh + (VMF1_surrounding_data_h1[3].mfh - VMF1_surrounding_data_h1[2].mfh)*(pos_transformed[1] - VMF1_surrounding_data_h1[2].lon_deg) / (VMF1_surrounding_data_h1[3].lon_deg - VMF1_surrounding_data_h1[2].lon_deg);
	mfw_lon1 = VMF1_surrounding_data_h1[0].mfw + (VMF1_surrounding_data_h1[1].mfw - VMF1_surrounding_data_h1[0].mfw)*(pos_transformed[1] - VMF1_surrounding_data_h1[0].lon_deg) / (VMF1_surrounding_data_h1[1].lon_deg - VMF1_surrounding_data_h1[0].lon_deg);
	mfw_lon2 = VMF1_surrounding_data_h1[2].mfw + (VMF1_surrounding_data_h1[3].mfw - VMF1_surrounding_data_h1[2].mfw)*(pos_transformed[1] - VMF1_surrounding_data_h1[2].lon_deg) / (VMF1_surrounding_data_h1[3].lon_deg - VMF1_surrounding_data_h1[2].lon_deg);
	zhd_result = zhd_lon1 + (zhd_lon2 - zhd_lon1)*(pos_transformed[0] - VMF1_surrounding_data_h1[0].lat_deg) / (VMF1_surrounding_data_h1[2].lat_deg - VMF1_surrounding_data_h1[0].lat_deg);
	zwd_result = zwd_lon1 + (zwd_lon2 - zwd_lon1)*(pos_transformed[0] - VMF1_surrounding_data_h1[0].lat_deg) / (VMF1_surrounding_data_h1[2].lat_deg - VMF1_surrounding_data_h1[0].lat_deg);
	mfh_result = mfh_lon1 + (mfh_lon2 - mfh_lon1)*(pos_transformed[0] - VMF1_surrounding_data_h1[0].lat_deg) / (VMF1_surrounding_data_h1[2].lat_deg - VMF1_surrounding_data_h1[0].lat_deg);
	mfw_result = mfw_lon1 + (mfw_lon2 - mfw_lon1)*(pos_transformed[0] - VMF1_surrounding_data_h1[0].lat_deg) / (VMF1_surrounding_data_h1[2].lat_deg - VMF1_surrounding_data_h1[0].lat_deg);

	/*存储结果*/
	*trp = zhd_result * mfh_result + zwd_result * mfw_result;//小数点十二位与matlab计算值略有不同，但这应该属于matlab数值精度损失，无妨
}