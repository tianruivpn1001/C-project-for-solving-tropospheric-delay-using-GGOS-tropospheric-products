#include "rtklib.h"
#include <stdio.h>
#include <string.h>
#include <iostream>

/*********测试用主程序(注意：采用的是VMFG_FC数据)********/
int main()
{
	gtime_t time_test = transform_time("2020/1/1/3:0:0");//测试用历元
	gtime_t time_former_test= transform_time("2020/1/1/0:0:0");
	//下面用的相对路径，一定要保证数据文件位置正确！
	char *infile[] = { ".\\orography_ell.txt",".\\VMFG_20200101.H00.txt",".\\VMFG_20200101.H06.txt" };
	const double azel_test[2] = { 0,90 * D2R };//第一个是方位角，第二个是仰角
	const double pos_BJFS[3] = { (39 + 36.0 / 60.0 + 31.0 / 3600.0)*D2R,(115 + 53.0 / 60.0 + 33.0 / 3600.0)*D2R,87.5 };//用BJFS站进行测试（下面坐标取自下载的官方SNX文件）
	double result_test=0.0;//测试结果
	Demo_tropcorr(time_test, time_former_test,infile, pos_BJFS, azel_test, &result_test);
	printf("测试：2020年1月1日3时BJFS站的天顶对流层延迟解算结果为：%f m,此时的精确天顶对流层延迟为2.3606m（依据IGS发布的全球跟踪站数据）",result_test);
}

extern gtime_t transform_time(const char *time) {
	/*自己写的转换函数，与rtklib的epoch2time共同使用*/
	double result_d[6];//存储年月日时分秒六大元素
	gtime_t result_time;//最终返回结果

	//sscanf函数用法注意，非常好用的字符串分解函数
	sscanf(time, "%lf/%lf/%lf/%lf:%lf:%lf", result_d, result_d + 1, result_d + 2, result_d + 3, result_d + 4, result_d + 5);
	result_time = epoch2time(result_d);//epoch2time是rtklib中已经写好的函数，用于将六元素时间数组转化为rtklib中用的时间，rtklib用的时间是1970年至今的积秒
	return result_time;
}

