#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "rt_nonfinite.h"
#include "mean.h"
#include "diff.h"
#include "pqrstDetection_22nov2018_emxutil.h"
#include "fprintf.h"
#include "all.h"
#include "findpeaks.h"
#include "conv.h"
#include "power.h"
#include "filtfilt.h"
#include "abs.h"
#include "filter.h"
#include "rtGetNaN.h"
#include "rtGetInf.h"
#include "eml_setop.h"
#include "sort1.h"
#include "fileManager.h"
#include "sortIdx.h"
#include <stdbool.h>
#include <stdio.h>

//extern int errno ;
#define PI 3.14
#define BITS sizeof(int) * 8 // Total bits required to represent integer

extern  void pqrstDetection_22nov2018(const double ecg[6250], double fs, double
*pr_interval, double *qrs_interval, double *qt_interval, double *heart, double
                                      *qtc_interval, double *st_interval,bool change,int *count1,double *minVal,double *ramp,double *samp);

extern void b_filtfilt(const double x_in[6250], double y_out[6250]);
extern void filtfilt(const double x_in[6250], double y_out[6250]);
extern void flipud(double x[6286]);
extern void b_abs(const double x[6250], double y[6250]);




extern real_T rtInf;
extern real_T rtMinusInf;
extern real_T rtNaN;
extern real32_T rtInfF;
extern real32_T rtMinusInfF;
extern real32_T rtNaNF;
extern void rt_InitInfAndNaN(size_t realSize);
extern boolean_T rtIsInf(real_T value);
extern boolean_T rtIsInfF(real32_T value);
extern boolean_T rtIsNaN(real_T value);
extern boolean_T rtIsNaNF(real32_T value);



/* Function Declarations */
extern void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int
elementSize);
extern void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);
extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions);
extern void emxInit_boolean_T1(emxArray_boolean_T **pEmxArray, int numDimensions);
extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
extern void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);


extern void power(const double a[6250], double y[6250]);


/* Function Declarations */
extern double b_mean(const double x[8]);
extern double mean(const emxArray_real_T *x);


/* Function Declarations */
extern void findpeaks(const emxArray_real_T *Yin, double varargin_2,
                      emxArray_real_T *Ypk, emxArray_real_T *Xpk);
/* Function Declarations */
extern void conv(const double A[6250], const emxArray_real_T *B, emxArray_real_T
*C);

/* Function Declarations */
extern void b_diff(const double x[9], double y[8]);
extern void diff(const emxArray_real_T *x, emxArray_real_T *y);

/* Function Declarations */
extern boolean_T all(const emxArray_boolean_T *x);


/* Function Declarations */
extern int b_cfprintf(void);
extern int c_cfprintf(void);
extern int cfprintf(void);
extern int d_cfprintf(void);
extern int e_cfprintf(double varargin_1);
extern int f_cfprintf(void);
extern int g_cfprintf(double varargin_1);
extern int h_cfprintf(void);
extern int i_cfprintf(double varargin_1);
extern int j_cfprintf(void);
extern int k_cfprintf(double varargin_1);
extern int l_cfprintf(double varargin_1);
extern int m_cfprintf(void);
extern int n_cfprintf(double varargin_1);
extern int o_cfprintf(void);
extern int p_cfprintf(double varargin_1);
extern int q_cfprintf(void);
extern int r_cfprintf(double varargin_1);
extern int s_cfprintf(double varargin_1);
extern int t_cfprintf(double varargin_1);
extern int u_cfprintf(double varargin_1);
extern int v_cfprintf(double varargin_1);
extern int w_cfprintf(void);


/* Function Declarations */
extern void b_filter(const double b[7], const double x[6286], const double zi[6],
                     double y[6286]);
extern void filter(const double b[7], const double a[7], const double x[6286],
                   const double zi[6], double y[6286]);

/* Function Declarations */
extern void c_sort(emxArray_real_T *x);
extern void sort(emxArray_real_T *x, emxArray_int32_T *idx);

/* Function Declarations */
extern void b_fileManager(FILE * *f, boolean_T *a);
extern void c_fileManager(FILE * *f, boolean_T *a);
extern void d_fileManager(FILE * *f, boolean_T *a);
extern void e_fileManager(FILE * *f, boolean_T *a);
extern void f_fileManager(FILE * *f, boolean_T *a);
extern void fileManager(FILE * *f, boolean_T *a);
extern void g_fileManager(FILE * *f, boolean_T *a);
extern void h_fileManager(FILE * *f, boolean_T *a);
extern void i_fileManager(FILE * *f, boolean_T *a);
extern void j_fileManager(FILE * *f, boolean_T *a);
extern void k_fileManager(FILE * *f, boolean_T *a);
extern void l_fileManager(FILE * *f, boolean_T *a);
extern void m_fileManager(FILE * *f, boolean_T *a);
extern void n_fileManager(FILE * *f, boolean_T *a);
extern void o_fileManager(FILE * *f, boolean_T *a);
extern void p_fileManager(FILE * *f, boolean_T *a);
extern void q_fileManager(FILE * *f, boolean_T *a);
extern void r_fileManager(FILE * *f, boolean_T *a);
extern void s_fileManager(FILE * *f, boolean_T *a);
extern void t_fileManager(FILE * *f, boolean_T *a);
extern void u_fileManager(FILE * *f, boolean_T *a);
extern void v_fileManager(FILE * *f, boolean_T *a);
extern void w_fileManager(FILE * *f, boolean_T *a);

extern real_T rtGetNaN(void);
extern real32_T rtGetNaNF(void);

/* Function Declarations */
extern void do_vectors(const emxArray_real_T *a, const emxArray_real_T *b,
                       emxArray_real_T *c, emxArray_int32_T *ia, emxArray_int32_T *ib);

/* Function Declarations */
extern void b_sortIdx(emxArray_real_T *x, emxArray_int32_T *idx);
extern void sortIdx(emxArray_real_T *x, emxArray_int32_T *idx);


extern real_T rtGetInf(void);
extern real32_T rtGetInfF(void);
extern real_T rtGetMinusInf(void);
extern real32_T rtGetMinusInfF(void);


//const char *filename="\D:\Priyamvad_projects\ECG\ecg_signal_generation\ecg_6250_sample.txt";
//D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\ecg_6250_sample7.txt
//D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\ecg_individual_channels_samples\\ecg_pd1\\lead1.txt
//D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\ecg_individual_channels_samples\\ecg_pd2\\lead1.txt
//D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\ecg_individual_channels_samples\\hopsisoft_irenbhai_29oct2018\\lead1.txt
//D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\ecg_individual_channels_samples\\hopsisoft_nirzar_29oct2018\\lead1.txt
//D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\ecg_individual_channels_samples\\hopsisoft_nitin_29oct2018\\lead1.txt
//D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\simulator_18nov2018\\lead2.txt
int main (void)
{
    FILE *pf1,*pf2,*pf3,*pf4,*pf5,*pf6,*pf7,*pf8,*pf9,*pf10,*pf11,*pf12;
    double heart1=0.0,heart2=0.0,heart3=0.0,heart4=0.0,heart5=0.0,heart6=0.0,heart7=0.0,heart8=0.0;
    bool heartb1=false,heartb2=false,heartb3=false,heartb4=false,heartb5=false,heartb6=false,heartb7=false,heartb8=false,change=false;
    bool heartb9=false,heartb10=false,heartb11=false,heartb12=false,ld=false,avf=false,qrs=false;
    double ramp=0.0,ramp1=0.0,ramp2=0.0,ramp3=0.0,ramp4=0.0,ramp5=0.0,ramp6=0.0;
    double samp=0.0,samp1=0.0,samp2=0.0,samp3=0.0,samp4=0.0,samp5=0.0,samp6=0.0;
    double minValues[5],minVal1=0.0,minVal6=0.0;
    int checkHeart=0,count=0;
    double sampleFrequency=250;
    double minVal=0.0;
    double pr_interval=0.0;
    double qrs_interval=0.0;
    double qt_interval=0.0;
    double heart=0.0;
    double hr_avg=0.0;
    double qtc_interval=0.0;
    double st_interval=0.0;
    static double num1[6250];
    char path1[500],path2[500],path3[500],path4[500],path5[500],path6[500],path7[500],path8[500],path9[500],path10[500],path11[500],path12[500];
    printf("How to give path name:->  <main directory name>:\\\\<sub directory name>\\\\<file name><.txt>\n");
    printf("Example1:->  D:\\\\epsilon\\\\mysigmoid.txt\n");
    printf("Example2:->  D:\\\\epsilon\\\\hello\\\\project\\\\dft.txt\n\n");
    printf("****Average heartrate******\n");

    printf("Enter the proper path of file as given in above examples for lead1 channel:->");
    //fgets(path2,200,stdin);
    scanf("%499[^\n]%*c", path1);

    printf("Enter the proper path of file as given in above examples for lead2 channel:->");
    //fgets(path2,200,stdin);
    scanf("%499[^\n]%*c", path2);

    printf("Enter the proper path of file as given in above examples for channel3:->");
    //fgets(path3,200,stdin);
    scanf("%499[^\n]%*c", path3);

    printf("Enter the proper path of file as given in above examples for channel4:->");

    scanf("%499[^\n]%*c", path4);

    printf("Enter the proper path of file as given in above examples for channel5:->");

    scanf("%499[^\n]%*c", path5);

    printf("Enter the proper path of file as given in above examples for channel6:->");

    scanf("%499[^\n]%*c", path6);

    printf("Enter the proper path of file as given in above examples for channel7:->");

    scanf("%499[^\n]%*c", path7);

    printf("Enter the proper path of file as given in above examples for channel8:->");

    scanf("%499[^\n]%*c", path8);

    printf("Enter the proper path of file as given in above examples for channelaVF:->");

    scanf("%499[^\n]%*c", path12);


    memset(&minValues[0], 0, 6U * sizeof(double));
    pf1 = fopen(path1, "r");
    //D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\ecg_6250_sample.txt
    if (pf1 == NULL)
    {
        printf("Value of errno: %d\n", errno);
        perror("Lead1 file opening error:->");
        printf("\n");
        //printf("Error opening file: %s\n", strerror( errnum ));
    }
    else
    {
        printf("Lead 1 File opening successfull\n");
        //open file for reading
        //printf("Reading numbers from file\n");
        for(int i=0; i<2500; i++)
        {
            fscanf(pf1,"%lf", &num1[i]);
            //printf("value %d:->%lf\n",i,num1[i]);
        }
        if(num1[1])
        {
            pqrstDetection_22nov2018(num1,sampleFrequency,&pr_interval,&qrs_interval,&qt_interval,&heart,&qtc_interval,&st_interval,change,&checkHeart,&minVal,&ramp,&samp);
            change=true;
            heart1= heart;
            minVal1=minVal;
            ramp1=ramp;
            samp1=samp;
            //printf("Minvalues1:->%lf\n",minVal1);
            //printf("R amplitude:->%lf\n",ramp1);
            //printf("S amplitude:->%lf\n",samp1);
            printf("\n");
            //printf("heart1:%lf\n",heart1);
            heartb1=true;
            if(checkHeart==1)
            {
                ++count;
            }
            //printf("count:->%d\n",checkHeart);
            memset(&num1[0], 0, 2500U * sizeof(double));
        }
        else
        {
            printf("Please provide proper text file for lead 1 containing samples i.e samples values should be number!!!\n");
            heartb1=false;
        }
        //printf("\nEnd of file.\n");
    }
    fclose(pf1);

    pf2 = fopen(path2, "r");
    //D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\mysigmoid.m
    if (pf2 == NULL)
    {
        printf("Value of errno: %d\n", errno);
        perror("Lead 2 file opening error:->");
        printf("\n");
        //printf("Error opening file: %s\n", strerror( errnum ));
    }
    else
    {
        printf("Lead2 File opening successfull\n");
        //open file for reading
        //printf("Reading numbers from file\n");
        for(int i=0; i<6250; i++)
        {
            fscanf(pf2,"%lf", &num1[i]);
            //printf("value %d:->%lf\n",i,num1[i]);
        }
        if(num1[1])
        {

            pqrstDetection_22nov2018(num1,sampleFrequency,&pr_interval,&qrs_interval,&qt_interval,&heart,&qtc_interval,&st_interval,change,&checkHeart,&minVal,&ramp,&samp);
            change=false;
            heart2= heart;
            //printf("heart2:%lf\n",heart2);
            //minVal2=minVal;
            //ramp2=ramp;
            //samp2=samp;
            printf("\n");
            heartb2=true;
            if(checkHeart==1)
            {
                ++count;
            }
            //printf("count:->%d\n",checkHeart);
            memset(&num1[0], 0, 6250U * sizeof(double));
        }
        else
        {
            printf("Please provide proper text file for Lead2 containing samples i.e samples values should be number!!!\n");
            heartb2=false;
        }
        //printf("\nEnd of file.\n");
    }
    fclose(pf2);


    pf3 = fopen(path3, "r");
    //D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\mysigmoid.m
    if (pf3 == NULL)
    {
        printf("Value of errno: %d\n", errno);
        perror("Channel3 file opening error:->");
        printf("\n");
        //printf("Error opening file: %s\n", strerror( errnum ));
    }
    else
    {
        printf("Channel3 File opening successfull\n");
        //open file for reading
        //printf("Reading numbers from file\n");
        for(int i=0; i<6250; i++)
        {
            fscanf(pf3,"%lf", &num1[i]);
            //printf("value %d:->%lf\n",i,num1[i]);
        }
        if(num1[1])
        {

            pqrstDetection_22nov2018(num1,sampleFrequency,&pr_interval,&qrs_interval,&qt_interval,&heart,&qtc_interval,&st_interval,change,&checkHeart,&minVal,&ramp,&samp);
            change=false;
            heart3= heart;
            //printf("heart3:%lf\n",heart3);
            printf("\n");
            heartb3=true;
            if(checkHeart==1)
            {
                ++count;
            }
            //printf("count:->%d\n",checkHeart);
            memset(&num1[0], 0, 6250U * sizeof(double));
        }
        else
        {
            printf("Please provide proper text file for channel3 containing samples i.e samples values should be number!!!\n");
            heartb3=false;
        }
        //printf("\nEnd of file.\n");
    }
    fclose(pf3);

    pf4 = fopen(path4, "r");
    //D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\mysigmoid.m
    if (pf4 == NULL)
    {
        printf("Value of errno: %d\n", errno);
        perror("Channel4 file opening error:->");
        printf("\n");
        //printf("Error opening file: %s\n", strerror( errnum ));
    }
    else
    {
        printf("Channel4 File opening successfull\n");
        //open file for reading
        //printf("Reading numbers from file\n");
        for(int i=0; i<6250; i++)
        {
            fscanf(pf4,"%lf", &num1[i]);
            //printf("value %d:->%lf\n",i,num1[i]);
        }
        if(num1[1])
        {

            pqrstDetection_22nov2018(num1,sampleFrequency,&pr_interval,&qrs_interval,&qt_interval,&heart,&qtc_interval,&st_interval,change,&checkHeart,&minVal,&ramp,&samp);
            change=false;
            heart4= heart;
            //printf("heart4:%lf\n",heart4);
            printf("\n");
            heartb4=true;
            if(checkHeart==1)
            {
                ++count;
            }
            //printf("count:->%d\n",checkHeart);
            memset(&num1[0], 0, 6250U * sizeof(double));
        }
        else
        {
            printf("Please provide proper text file for channel4 containing samples i.e samples values should be number!!!\n");
            heartb4=false;
        }
        //printf("\nEnd of file.\n");
    }
    fclose(pf4);

    pf5 = fopen(path5, "r");
    //D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\mysigmoid.m
    if (pf5 == NULL)
    {
        printf("Value of errno: %d\n", errno);
        perror("Channel5 file opening error:->");
        printf("\n");
        //printf("Error opening file: %s\n", strerror( errnum ));
    }
    else
    {
        printf("Channel5 File opening successfull\n");
        //open file for reading
        //printf("Reading numbers from file\n");
        for(int i=0; i<6250; i++)
        {
            fscanf(pf5,"%lf", &num1[i]);
            //printf("value %d:->%lf\n",i,num1[i]);
        }
        if(num1[1])
        {
            pqrstDetection_22nov2018(num1,sampleFrequency,&pr_interval,&qrs_interval,&qt_interval,&heart,&qtc_interval,&st_interval,change,&checkHeart,&minVal,&ramp,&samp);
            change=false;
            heart5= heart;
            //printf("heart5:%lf\n",heart5);
            printf("\n");
            heartb5=true;
            if(checkHeart==1)
            {
                ++count;
            }
            //printf("count:->%d\n",checkHeart);
            memset(&num1[0], 0, 6250U * sizeof(double));
        }
        else
        {
            printf("Please provide proper text file for channel5 containing samples i.e samples values should be number!!!\n");
            heartb5=false;
        }
        //printf("\nEnd of file.\n");
    }
    fclose(pf5);

    pf6 = fopen(path6, "r");
    //D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\mysigmoid.m
    if (pf6 == NULL)
    {
        printf("Value of errno: %d\n", errno);
        perror("Channel6 file opening error:->");
        printf("\n");
        //printf("Error opening file: %s\n", strerror( errnum ));
    }
    else
    {
        printf("Channel6 File opening successfull\n");
        //open file for reading
        //printf("Reading numbers from file\n");
        for(int i=0; i<6250; i++)
        {
            fscanf(pf6,"%lf", &num1[i]);
            //printf("value %d:->%lf\n",i,num1[i]);
        }
        if(num1[1])
        {

            pqrstDetection_22nov2018(num1,sampleFrequency,&pr_interval,&qrs_interval,&qt_interval,&heart,&qtc_interval,&st_interval,change,&checkHeart,&minVal,&ramp,&samp);
            change=false;
            heart6= heart;
            printf("\n");
            //printf("heart6:%lf\n",heart6);
            heartb6=true;
            if(checkHeart==1)
            {
                ++count;
            }
            //printf("count:->%d\n",checkHeart);
            memset(&num1[0], 0, 6250U * sizeof(double));
        }
        else
        {
            printf("Please provide proper text file for channel6 containing samples i.e samples values should be number!!!\n");
            heartb6=false;
        }
        //printf("\nEnd of file.\n");
    }
    fclose(pf6);

    pf7 = fopen(path7, "r");
    //D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\mysigmoid.m
    if (pf7 == NULL)
    {
        printf("Value of errno: %d\n", errno);
        perror("Channel7 file opening error:->");
        printf("\n");
        //printf("Error opening file: %s\n", strerror( errnum ));
    }
    else
    {
        printf("Channel7 File opening successfull\n");
        //open file for reading
        //printf("Reading numbers from file\n");
        for(int i=0; i<6250; i++)
        {
            fscanf(pf7,"%lf", &num1[i]);
            //printf("value %d:->%lf\n",i,num1[i]);
        }
        if(num1[1])
        {

            pqrstDetection_22nov2018(num1,sampleFrequency,&pr_interval,&qrs_interval,&qt_interval,&heart,&qtc_interval,&st_interval,change,&checkHeart,&minVal,&ramp,&samp);
            change=true;
            heart7= heart;
            printf("\n");
            //printf("heart7:%lf\n",heart7);
            heartb7=true;
            if(checkHeart==1)
            {
                ++count;
            }
            //printf("count:->%d\n",checkHeart);
            memset(&num1[0], 0, 6250U * sizeof(double));
        }
        else
        {
            printf("Please provide proper text file for channel7 containing samples i.e samples values should be number!!!\n");
            heartb7=false;
        }
        //printf("\nEnd of file.\n");
    }
    fclose(pf7);

    pf8 = fopen(path8, "r");
    //D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\mysigmoid.m
    if (pf8 == NULL)
    {
        printf("Value of errno: %d\n", errno);
        perror("Channel8 file opening error:->");
        printf("\n");
        //printf("Error opening file: %s\n", strerror( errnum ));
    }
    else
    {

        printf("Channel8 File opening successfull\n");
        //open file for reading
        //printf("Reading numbers from file\n");
        for(int i=0; i<6250; i++)
        {
            fscanf(pf8,"%lf", &num1[i]);
            //printf("value %d:->%lf\n",i,num1[i]);
        }
        if(num1[1])
        {
            pqrstDetection_22nov2018(num1,sampleFrequency,&pr_interval,&qrs_interval,&qt_interval,&heart,&qtc_interval,&st_interval,change,&checkHeart,&minVal,&ramp,&samp);
            change=false;
            printf("\n");
            heart8= heart;
            //printf("heart8:%lf\n",heart8);
            heartb8=true;
            if(checkHeart==1)
            {
                ++count;
            }
            //printf("count:->%d\n",checkHeart);
            memset(&num1[0], 0, 6250U * sizeof(double));
        }
        else
        {
            printf("Please provide proper text file for channel8 containing samples i.e samples values should be number!!!\n");
            heartb8=false;
        }
        //printf("\nEnd of file.\n");
    }
    fclose(pf8);

    pf12 = fopen(path12, "r");
    //D:\\Priyamvad_projects\\ECG\\ecg_signal_generation\\mysigmoid.m
    if (pf12 == NULL)
    {
        printf("Value of errno: %d\n", errno);
        perror("aVF channel file opening error:->");
        printf("\n");
        //printf("Error opening file: %s\n", strerror( errnum ));
    }
    else
    {
        printf("aVF channel File opening successfull\n");
        //open file for reading
        //printf("Reading numbers from file\n");
        for(int i=0; i<6250; i++)
        {
            fscanf(pf12,"%lf", &num1[i]);
            //printf("value %d:->%lf\n",i,num1[i]);
        }
        if(num1[1])
        {
            pqrstDetection_22nov2018(num1,sampleFrequency,&pr_interval,&qrs_interval,&qt_interval,&heart,&qtc_interval,&st_interval,change,&checkHeart,&minVal,&ramp,&samp);
            change=false;
            //printf("minval:->%lf\n",minVal);
            minVal6=minVal;
            ramp6=ramp;
            samp6=samp;
            //printf("Minvalues6:->%lf\n",minVal6);
            //printf("R amplitude:->%lf\n",ramp6);
            //printf("S amplitude:->%lf\n",samp6);
            heartb12=true;
            printf("\n");
            //printf("heart8:%lf\n",heart8);
            memset(&num1[0], 0, 6250U * sizeof(double));
        }
        else
        {
            printf("Please provide proper text file for aVF channel containing samples i.e samples values should be number!!!\n");
            heartb12=false;
        }
        //printf("\nEnd of file.\n");
    }
    fclose(pf12);

    if((heartb1==true)&&(heartb2==true)&&(heartb3==true)&&(heartb4==true)&&(heartb5==true)&&(heartb6==true)&&(heartb7==true)&&(heartb8==true))
    {
        //printf("count:->%d\n",count);
        hr_avg=(heart1+heart2+heart3+heart4+heart5+heart6+heart7+heart8)/count;
        printf("Average heart rate:-> %lf\n",hr_avg);
    }
    else
    {
        printf("Average heartrate cannot be calculated as ");
        printf("there is problem in opening file from any of the eight channels\nPlease check whether file path is correct or data in file is numeric or not!!!\n");
    }

    /***Calculation of equivalent voltage(mV) generated by minimum R-S wave difference value from lead1 and aVF****/

    if((heartb1==true)&&(heartb12==true)&&(ramp1 > 0)&&(ramp6 > 0))

    {

        float result=0.0,div=0.0,mul;
        mul=180.0/PI;
        //Calculation of QRS axis
        div=minVal6/minVal1;
        result = atan(div);

        printf("QRS axis in radians:->%lf\n",result);

        // Converting radians to degrees
        result = result * mul;
        printf("QRS axis in degrees:->%lf\n",result);

    }
    else
    {
        printf("Error in handling either file containing lead1 samples or lead aVF samples.Please check path of file or data in file\n");
    }


    printf("End of program....Thank you...\n");
    printf("Press any key to come out of program\n");

    for(int i=0; i<20; i++)
    {
        printf(".");
    }
    getchar();

}




/* pqrst file Function Declarations */
//static double rt_powd_snf(double u0, double u1);
//static double rt_roundd_snf(double u);

/* Function Definitions */

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_powd_snf(double u0, double u1)
{
    double y;
    double d4;
    double d5;
    if (rtIsNaN(u0) || rtIsNaN(u1))
    {
        y = rtNaN;
    }
    else
    {
        d4 = fabs(u0);
        d5 = fabs(u1);
        if (rtIsInf(u1))
        {
            if (d4 == 1.0)
            {
                y = rtNaN;
            }
            else if (d4 > 1.0)
            {
                if (u1 > 0.0)
                {
                    y = rtInf;
                }
                else
                {
                    y = 0.0;
                }
            }
            else if (u1 > 0.0)
            {
                y = 0.0;
            }
            else
            {
                y = rtInf;
            }
        }
        else if (d5 == 0.0)
        {
            y = 1.0;
        }
        else if (d5 == 1.0)
        {
            if (u1 > 0.0)
            {
                y = u0;
            }
            else
            {
                y = 1.0 / u0;
            }
        }
        else if (u1 == 2.0)
        {
            y = u0 * u0;
        }
        else if ((u1 == 0.5) && (u0 >= 0.0))
        {
            y = sqrt(u0);
        }
        else if ((u0 < 0.0) && (u1 > floor(u1)))
        {
            y = rtNaN;
        }
        else
        {
            y = pow(u0, u1);
        }
    }

    return y;
}

/*
 * Arguments    : double u
 * Return Type  : double
 */
static double rt_roundd_snf(double u)
{
    double y;
    if (fabs(u) < 4.503599627370496E+15)
    {
        if (u >= 0.5)
        {
            y = floor(u + 0.5);
        }
        else if (u > -0.5)
        {
            y = u * 0.0;
        }
        else
        {
            y = ceil(u - 0.5);
        }
    }
    else
    {
        y = u;
    }

    return y;
}

/*
 * % function [R_i,R_amp,S_i,S_amp,T_i,T_amp]=peakdetect(ecg,fs,view)
 * % =================== Online Adaptive QRS detector ==================== %%
 * % ========================== Description ============================= %%
 *  QRS detection
 *  Detects Q , R and S waves,T Waves
 *  Uses the state-machine logic to determine different peaks in an ECG
 *  signal. It has the ability to confront noise by canceling out the noise
 *  by high pass filtering and baseline wander by low pass. Besides, check
 *  out criterion to stop detection of spikes.
 *  The code is written in a way for future online implementation.
 * % Inputs
 *  ecg : raw ecg vector
 *  fs : sampling frequency
 *  view : display results? (0: no, 1: Yes)
 * Arguments    : const double ecg[6250]
 *                double fs
 *                double *pr_interval
 *                double *qrs_interval
 *                double *qt_interval
 *                double *heart
 *                double *qtc_interval
 *                double *st_interval
 * Return Type  : void
 */
 void pqrstDetection_22nov2018(const double ecg[6250], double fs, double
*pr_interval, double *qrs_interval, double *qt_interval, double *heart, double
                              *qtc_interval, double *st_interval,bool change,int *count1,double *minVal,double *ramp,double *samp)
{
    double rr_interval;
    int skip;
    double m_selected_RR;
    double mean_RR;
    int ser_back;
    static double ecg_bpass[6250];
    static double ecg_d2[6250];
    int ixstart;
    double ndbl;
    int idx;
    boolean_T exitg18;
    int i0;
    boolean_T exitg17;
    emxArray_real_T *r0;
    double anew;
    int nm1d2;
    emxArray_real_T *ecg_m;
    emxArray_real_T *pks;
    emxArray_real_T *locs;
    static double dv0[6250];
    int n;
    emxArray_real_T *qrs_i;
    emxArray_real_T *qrs_i_raw;
    emxArray_real_T *qrs_amp_raw;
    double apnd;
    double THR_SIG;
    boolean_T exitg16;
    emxArray_real_T *b_ecg_m;
    double THR_NOISE;
    double SIG_LEV;
    double NOISE_LEV;
    double THR_SIG1;
    boolean_T exitg15;
    static double ecg_bpass_data[6250];
    int ecg_bpass_size[1];
    emxArray_real_T b_ecg_bpass_data;
    double THR_NOISE1;
    double SIG_LEV1;
    double NOISE_LEV1;
    unsigned int Beat_C;
    unsigned int Beat_C1;
    double y_i;
    int x_i;
    int i;
    emxArray_real_T *r1;
    emxArray_real_T *c_ecg_m;
    emxArray_real_T *d_ecg_m;
    int i1;
    int nx;
    double diffRR[8];
    boolean_T exitg13;
    boolean_T exitg14;
    emxArray_real_T *nois_c;
    double test_m;
    double comp;
    boolean_T exitg12;
    emxArray_real_T *nois_i;
    emxArray_real_T *SIGL_buf;
    boolean_T exitg11;
    double locs_temp;
    emxArray_real_T *NOISL_buf;
    double Slope1;
    emxArray_real_T *SIGL_buf1;
    emxArray_real_T *NOISL_buf1;
    double Slope2;
    double y_i_t;
    int x_i_t;
    boolean_T exitg10;
    boolean_T exitg9;
    emxArray_real_T *THRS_buf1;
    emxArray_real_T *THRS_buf;
    emxArray_real_T *Tamp;
    int j;
    emxArray_real_T *qrs_c;
    emxArray_boolean_T *x;
    emxArray_int32_T *ii;
    emxArray_boolean_T *b_qrs_c;
    emxArray_boolean_T *c_qrs_c;
    emxArray_boolean_T *d_qrs_c;
    emxArray_boolean_T *e_qrs_c;
    emxArray_int32_T *f_qrs_c;
    emxArray_int32_T *g_qrs_c;
    emxArray_int32_T *h_qrs_c;
    emxArray_int32_T *i_qrs_c;
    emxArray_int32_T *j_qrs_c;
    emxArray_int32_T *k_qrs_c;
    emxArray_int32_T *l_qrs_c;
    emxArray_int32_T *m_qrs_c;
    emxArray_int32_T *n_qrs_c;
    emxArray_int32_T *o_qrs_c;
    emxArray_int32_T *p_qrs_c;
    emxArray_int32_T *q_qrs_c;
    double cdiff;
    double absa;
    double absb;
    boolean_T exitg8;
    boolean_T exitg7;
    boolean_T guard4 = false;
    boolean_T exitg6;
    boolean_T exitg5;
    boolean_T guard3 = false;
    boolean_T exitg4;
    boolean_T exitg3;
    boolean_T guard2 = false;
    boolean_T exitg2;
    boolean_T exitg1;
    boolean_T guard1 = false;
    int rcount=0;
    int scount=0;

    /* % Outputs */
    /*  indexes and amplitudes of R_i, R_amp, etc */
    /*  heart_rate computed heart rate */
    /*  buffer_plot : processed signal */
    /* % ============== Licensce ========================================== %% */
    /*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS */
    /*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT */
    /*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS */
    /*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT */
    /*  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, */
    /*  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED */
    /*  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR */
    /*  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF */
    /*  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING */
    /*  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS */
    /*  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */
    /*  Author : */
    /*  Hooman Sedghamiz, Feb, 2018 */
    /*  MSc. Biomedical Engineering, Linkoping University */
    /*  Email : Hooman.sedghamiz@gmail.com */
    /* % Updates : */
    /*   Feb, 2018 : Clean up and fixes. */
    /* % ========================= initialize ============================ %% */
    /* leftold=zeros(1,length(ecg));%contains p wave values */
    /* 'pqrstDetection_22nov2018:44' pr_interval=0; */
    *pr_interval = 0.0;

    /* 'pqrstDetection_22nov2018:45' qrs_interval=0; */
    *qrs_interval = 0.0;

    /* 'pqrstDetection_22nov2018:46' qt_interval=0; */
    *qt_interval = 0.0;

    /* 'pqrstDetection_22nov2018:47' qtc_interval=0; */
    *qtc_interval = 0.0;

    /* 'pqrstDetection_22nov2018:48' st_interval=0; */
    *st_interval = 0.0;

    /* 'pqrstDetection_22nov2018:49' rr_interval=0; */
    rr_interval = 0.0;

    /* 'pqrstDetection_22nov2018:50' delay = 0; */
    /* 'pqrstDetection_22nov2018:51' skip = 0; */
    skip = 0;

    /*  becomes one when a T wave is detected */
    /* 'pqrstDetection_22nov2018:52' m_selected_RR = 0; */
    m_selected_RR = 0.0;

    /* 'pqrstDetection_22nov2018:53' mean_RR = 0; */
    mean_RR = 0.0;

    /* 'pqrstDetection_22nov2018:54' ser_back = 0; */
    ser_back = 0;

    /* % ==================== Noise cancelation(Filtering) =================== %% */
    /* [a,b] = butter(3,[0.004 0.36],'bandpass');     % bandpass filtering */
    /* ecg_h = filtfilt(a,b,ecg); */
    /* %  bandpass filter for Noise cancelation of other sampling frequencies(Filtering) */
    /* f1=1;                                                                      % cuttoff low frequency to get rid of baseline wander */
    /* f2=15;                                                                     % cuttoff frequency to discard high frequency noise */
    /* Wn=[f1 f2]*2/fs;%cuttoff frequencies                                       % cutt off based on fs */
    /* N = 3;                                                                     % order of 3 less processing */
    /* 'pqrstDetection_22nov2018:65' [a,b] = butter(3,[0.0080 0.1200],'bandpass'); */
    /* determine the cut off frequency                                                      % bandpass filtering */
    /* 'pqrstDetection_22nov2018:66' ecg_bpass = filtfilt(a,b,ecg); */
    filtfilt(ecg, ecg_bpass);

    /* filtering of ecg signal */
    /* 'pqrstDetection_22nov2018:67' ecg_bpass = ecg_bpass/ max( abs(ecg_bpass)); */
    b_abs(ecg_bpass, ecg_d2);
    ixstart = 1;
    ndbl = ecg_d2[0];
    if (rtIsNaN(ecg_d2[0]))
    {
        idx = 2;
        exitg18 = false;
        while ((!exitg18) && (idx < 6251))
        {
            ixstart = idx;
            if (!rtIsNaN(ecg_d2[idx - 1]))
            {
                ndbl = ecg_d2[idx - 1];
                exitg18 = true;
            }
            else
            {
                idx++;
            }
        }
    }

    if (ixstart < 6250)
    {
        while (ixstart + 1 < 6251)
        {
            if (ecg_d2[ixstart] > ndbl)
            {
                ndbl = ecg_d2[ixstart];
            }

            ixstart++;
        }
    }

    for (i0 = 0; i0 < 6250; i0++)
    {
        ecg_bpass[i0] /= ndbl;
    }

    /* max function returns largest element of array */
    /* sig=ecg/ max( abs(ecg)); */
    /* % ==================== derivative filter ========================== %% */
    /*  ------ H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2)) --------- % */
    /* ecg_d1 =filtfilt([31.2500   51.2500   45.0000    5.0000  -35.0000  -56.2500  -36.2500],1,ecg_h); */
    /* ecg_d1 = ecg_d1/max(ecg_d1); */
    /* 'pqrstDetection_22nov2018:76' ecg_d2=filtfilt([31.2500   51.2500   45.0000    5.0000  -35.0000  -56.2500  -36.2500],1,ecg_bpass); */
    b_filtfilt(ecg_bpass, ecg_d2);

    /* 'pqrstDetection_22nov2018:77' ecg_d2=ecg_d2/max(ecg_d2); */
    ixstart = 1;
    ndbl = ecg_d2[0];
    if (rtIsNaN(ecg_d2[0]))
    {
        idx = 2;
        exitg17 = false;
        while ((!exitg17) && (idx < 6251))
        {
            ixstart = idx;
            if (!rtIsNaN(ecg_d2[idx - 1]))
            {
                ndbl = ecg_d2[idx - 1];
                exitg17 = true;
            }
            else
            {
                idx++;
            }
        }
    }

    if (ixstart < 6250)
    {
        while (ixstart + 1 < 6251)
        {
            if (ecg_d2[ixstart] > ndbl)
            {
                ndbl = ecg_d2[ixstart];
            }

            ixstart++;
        }
    }

    for (i0 = 0; i0 < 6250; i0++)
    {
        ecg_d2[i0] /= ndbl;
    }

    emxInit_real_T1(&r0, 2);

    /* % ========== Squaring nonlinearly enhance the dominant peaks ========== %% */
    /* ecg_s1 = ecg_d1.^2; */
    /* 'pqrstDetection_22nov2018:82' ecg_s2= ecg_d2.^2; */
    /* % ============  Moving average ================== %% */
    /* -------Y(nt) = (1/N)[x(nT-(N - 1)T)+ x(nT - (N - 2)T)+...+x(nT)]---------% */
    /* 'pqrstDetection_22nov2018:87' ecg_m = conv(ecg_s2 ,ones(1 ,round(0.12*fs))/round(0.12*fs)); */
    ndbl = rt_roundd_snf(0.12 * fs);
    anew = rt_roundd_snf(0.12 * fs);
    i0 = r0->size[0] * r0->size[1];
    r0->size[0] = 1;
    r0->size[1] = (int)ndbl;
    emxEnsureCapacity((emxArray__common *)r0, i0, (int)sizeof(double));
    nm1d2 = (int)ndbl;
    for (i0 = 0; i0 < nm1d2; i0++)
    {
        r0->data[i0] = 1.0 / anew;
    }

    emxInit_real_T(&ecg_m, 1);
    emxInit_real_T1(&pks, 2);
    emxInit_real_T1(&locs, 2);
    power(ecg_d2, dv0);
    conv(dv0, r0, ecg_m);

    /* 'pqrstDetection_22nov2018:88' delay = delay + round(0.12*fs)/2; */
    /* % ===================== Fiducial Marks ============================== %% */
    /*  Note : a minimum distance of 40 samples is considered between each R wave */
    /*  since in physiological point of view no RR wave can occur in less than */
    /*  200 msec distance */
    /* 'pqrstDetection_22nov2018:95' [pks,locs] = findpeaks(ecg_m,'MINPEAKDISTANCE',round(0.16*fs)); */
    findpeaks(ecg_m, rt_roundd_snf(0.16 * fs), pks, locs);

    /* 1 */
    /* % =================== Initialize Some Other Parameters =============== %% */
    /* 'pqrstDetection_22nov2018:97' LLp = length(pks); */
    emxFree_real_T(&r0);
    nm1d2 = pks->size[0];
    if (nm1d2 >= 1)
    {
    }
    else
    {
        nm1d2 = 1;
    }

    if (pks->size[0] == 0)
    {
        n = 0;
    }
    else
    {
        n = nm1d2;
    }

    emxInit_real_T1(&qrs_i, 2);

    /*  ---------------- Stores QRS wrt Sig and Filtered Sig ------------------% */
    /* 'pqrstDetection_22nov2018:100' qrs_c = zeros(1,LLp); */
    /*  amplitude of R */
    /* 'pqrstDetection_22nov2018:101' qrs_i = zeros(1,LLp); */
    i0 = qrs_i->size[0] * qrs_i->size[1];
    qrs_i->size[0] = 1;
    qrs_i->size[1] = n;
    emxEnsureCapacity((emxArray__common *)qrs_i, i0, (int)sizeof(double));
    for (i0 = 0; i0 < n; i0++)
    {
        qrs_i->data[i0] = 0.0;
    }

    emxInit_real_T1(&qrs_i_raw, 2);

    /*  index */
    /* 'pqrstDetection_22nov2018:102' qrs_i_raw =zeros(1,LLp) ; */
    i0 = qrs_i_raw->size[0] * qrs_i_raw->size[1];
    qrs_i_raw->size[0] = 1;
    qrs_i_raw->size[1] = n;
    emxEnsureCapacity((emxArray__common *)qrs_i_raw, i0, (int)sizeof(double));
    for (i0 = 0; i0 < n; i0++)
    {
        qrs_i_raw->data[i0] = 0.0;
    }

    emxInit_real_T1(&qrs_amp_raw, 2);

    /*  amplitude of R */
    /* 'pqrstDetection_22nov2018:103' qrs_amp_raw= zeros(1,LLp); */
    i0 = qrs_amp_raw->size[0] * qrs_amp_raw->size[1];
    qrs_amp_raw->size[0] = 1;
    qrs_amp_raw->size[1] = n;
    emxEnsureCapacity((emxArray__common *)qrs_amp_raw, i0, (int)sizeof(double));
    for (i0 = 0; i0 < n; i0++)
    {
        qrs_amp_raw->data[i0] = 0.0;
    }

    /*  Index */
    /* coder.varsize('qrs_i_raw','qrs_amp_raw'); */
    /*  ------------------- Noise Buffers ---------------------------------% */
    /* 'pqrstDetection_22nov2018:107' nois_c = zeros(1,LLp); */
    /* 'pqrstDetection_22nov2018:108' nois_i = zeros(1,LLp); */
    /*  ------------------- Buffers for Signal and Noise ----------------- % */
    /* 'pqrstDetection_22nov2018:110' SIGL_buf = zeros(1,LLp); */
    /* 'pqrstDetection_22nov2018:111' NOISL_buf = zeros(1,LLp); */
    /* 'pqrstDetection_22nov2018:112' SIGL_buf1 = zeros(1,LLp); */
    /* 'pqrstDetection_22nov2018:113' NOISL_buf1 = zeros(1,LLp); */
    /* 'pqrstDetection_22nov2018:114' THRS_buf1 = zeros(1,LLp); */
    /* 'pqrstDetection_22nov2018:115' THRS_buf = zeros(1,LLp); */
    /* % initialize the training phase (2 seconds of the signal) to determine the THR_SIG and THR_NOISE */
    /* 'pqrstDetection_22nov2018:118' THR_SIG = max(ecg_m(1:2*fs))*1/3; */
    apnd = 2.0 * fs;
    if (1.0 > apnd)
    {
        nm1d2 = 0;
    }
    else
    {
        nm1d2 = (int)apnd;
    }

    ixstart = 1;
    ndbl = ecg_m->data[0];
    if (nm1d2 > 1)
    {
        if (rtIsNaN(ndbl))
        {
            idx = 2;
            exitg16 = false;
            while ((!exitg16) && (idx <= nm1d2))
            {
                ixstart = idx;
                if (!rtIsNaN(ecg_m->data[idx - 1]))
                {
                    ndbl = ecg_m->data[idx - 1];
                    exitg16 = true;
                }
                else
                {
                    idx++;
                }
            }
        }

        if (ixstart < nm1d2)
        {
            while (ixstart + 1 <= nm1d2)
            {
                if (ecg_m->data[ixstart] > ndbl)
                {
                    ndbl = ecg_m->data[ixstart];
                }

                ixstart++;
            }
        }
    }

    THR_SIG = ndbl / 3.0;

    /* 2                                          % 0.25 of the max amplitude */
    /* 'pqrstDetection_22nov2018:119' THR_NOISE = mean(ecg_m(1:2*fs))*1/2; */
    apnd = 2.0 * fs;
    if (1.0 > apnd)
    {
        nm1d2 = 0;
    }
    else
    {
        nm1d2 = (int)apnd;
    }

    emxInit_real_T(&b_ecg_m, 1);
    i0 = b_ecg_m->size[0];
    b_ecg_m->size[0] = nm1d2;
    emxEnsureCapacity((emxArray__common *)b_ecg_m, i0, (int)sizeof(double));
    for (i0 = 0; i0 < nm1d2; i0++)
    {
        b_ecg_m->data[i0] = ecg_m->data[i0];
    }

    ndbl = mean(b_ecg_m);
    THR_NOISE = ndbl / 2.0;

    /* 3                                       % 0.5 of the mean signal is considered to be noise */
    /* 'pqrstDetection_22nov2018:120' SIG_LEV= THR_SIG; */
    SIG_LEV = THR_SIG;

    /* 'pqrstDetection_22nov2018:121' NOISE_LEV = THR_NOISE; */
    NOISE_LEV = THR_NOISE;

    /* % Initialize bandpath filter threshold(2 seconds of the bandpass signal) */
    /* 'pqrstDetection_22nov2018:124' THR_SIG1 = max(ecg_bpass(1:2*fs))*1/3; */
    apnd = 2.0 * fs;
    emxFree_real_T(&b_ecg_m);
    if (1.0 > apnd)
    {
        nm1d2 = 0;
    }
    else
    {
        nm1d2 = (int)apnd;
    }

    ixstart = 1;
    ndbl = ecg_bpass[0];
    if (nm1d2 > 1)
    {
        if (rtIsNaN(ecg_bpass[0]))
        {
            idx = 2;
            exitg15 = false;
            while ((!exitg15) && (idx <= nm1d2))
            {
                ixstart = idx;
                if (!rtIsNaN(ecg_bpass[idx - 1]))
                {
                    ndbl = ecg_bpass[idx - 1];
                    exitg15 = true;
                }
                else
                {
                    idx++;
                }
            }
        }

        if (ixstart < nm1d2)
        {
            while (ixstart + 1 <= nm1d2)
            {
                if (ecg_bpass[ixstart] > ndbl)
                {
                    ndbl = ecg_bpass[ixstart];
                }

                ixstart++;
            }
        }
    }

    THR_SIG1 = ndbl / 3.0;

    /*  0.25 of the max amplitude */
    /* 'pqrstDetection_22nov2018:125' THR_NOISE1 = mean(ecg_bpass(1:2*fs))*1/2; */
    apnd = 2.0 * fs;
    if (1.0 > apnd)
    {
        nm1d2 = 0;
    }
    else
    {
        nm1d2 = (int)apnd;
    }

    ecg_bpass_size[0] = nm1d2;
    for (i0 = 0; i0 < nm1d2; i0++)
    {
        ecg_bpass_data[i0] = ecg_bpass[i0];
    }

    b_ecg_bpass_data.data = (double *)&ecg_bpass_data;
    b_ecg_bpass_data.size = (int *)&ecg_bpass_size;
    b_ecg_bpass_data.allocatedSize = 6250;
    b_ecg_bpass_data.numDimensions = 1;
    b_ecg_bpass_data.canFreeData = false;
    ndbl = mean(&b_ecg_bpass_data);
    THR_NOISE1 = ndbl / 2.0;

    /* 'pqrstDetection_22nov2018:126' SIG_LEV1 = THR_SIG1; */
    SIG_LEV1 = THR_SIG1;

    /*  Signal level in Bandpassed filter */
    /* 'pqrstDetection_22nov2018:127' NOISE_LEV1 = THR_NOISE1; */
    NOISE_LEV1 = THR_NOISE1;

    /*  Noise level in Bandpassed filter */
    /* % ============ Thresholding and desicion rule ============= %% */
    /* 'pqrstDetection_22nov2018:129' Beat_C = 0; */
    Beat_C = 0U;

    /*  Raw Beats */
    /* 'pqrstDetection_22nov2018:130' Beat_C1 = 0; */
    Beat_C1 = 0U;

    /*  Filtered Beats */
    /* 'pqrstDetection_22nov2018:131' Noise_Count = 0; */
    /*  Noise Counter */
    /* 'pqrstDetection_22nov2018:132' y_i=0; */
    y_i = 0.0;

    /* 'pqrstDetection_22nov2018:133' x_i=0; */
    x_i = 0;

    /* 'pqrstDetection_22nov2018:134' for i = 1 : LLp */
    i = 0;
    emxInit_real_T(&r1, 1);
    emxInit_real_T(&c_ecg_m, 1);
    emxInit_real_T(&d_ecg_m, 1);
    while (i <= n - 1)
    {
        /*     %% ===== locate the corresponding peak in the filtered signal === %% */
        /* 'pqrstDetection_22nov2018:136' if locs(i)-round(0.150*fs)>= 1 && locs(i)<= length(ecg_bpass) */
        if ((locs->data[i] - rt_roundd_snf(0.15 * fs) >= 1.0) && (locs->data[i] <=
                                                                  6250.0))
        {
            /* 'pqrstDetection_22nov2018:137' [y_i,x_i] = max(ecg_bpass(locs(i)-round(0.150*fs):locs(i))); */
            apnd = locs->data[i] - rt_roundd_snf(0.15 * fs);
            if (apnd > locs->data[i])
            {
                i0 = 0;
                i1 = 0;
            }
            else
            {
                i0 = (int)apnd - 1;
                i1 = (int)locs->data[i];
            }

            ixstart = 1;
            nm1d2 = i1 - i0;
            ndbl = ecg_bpass[i0];
            nx = 1;
            if (i1 - i0 > 1)
            {
                if (rtIsNaN(ecg_bpass[i0]))
                {
                    idx = 2;
                    exitg14 = false;
                    while ((!exitg14) && (idx <= nm1d2))
                    {
                        ixstart = idx;
                        if (!rtIsNaN(ecg_bpass[(i0 + idx) - 1]))
                        {
                            ndbl = ecg_bpass[(i0 + idx) - 1];
                            nx = idx;
                            exitg14 = true;
                        }
                        else
                        {
                            idx++;
                        }
                    }
                }

                if (ixstart < i1 - i0)
                {
                    while (ixstart + 1 <= nm1d2)
                    {
                        if (ecg_bpass[i0 + ixstart] > ndbl)
                        {
                            ndbl = ecg_bpass[i0 + ixstart];
                            nx = ixstart + 1;
                        }

                        ixstart++;
                    }
                }
            }

            y_i = ndbl;
            x_i = nx;
        }
        else
        {
            /* 'pqrstDetection_22nov2018:138' else */
            /* 'pqrstDetection_22nov2018:139' if i == 1 */
            if (1 + i == 1)
            {
                /* 'pqrstDetection_22nov2018:140' [y_i,x_i] = max(ecg_bpass(1:locs(i))); */
                if (1.0 > locs->data[0])
                {
                    nm1d2 = 0;
                }
                else
                {
                    nm1d2 = (int)locs->data[0];
                }

                ixstart = 1;
                ndbl = ecg_bpass[0];
                nx = 1;
                if (nm1d2 > 1)
                {
                    if (rtIsNaN(ecg_bpass[0]))
                    {
                        idx = 2;
                        exitg13 = false;
                        while ((!exitg13) && (idx <= nm1d2))
                        {
                            ixstart = idx;
                            if (!rtIsNaN(ecg_bpass[idx - 1]))
                            {
                                ndbl = ecg_bpass[idx - 1];
                                nx = idx;
                                exitg13 = true;
                            }
                            else
                            {
                                idx++;
                            }
                        }
                    }

                    if (ixstart < nm1d2)
                    {
                        while (ixstart + 1 <= nm1d2)
                        {
                            if (ecg_bpass[ixstart] > ndbl)
                            {
                                ndbl = ecg_bpass[ixstart];
                                nx = ixstart + 1;
                            }

                            ixstart++;
                        }
                    }
                }

                y_i = ndbl;
                x_i = nx;

                /* 'pqrstDetection_22nov2018:141' ser_back = 1; */
                ser_back = 1;
            }
            else
            {
                if (locs->data[i] >= 6250.0)
                {
                    /* 'pqrstDetection_22nov2018:142' elseif locs(i)>= length(ecg_bpass) */
                    /* 'pqrstDetection_22nov2018:143' [y_i,x_i] = max(ecg_bpass(locs(i)-round(0.150*fs):end)); */
                    apnd = locs->data[i] - rt_roundd_snf(0.15 * fs);
                    if (apnd > 6250.0)
                    {
                        i0 = 0;
                        i1 = 0;
                    }
                    else
                    {
                        i0 = (int)apnd - 1;
                        i1 = 6250;
                    }

                    ixstart = 1;
                    nm1d2 = i1 - i0;
                    ndbl = ecg_bpass[i0];
                    nx = 1;
                    if (i1 - i0 > 1)
                    {
                        if (rtIsNaN(ecg_bpass[i0]))
                        {
                            idx = 2;
                            exitg12 = false;
                            while ((!exitg12) && (idx <= nm1d2))
                            {
                                ixstart = idx;
                                if (!rtIsNaN(ecg_bpass[(i0 + idx) - 1]))
                                {
                                    ndbl = ecg_bpass[(i0 + idx) - 1];
                                    nx = idx;
                                    exitg12 = true;
                                }
                                else
                                {
                                    idx++;
                                }
                            }
                        }

                        if (ixstart < i1 - i0)
                        {
                            while (ixstart + 1 <= nm1d2)
                            {
                                if (ecg_bpass[i0 + ixstart] > ndbl)
                                {
                                    ndbl = ecg_bpass[i0 + ixstart];
                                    nx = ixstart + 1;
                                }

                                ixstart++;
                            }
                        }
                    }

                    y_i = ndbl;
                    x_i = nx;
                }
            }
        }

        /*     %% ================= update the heart_rate ==================== %% */
        /* 'pqrstDetection_22nov2018:147' if Beat_C >= 9 */
        if ((int)Beat_C >= 9)
        {
            /* 'pqrstDetection_22nov2018:148' diffRR = diff(qrs_i(Beat_C-8:Beat_C)); */
            b_diff(*(double (*)[9])&qrs_i->data[(int)Beat_C - 9], diffRR);

            /*  calculate RR interval */
            /* 'pqrstDetection_22nov2018:149' mean_RR = mean(diffRR); */
            mean_RR = b_mean(diffRR);

            /*  calculate the mean of 8 previous R waves interval */
            /* 'pqrstDetection_22nov2018:150' comp =qrs_i(Beat_C)-qrs_i(Beat_C-1); */
            comp = qrs_i->data[(int)Beat_C - 1] - qrs_i->data[(int)Beat_C - 2];

            /*  latest RR */
            /* 'pqrstDetection_22nov2018:152' if comp <= 0.92*mean_RR || comp >= 1.16*mean_RR */
            if ((comp <= 0.92 * mean_RR) || (comp >= 1.16 * mean_RR))
            {
                /*  ------ lower down thresholds to detect better in MVI -------- % */
                /* 'pqrstDetection_22nov2018:154' THR_SIG = 0.5*(THR_SIG); */
                THR_SIG *= 0.5;

                /* 'pqrstDetection_22nov2018:155' THR_SIG1 = 0.5*(THR_SIG1); */
                THR_SIG1 *= 0.5;
            }
            else
            {
                /* 'pqrstDetection_22nov2018:156' else */
                /* 'pqrstDetection_22nov2018:157' m_selected_RR = mean_RR; */
                m_selected_RR = mean_RR;

                /*  The latest regular beats mean */
            }
        }

        /*     %% == calculate the mean last 8 R waves to ensure that QRS is not ==== %% */
        /* 'pqrstDetection_22nov2018:163' if m_selected_RR */
        if (m_selected_RR != 0.0)
        {
            /* 'pqrstDetection_22nov2018:164' test_m = m_selected_RR; */
            test_m = m_selected_RR;

            /* if the regular RR availabe use it */
        }
        else if (mean_RR != 0.0)
        {
            /* 'pqrstDetection_22nov2018:165' elseif mean_RR && m_selected_RR == 0 */
            /* 'pqrstDetection_22nov2018:166' test_m = mean_RR; */
            test_m = mean_RR;
        }
        else
        {
            /* 'pqrstDetection_22nov2018:167' else */
            /* 'pqrstDetection_22nov2018:168' test_m = 0; */
            test_m = 0.0;
        }

        /* 'pqrstDetection_22nov2018:171' if test_m */
        if ((test_m != 0.0) && (locs->data[i] - qrs_i->data[(int)Beat_C - 1] >=
                                rt_roundd_snf(1.66 * test_m)))
        {
            /* 'pqrstDetection_22nov2018:172' if (locs(i) - qrs_i(Beat_C)) >= round(1.66*test_m) */
            /*  it shows a QRS is missed */
            /* 'pqrstDetection_22nov2018:173' [pks_temp,locs_temp] = max(ecg_m(qrs_i(Beat_C)+ round(0.200*fs):locs(i)-round(0.200*fs))); */
            apnd = qrs_i->data[(int)Beat_C - 1] + rt_roundd_snf(0.2 * fs);
            anew = locs->data[i] - rt_roundd_snf(0.2 * fs);
            if (apnd > anew)
            {
                i0 = 1;
                i1 = 1;
            }
            else
            {
                i0 = (int)apnd;
                i1 = (int)anew + 1;
            }

            ixstart = 1;
            nm1d2 = i1 - i0;
            ndbl = ecg_m->data[i0 - 1];
            nx = 1;
            if (i1 - i0 > 1)
            {
                if (rtIsNaN(ndbl))
                {
                    idx = 2;
                    exitg11 = false;
                    while ((!exitg11) && (idx <= nm1d2))
                    {
                        ixstart = idx;
                        if (!rtIsNaN(ecg_m->data[(i0 + idx) - 2]))
                        {
                            ndbl = ecg_m->data[(i0 + idx) - 2];
                            nx = idx;
                            exitg11 = true;
                        }
                        else
                        {
                            idx++;
                        }
                    }
                }

                if (ixstart < i1 - i0)
                {
                    for (idx = ixstart + 1; idx <= nm1d2; idx++)
                    {
                        if (ecg_m->data[(i0 + idx) - 2] > ndbl)
                        {
                            ndbl = ecg_m->data[(i0 + idx) - 2];
                            nx = idx;
                        }
                    }
                }
            }

            /*  search back and locate the max in this interval */
            /* 'pqrstDetection_22nov2018:174' locs_temp = qrs_i(Beat_C)+ round(0.200*fs) + locs_temp -1; */
            locs_temp = ((qrs_i->data[(int)Beat_C - 1] + rt_roundd_snf(0.2 * fs)) +
                         (double)nx) - 1.0;

            /*  location */
            /* 'pqrstDetection_22nov2018:176' if pks_temp > THR_NOISE */
            if (ndbl > THR_NOISE)
            {
                /* 'pqrstDetection_22nov2018:177' Beat_C = Beat_C + 1; */
                Beat_C++;

                /* 'pqrstDetection_22nov2018:178' qrs_c(Beat_C) = pks_temp; */
                /* 'pqrstDetection_22nov2018:179' qrs_i(Beat_C) = locs_temp; */
                qrs_i->data[(int)Beat_C - 1] = locs_temp;

                /*  ------------- Locate in Filtered Sig ------------- % */
                /* 'pqrstDetection_22nov2018:181' if locs_temp <= length(ecg_bpass) */
                if (locs_temp <= 6250.0)
                {
                    /* 'pqrstDetection_22nov2018:182' [y_i_t,x_i_t] = max(ecg_bpass(locs_temp-round(0.150*fs):locs_temp)); */
                    apnd = locs_temp - rt_roundd_snf(0.15 * fs);
                    if (apnd > locs_temp)
                    {
                        i0 = 0;
                        i1 = 0;
                    }
                    else
                    {
                        i0 = (int)apnd - 1;
                        i1 = (int)locs_temp;
                    }

                    ixstart = 1;
                    nm1d2 = i1 - i0;
                    anew = ecg_bpass[i0];
                    nx = 1;
                    if (i1 - i0 > 1)
                    {
                        if (rtIsNaN(ecg_bpass[i0]))
                        {
                            idx = 2;
                            exitg10 = false;
                            while ((!exitg10) && (idx <= nm1d2))
                            {
                                ixstart = idx;
                                if (!rtIsNaN(ecg_bpass[(i0 + idx) - 1]))
                                {
                                    anew = ecg_bpass[(i0 + idx) - 1];
                                    nx = idx;
                                    exitg10 = true;
                                }
                                else
                                {
                                    idx++;
                                }
                            }
                        }

                        if (ixstart < i1 - i0)
                        {
                            while (ixstart + 1 <= nm1d2)
                            {
                                if (ecg_bpass[i0 + ixstart] > anew)
                                {
                                    anew = ecg_bpass[i0 + ixstart];
                                    nx = ixstart + 1;
                                }

                                ixstart++;
                            }
                        }
                    }

                    y_i_t = anew;
                    x_i_t = nx;
                }
                else
                {
                    /* 'pqrstDetection_22nov2018:183' else */
                    /* 'pqrstDetection_22nov2018:184' [y_i_t,x_i_t] = max(ecg_bpass(locs_temp-round(0.150*fs):end)); */
                    apnd = locs_temp - rt_roundd_snf(0.15 * fs);
                    if (apnd > 6250.0)
                    {
                        i0 = 0;
                        i1 = 0;
                    }
                    else
                    {
                        i0 = (int)apnd - 1;
                        i1 = 6250;
                    }

                    ixstart = 1;
                    nm1d2 = i1 - i0;
                    anew = ecg_bpass[i0];
                    nx = 1;
                    if (i1 - i0 > 1)
                    {
                        if (rtIsNaN(ecg_bpass[i0]))
                        {
                            idx = 2;
                            exitg9 = false;
                            while ((!exitg9) && (idx <= nm1d2))
                            {
                                ixstart = idx;
                                if (!rtIsNaN(ecg_bpass[(i0 + idx) - 1]))
                                {
                                    anew = ecg_bpass[(i0 + idx) - 1];
                                    nx = idx;
                                    exitg9 = true;
                                }
                                else
                                {
                                    idx++;
                                }
                            }
                        }

                        if (ixstart < i1 - i0)
                        {
                            while (ixstart + 1 <= nm1d2)
                            {
                                if (ecg_bpass[i0 + ixstart] > anew)
                                {
                                    anew = ecg_bpass[i0 + ixstart];
                                    nx = ixstart + 1;
                                }

                                ixstart++;
                            }
                        }
                    }

                    y_i_t = anew;
                    x_i_t = nx;
                }

                /*  ----------- Band pass Sig Threshold ------------------% */
                /* 'pqrstDetection_22nov2018:187' if y_i_t > THR_NOISE1 */
                if (y_i_t > THR_NOISE1)
                {
                    /* 'pqrstDetection_22nov2018:188' Beat_C1 = Beat_C1 + 1; */
                    Beat_C1++;

                    /* 'pqrstDetection_22nov2018:189' qrs_i_raw(Beat_C1) = locs_temp-round(0.150*fs)+ (x_i_t - 1); */
                    qrs_i_raw->data[(int)Beat_C1 - 1] = (locs_temp - rt_roundd_snf(0.15 *
                                                                                   fs)) + ((double)x_i_t - 1.0);

                    /*  save index of bandpass */
                    /* 'pqrstDetection_22nov2018:190' qrs_amp_raw(Beat_C1) = y_i_t; */
                    qrs_amp_raw->data[(int)Beat_C1 - 1] = y_i_t;

                    /*  save amplitude of bandpass */
                    /* 'pqrstDetection_22nov2018:191' SIG_LEV1 = 0.25*y_i_t + 0.75*SIG_LEV1; */
                    SIG_LEV1 = 0.25 * y_i_t + 0.75 * SIG_LEV1;

                    /*  when found with the second thres */
                }

                /* 'pqrstDetection_22nov2018:194' not_nois = 1; */
                /* 'pqrstDetection_22nov2018:195' SIG_LEV = 0.25*pks_temp + 0.75*SIG_LEV ; */
                SIG_LEV = 0.25 * ndbl + 0.75 * SIG_LEV;

                /*  when found with the second threshold */
            }
        }
        else
        {
            /* 'pqrstDetection_22nov2018:197' else */
            /* 'pqrstDetection_22nov2018:198' not_nois = 0; */
        }

        /*     %% ===================  find noise and QRS peaks ================== %% */
        /* 'pqrstDetection_22nov2018:203' if pks(i) >= THR_SIG */
        if (pks->data[i] >= THR_SIG)
        {
            /*  ------ if No QRS in 360ms of the previous QRS See if T wave ------% */
            /* 'pqrstDetection_22nov2018:205' if Beat_C >= 3 */
            if (((int)Beat_C >= 3) && (locs->data[i] - qrs_i->data[(int)Beat_C - 1] <=
                                       rt_roundd_snf(0.36 * fs)))
            {
                /* 'pqrstDetection_22nov2018:206' if (locs(i)-qrs_i(Beat_C)) <= round(0.3600*fs) */
                /* 'pqrstDetection_22nov2018:207' Slope1 = mean(diff(ecg_m(locs(i)-round(0.075*fs):locs(i)))); */
                apnd = locs->data[i] - rt_roundd_snf(0.075 * fs);
                if (apnd > locs->data[i])
                {
                    i0 = 0;
                    i1 = 0;
                }
                else
                {
                    i0 = (int)apnd - 1;
                    i1 = (int)locs->data[i];
                }

                nm1d2 = d_ecg_m->size[0];
                d_ecg_m->size[0] = i1 - i0;
                emxEnsureCapacity((emxArray__common *)d_ecg_m, nm1d2, (int)sizeof(double));
                nm1d2 = i1 - i0;
                for (i1 = 0; i1 < nm1d2; i1++)
                {
                    d_ecg_m->data[i1] = ecg_m->data[i0 + i1];
                }

                diff(d_ecg_m, r1);
                Slope1 = mean(r1);

                /* 4       % mean slope of the waveform at that position */
                /* 'pqrstDetection_22nov2018:208' Slope2 = mean(diff(ecg_m(qrs_i(Beat_C)-round(0.075*fs):qrs_i(Beat_C)))); */
                apnd = qrs_i->data[(int)Beat_C - 1] - rt_roundd_snf(0.075 * fs);
                if (apnd > qrs_i->data[(int)Beat_C - 1])
                {
                    i0 = 0;
                    i1 = 0;
                }
                else
                {
                    i0 = (int)apnd - 1;
                    i1 = (int)qrs_i->data[(int)Beat_C - 1];
                }

                nm1d2 = c_ecg_m->size[0];
                c_ecg_m->size[0] = i1 - i0;
                emxEnsureCapacity((emxArray__common *)c_ecg_m, nm1d2, (int)sizeof(double));
                nm1d2 = i1 - i0;
                for (i1 = 0; i1 < nm1d2; i1++)
                {
                    c_ecg_m->data[i1] = ecg_m->data[i0 + i1];
                }

                diff(c_ecg_m, r1);
                Slope2 = mean(r1);

                /* 5 % mean slope of previous R wave */
                /* 'pqrstDetection_22nov2018:209' if abs(Slope1) <= abs(0.5*(Slope2)) */
                if (fabs(Slope1) <= fabs(0.5 * Slope2))
                {
                    /*  slope less then 0.5 of previous R */
                    /* 'pqrstDetection_22nov2018:210' Noise_Count = Noise_Count + 1; */
                    /* 'pqrstDetection_22nov2018:211' nois_c(Noise_Count) = pks(i); */
                    /* 'pqrstDetection_22nov2018:212' nois_i(Noise_Count) = locs(i); */
                    /* 'pqrstDetection_22nov2018:213' skip = 1; */
                    skip = 1;

                    /*  T wave identification */
                    /*  ----- adjust noise levels ------ % */
                    /* 'pqrstDetection_22nov2018:215' NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1; */
                    NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1;

                    /* 'pqrstDetection_22nov2018:216' NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV; */
                    NOISE_LEV = 0.125 * pks->data[i] + 0.875 * NOISE_LEV;
                }
                else
                {
                    /* 'pqrstDetection_22nov2018:217' else */
                    /* 'pqrstDetection_22nov2018:218' skip = 0; */
                }
            }

            /* ---------- skip is 1 when a T wave is detected -------------- % */
            /* 'pqrstDetection_22nov2018:224' if skip == 0 */
            if (skip == 0)
            {
                /* 'pqrstDetection_22nov2018:225' Beat_C = Beat_C + 1; */
                Beat_C++;

                /* 'pqrstDetection_22nov2018:226' qrs_c(Beat_C) = pks(i); */
                /* 'pqrstDetection_22nov2018:227' qrs_i(Beat_C) = locs(i); */
                qrs_i->data[(int)Beat_C - 1] = locs->data[i];

                /* --------------- bandpass filter check threshold --------------- % */
                /* 'pqrstDetection_22nov2018:230' if y_i >= THR_SIG1 */
                if (y_i >= THR_SIG1)
                {
                    /* 'pqrstDetection_22nov2018:231' Beat_C1 = Beat_C1 + 1; */
                    Beat_C1++;

                    /* 'pqrstDetection_22nov2018:232' if ser_back */
                    if (ser_back != 0)
                    {
                        /* 'pqrstDetection_22nov2018:233' qrs_i_raw(Beat_C1) = x_i; */
                        qrs_i_raw->data[(int)Beat_C1 - 1] = x_i;

                        /*  save index of bandpass */
                    }
                    else
                    {
                        /* 'pqrstDetection_22nov2018:234' else */
                        /* 'pqrstDetection_22nov2018:235' qrs_i_raw(Beat_C1)= locs(i)-round(0.150*fs)+ (x_i - 1); */
                        qrs_i_raw->data[(int)Beat_C1 - 1] = (locs->data[i] - rt_roundd_snf
                                (0.15 * fs)) + ((double)x_i - 1.0);

                        /*  save index of bandpass */
                    }

                    /* 'pqrstDetection_22nov2018:237' qrs_amp_raw(Beat_C1) =  y_i; */
                    qrs_amp_raw->data[(int)Beat_C1 - 1] = y_i;

                    /*  save amplitude of bandpass */
                    /* 'pqrstDetection_22nov2018:238' SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1; */
                    SIG_LEV1 = 0.125 * y_i + 0.875 * SIG_LEV1;

                    /*  adjust threshold for bandpass filtered sig */
                }

                /* 'pqrstDetection_22nov2018:240' SIG_LEV = 0.125*pks(i) + 0.875*SIG_LEV ; */
                SIG_LEV = 0.125 * pks->data[i] + 0.875 * SIG_LEV;

                /*  adjust Signal level */
            }
        }
        else if ((THR_NOISE <= pks->data[i]) && (pks->data[i] < THR_SIG))
        {
            /* 'pqrstDetection_22nov2018:243' elseif (THR_NOISE <= pks(i)) && (pks(i) < THR_SIG) */
            /* 'pqrstDetection_22nov2018:244' NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1; */
            NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1;

            /*  adjust Noise level in filtered sig */
            /* 'pqrstDetection_22nov2018:245' NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV; */
            NOISE_LEV = 0.125 * pks->data[i] + 0.875 * NOISE_LEV;

            /*  adjust Noise level in MVI */
        }
        else
        {
            if (pks->data[i] < THR_NOISE)
            {
                /* 'pqrstDetection_22nov2018:246' elseif pks(i) < THR_NOISE */
                /* 'pqrstDetection_22nov2018:247' Noise_Count = Noise_Count + 1; */
                /* 'pqrstDetection_22nov2018:248' nois_c(Noise_Count) = pks(i); */
                /* 'pqrstDetection_22nov2018:249' nois_i(Noise_Count) = locs(i); */
                /* 'pqrstDetection_22nov2018:250' NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1; */
                NOISE_LEV1 = 0.125 * y_i + 0.875 * NOISE_LEV1;

                /*  noise level in filtered signal */
                /* 'pqrstDetection_22nov2018:251' NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV; */
                NOISE_LEV = 0.125 * pks->data[i] + 0.875 * NOISE_LEV;

                /*  adjust Noise level in MVI */
            }
        }

        /*     %% ================== adjust the threshold with SNR ============= %% */
        /* 'pqrstDetection_22nov2018:255' if NOISE_LEV ~= 0 || SIG_LEV ~= 0 */
        if ((NOISE_LEV != 0.0) || (SIG_LEV != 0.0))
        {
            /* 'pqrstDetection_22nov2018:256' THR_SIG = NOISE_LEV + 0.25*(abs(SIG_LEV - NOISE_LEV)); */
            THR_SIG = NOISE_LEV + 0.25 * fabs(SIG_LEV - NOISE_LEV);

            /* 'pqrstDetection_22nov2018:257' THR_NOISE = 0.5*(THR_SIG); */
            THR_NOISE = 0.5 * THR_SIG;
        }

        /* ------ adjust the threshold with SNR for bandpassed signal -------- % */
        /* 'pqrstDetection_22nov2018:261' if NOISE_LEV1 ~= 0 || SIG_LEV1 ~= 0 */
        if ((NOISE_LEV1 != 0.0) || (SIG_LEV1 != 0.0))
        {
            /* 'pqrstDetection_22nov2018:262' THR_SIG1 = NOISE_LEV1 + 0.25*(abs(SIG_LEV1 - NOISE_LEV1)); */
            THR_SIG1 = NOISE_LEV1 + 0.25 * fabs(SIG_LEV1 - NOISE_LEV1);

            /* 'pqrstDetection_22nov2018:263' THR_NOISE1 = 0.5*(THR_SIG1); */
            THR_NOISE1 = 0.5 * THR_SIG1;
        }

        /* --------- take a track of thresholds of smoothed signal -------------% */
        /* 'pqrstDetection_22nov2018:268' SIGL_buf(i) = SIG_LEV; */
        /* 'pqrstDetection_22nov2018:269' NOISL_buf(i) = NOISE_LEV; */
        /* 'pqrstDetection_22nov2018:270' THRS_buf(i) = THR_SIG; */
        /* -------- take a track of thresholds of filtered signal ----------- % */
        /* 'pqrstDetection_22nov2018:272' SIGL_buf1(i) = SIG_LEV1; */
        /* 'pqrstDetection_22nov2018:273' NOISL_buf1(i) = NOISE_LEV1; */
        /* 'pqrstDetection_22nov2018:274' THRS_buf1(i) = THR_SIG1; */
        /*  ----------------------- reset parameters -------------------------- % */
        /* 'pqrstDetection_22nov2018:276' skip = 0; */
        skip = 0;

        /* 'pqrstDetection_22nov2018:277' not_nois = 0; */
        /* 'pqrstDetection_22nov2018:278' ser_back = 0; */
        ser_back = 0;
        i++;
    }

    emxFree_real_T(&d_ecg_m);
    emxFree_real_T(&c_ecg_m);
    emxFree_real_T(&r1);
    emxFree_real_T(&locs);
    emxFree_real_T(&pks);

    /* disp(Beat_C); */
    /* % ==============  define two buffers ================= %% */
    /* buffer_mean=mean(abs(ecg_h(1:2*fs)-mean(ecg_h(1:2*fs)))); %1.9803e+03     % adaptive threshold DC corrected (baseline removed) */
    /* buffer_T = mean(ecg_h(1:2*fs));  % -22.8507                      % second adaptive threshold to be used for T wave detection */
    /* 'pqrstDetection_22nov2018:286' Rloc=zeros(1,Beat_C-2); */
    i0 = qrs_i->size[0] * qrs_i->size[1];
    qrs_i->size[0] = 1;
    qrs_i->size[1] = (int)Beat_C - 2;
    emxEnsureCapacity((emxArray__common *)qrs_i, i0, (int)sizeof(double));
    nm1d2 = (int)Beat_C - 2;
    for (i0 = 0; i0 < nm1d2; i0++)
    {
        qrs_i->data[i0] = 0.0;
    }

    emxInit_real_T1(&nois_c, 2);

    /* 'pqrstDetection_22nov2018:287' Ramp=zeros(1,Beat_C-2); */
    i0 = nois_c->size[0] * nois_c->size[1];
    nois_c->size[0] = 1;
    nois_c->size[1] = (int)Beat_C - 2;
    emxEnsureCapacity((emxArray__common *)nois_c, i0, (int)sizeof(double));
    nm1d2 = (int)Beat_C - 2;
    for (i0 = 0; i0 < nm1d2; i0++)
    {
        nois_c->data[i0] = 0.0;
    }

    emxInit_real_T1(&nois_i, 2);

    /* 'pqrstDetection_22nov2018:288' Ploc=zeros(1,Beat_C-2); */
    i0 = nois_i->size[0] * nois_i->size[1];
    nois_i->size[0] = 1;
    nois_i->size[1] = (int)Beat_C - 2;
    emxEnsureCapacity((emxArray__common *)nois_i, i0, (int)sizeof(double));
    nm1d2 = (int)Beat_C - 2;
    for (i0 = 0; i0 < nm1d2; i0++)
    {
        nois_i->data[i0] = 0.0;
    }

    emxInit_real_T1(&SIGL_buf, 2);

    /* 'pqrstDetection_22nov2018:289' Pamp=zeros(1,Beat_C-2); */
    i0 = SIGL_buf->size[0] * SIGL_buf->size[1];
    SIGL_buf->size[0] = 1;
    SIGL_buf->size[1] = (int)Beat_C - 2;
    emxEnsureCapacity((emxArray__common *)SIGL_buf, i0, (int)sizeof(double));
    nm1d2 = (int)Beat_C - 2;
    for (i0 = 0; i0 < nm1d2; i0++)
    {
        SIGL_buf->data[i0] = 0.0;
    }

    emxInit_real_T1(&NOISL_buf, 2);

    /* 'pqrstDetection_22nov2018:290' Qloc=zeros(1,Beat_C-2); */
    i0 = NOISL_buf->size[0] * NOISL_buf->size[1];
    NOISL_buf->size[0] = 1;
    NOISL_buf->size[1] = (int)Beat_C - 2;
    emxEnsureCapacity((emxArray__common *)NOISL_buf, i0, (int)sizeof(double));
    nm1d2 = (int)Beat_C - 2;
    for (i0 = 0; i0 < nm1d2; i0++)
    {
        NOISL_buf->data[i0] = 0.0;
    }

    emxInit_real_T1(&SIGL_buf1, 2);

    /* 'pqrstDetection_22nov2018:291' Qamp=zeros(1,Beat_C-2); */
    i0 = SIGL_buf1->size[0] * SIGL_buf1->size[1];
    SIGL_buf1->size[0] = 1;
    SIGL_buf1->size[1] = (int)Beat_C - 2;
    emxEnsureCapacity((emxArray__common *)SIGL_buf1, i0, (int)sizeof(double));
    nm1d2 = (int)Beat_C - 2;
    for (i0 = 0; i0 < nm1d2; i0++)
    {
        SIGL_buf1->data[i0] = 0.0;
    }

    emxInit_real_T1(&NOISL_buf1, 2);

    /* 'pqrstDetection_22nov2018:292' Sloc=zeros(1,Beat_C-2); */
    i0 = NOISL_buf1->size[0] * NOISL_buf1->size[1];
    NOISL_buf1->size[0] = 1;
    NOISL_buf1->size[1] = (int)Beat_C - 2;
    emxEnsureCapacity((emxArray__common *)NOISL_buf1, i0, (int)sizeof(double));
    nm1d2 = (int)Beat_C - 2;
    for (i0 = 0; i0 < nm1d2; i0++)
    {
        NOISL_buf1->data[i0] = 0.0;
    }

    emxInit_real_T1(&THRS_buf1, 2);

    /* 'pqrstDetection_22nov2018:293' Samp=zeros(1,Beat_C-2); */
    i0 = THRS_buf1->size[0] * THRS_buf1->size[1];
    THRS_buf1->size[0] = 1;
    THRS_buf1->size[1] = (int)Beat_C - 2;
    emxEnsureCapacity((emxArray__common *)THRS_buf1, i0, (int)sizeof(double));
    nm1d2 = (int)Beat_C - 2;
    for (i0 = 0; i0 < nm1d2; i0++)
    {
        THRS_buf1->data[i0] = 0.0;
    }

    emxInit_real_T1(&THRS_buf, 2);

    /* 'pqrstDetection_22nov2018:294' Tloc=zeros(1,Beat_C-2); */
    i0 = THRS_buf->size[0] * THRS_buf->size[1];
    THRS_buf->size[0] = 1;
    THRS_buf->size[1] = (int)Beat_C - 2;
    emxEnsureCapacity((emxArray__common *)THRS_buf, i0, (int)sizeof(double));
    nm1d2 = (int)Beat_C - 2;
    for (i0 = 0; i0 < nm1d2; i0++)
    {
        THRS_buf->data[i0] = 0.0;
    }

    emxInit_real_T1(&Tamp, 2);

    /* 'pqrstDetection_22nov2018:295' Tamp=zeros(1,Beat_C-2); */
    i0 = Tamp->size[0] * Tamp->size[1];
    Tamp->size[0] = 1;
    Tamp->size[1] = (int)Beat_C - 2;
    emxEnsureCapacity((emxArray__common *)Tamp, i0, (int)sizeof(double));
    nm1d2 = (int)Beat_C - 2;
    for (i0 = 0; i0 < nm1d2; i0++)
    {
        Tamp->data[i0] = 0.0;
    }

    /* % ================== Counters ============================ %% */
    /* 'pqrstDetection_22nov2018:298' for i=1:Beat_C-2 */
    for (i = 0; i <= (int)Beat_C - 3; i++)
    {
        /* 'pqrstDetection_22nov2018:299' Rloc(i)=qrs_i_raw(i+1); */
        qrs_i->data[i] = qrs_i_raw->data[1 + i];

        /* 'pqrstDetection_22nov2018:300' Ramp(i)=qrs_amp_raw(i+1); */
        nois_c->data[i] = qrs_amp_raw->data[1 + i];
    }

    emxFree_real_T(&qrs_amp_raw);
    emxFree_real_T(&qrs_i_raw);

    /* disp((Rloc)); */
    /* disp((Ramp)); */
    /* 'pqrstDetection_22nov2018:305' X=Rloc; */
    /* 'pqrstDetection_22nov2018:306' y1=ecg_bpass; */
    /* disp(length(y1)); */
    /* 'pqrstDetection_22nov2018:308' for i=1:1:1 */
    /*  If you have a 12 lead data than, for(i=1:1:12) */
    /* 'pqrstDetection_22nov2018:310' for j=1:1:length(X) */
    j = 0;
    emxInit_real_T1(&qrs_c, 2);
    emxInit_boolean_T(&x, 1);
    emxInit_int32_T(&ii, 1);
    emxInit_boolean_T1(&b_qrs_c, 2);
    emxInit_boolean_T1(&c_qrs_c, 2);
    emxInit_boolean_T1(&d_qrs_c, 2);
    emxInit_boolean_T1(&e_qrs_c, 2);
    emxInit_int32_T(&f_qrs_c, 1);
    emxInit_int32_T(&g_qrs_c, 1);
    emxInit_int32_T(&h_qrs_c, 1);
    emxInit_int32_T(&i_qrs_c, 1);
    emxInit_int32_T(&j_qrs_c, 1);
    emxInit_int32_T(&k_qrs_c, 1);
    emxInit_int32_T(&l_qrs_c, 1);
    emxInit_int32_T(&m_qrs_c, 1);
    emxInit_int32_T(&n_qrs_c, 1);
    emxInit_int32_T(&o_qrs_c, 1);
    emxInit_int32_T(&p_qrs_c, 1);
    emxInit_int32_T(&q_qrs_c, 1);
    while (j <= qrs_i->size[1] - 1)
    {
        /* P detection */
        /* 'pqrstDetection_22nov2018:313' a=(Rloc(i,j)-55:Rloc(i,j)-20); */
        if (rtIsNaN(qrs_i->data[qrs_i->size[0] * j] - 55.0) || rtIsNaN(qrs_i->
                data[qrs_i->size[0] * j] - 20.0))
        {
            n = 1;
            anew = rtNaN;
            apnd = qrs_i->data[qrs_i->size[0] * j] - 20.0;
        }
        else if (qrs_i->data[qrs_i->size[0] * j] - 20.0 < qrs_i->data[qrs_i->size
                                                                      [0] * j] - 55.0)
        {
            n = 0;
            anew = qrs_i->data[qrs_i->size[0] * j] - 55.0;
            apnd = qrs_i->data[qrs_i->size[0] * j] - 20.0;
        }
        else if (rtIsInf(qrs_i->data[qrs_i->size[0] * j] - 55.0) || rtIsInf
                (qrs_i->data[qrs_i->size[0] * j] - 20.0))
        {
            n = 1;
            anew = rtNaN;
            apnd = qrs_i->data[qrs_i->size[0] * j] - 20.0;
        }
        else
        {
            anew = qrs_i->data[qrs_i->size[0] * j] - 55.0;
            ndbl = floor(((qrs_i->data[qrs_i->size[0] * j] - 20.0) - (qrs_i->
                    data[qrs_i->size[0] * j] - 55.0)) + 0.5);
            apnd = (qrs_i->data[qrs_i->size[0] * j] - 55.0) + ndbl;
            cdiff = apnd - (qrs_i->data[qrs_i->size[0] * j] - 20.0);
            absa = fabs(qrs_i->data[qrs_i->size[0] * j] - 55.0);
            absb = fabs(qrs_i->data[qrs_i->size[0] * j] - 20.0);
            if ((absa >= absb) || rtIsNaN(absb))
            {
                absb = absa;
            }

            if (fabs(cdiff) < 4.4408920985006262E-16 * absb)
            {
                ndbl++;
                apnd = qrs_i->data[qrs_i->size[0] * j] - 20.0;
            }
            else if (cdiff > 0.0)
            {
                apnd = (qrs_i->data[qrs_i->size[0] * j] - 55.0) + (ndbl - 1.0);
            }
            else
            {
                ndbl++;
            }

            if (ndbl >= 0.0)
            {
                n = (int)ndbl;
            }
            else
            {
                n = 0;
            }
        }

        i0 = qrs_c->size[0] * qrs_c->size[1];
        qrs_c->size[0] = 1;
        qrs_c->size[1] = n;
        emxEnsureCapacity((emxArray__common *)qrs_c, i0, (int)sizeof(double));
        if (n > 0)
        {
            qrs_c->data[0] = anew;
            if (n > 1)
            {
                qrs_c->data[n - 1] = apnd;
                nm1d2 = (n - 1) / 2;
                for (nx = 1; nx < nm1d2; nx++)
                {
                    qrs_c->data[nx] = anew + (double)nx;
                    qrs_c->data[(n - nx) - 1] = apnd - (double)nx;
                }

                if (nm1d2 << 1 == n - 1)
                {
                    qrs_c->data[nm1d2] = (anew + apnd) / 2.0;
                }
                else
                {
                    qrs_c->data[nm1d2] = anew + (double)nm1d2;
                    qrs_c->data[nm1d2 + 1] = apnd - (double)nm1d2;
                }
            }
        }

        /* 'pqrstDetection_22nov2018:315' if all(a < 6250) */
        i0 = e_qrs_c->size[0] * e_qrs_c->size[1];
        e_qrs_c->size[0] = 1;
        e_qrs_c->size[1] = qrs_c->size[1];
        emxEnsureCapacity((emxArray__common *)e_qrs_c, i0, (int)sizeof(boolean_T));
        nm1d2 = qrs_c->size[0] * qrs_c->size[1];
        for (i0 = 0; i0 < nm1d2; i0++)
        {
            e_qrs_c->data[i0] = (qrs_c->data[i0] < 6250.0);
        }

        if (all(e_qrs_c))
        {
            /* 'pqrstDetection_22nov2018:316' m=max(y1(a)); */
            ixstart = 1;
            i0 = f_qrs_c->size[0];
            f_qrs_c->size[0] = qrs_c->size[1];
            emxEnsureCapacity((emxArray__common *)f_qrs_c, i0, (int)sizeof(int));
            nm1d2 = qrs_c->size[1];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                f_qrs_c->data[i0] = (int)qrs_c->data[qrs_c->size[0] * i0];
            }

            n = f_qrs_c->size[0];
            ndbl = ecg_bpass[(int)qrs_c->data[0] - 1];
            i0 = g_qrs_c->size[0];
            g_qrs_c->size[0] = qrs_c->size[1];
            emxEnsureCapacity((emxArray__common *)g_qrs_c, i0, (int)sizeof(int));
            nm1d2 = qrs_c->size[1];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                g_qrs_c->data[i0] = (int)qrs_c->data[qrs_c->size[0] * i0];
            }

            if (g_qrs_c->size[0] > 1)
            {
                if (rtIsNaN(ndbl))
                {
                    idx = 2;
                    exitg8 = false;
                    while ((!exitg8) && (idx <= n))
                    {
                        ixstart = idx;
                        if (!rtIsNaN(ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * (idx - 1)]
                                               - 1]))
                        {
                            ndbl = ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * (idx - 1)] - 1];
                            exitg8 = true;
                        }
                        else
                        {
                            idx++;
                        }
                    }
                }

                i0 = h_qrs_c->size[0];
                h_qrs_c->size[0] = qrs_c->size[1];
                emxEnsureCapacity((emxArray__common *)h_qrs_c, i0, (int)sizeof(int));
                nm1d2 = qrs_c->size[1];
                for (i0 = 0; i0 < nm1d2; i0++)
                {
                    h_qrs_c->data[i0] = (int)qrs_c->data[qrs_c->size[0] * i0];
                }

                if (ixstart < h_qrs_c->size[0])
                {
                    while (ixstart + 1 <= n)
                    {
                        if (ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * ixstart] - 1] > ndbl)
                        {
                            ndbl = ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * ixstart] - 1];
                        }

                        ixstart++;
                    }
                }
            }

            /* 'pqrstDetection_22nov2018:318' b=find(y1(a)==m); */
            i0 = x->size[0];
            x->size[0] = qrs_c->size[1];
            emxEnsureCapacity((emxArray__common *)x, i0, (int)sizeof(boolean_T));
            nm1d2 = qrs_c->size[1];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                x->data[i0] = (ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * i0] - 1] ==
                               ndbl);
            }

            nx = x->size[0];
            idx = 0;
            i0 = ii->size[0];
            ii->size[0] = x->size[0];
            emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
            nm1d2 = 1;
            exitg7 = false;
            while ((!exitg7) && (nm1d2 <= nx))
            {
                guard4 = false;
                if (x->data[nm1d2 - 1])
                {
                    idx++;
                    ii->data[idx - 1] = nm1d2;
                    if (idx >= nx)
                    {
                        exitg7 = true;
                    }
                    else
                    {
                        guard4 = true;
                    }
                }
                else
                {
                    guard4 = true;
                }

                if (guard4)
                {
                    nm1d2++;
                }
            }

            if (x->size[0] == 1)
            {
                if (idx == 0)
                {
                    i0 = ii->size[0];
                    ii->size[0] = 0;
                    emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
                }
            }
            else
            {
                i0 = ii->size[0];
                if (1 > idx)
                {
                    ii->size[0] = 0;
                }
                else
                {
                    ii->size[0] = idx;
                }

                emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
            }

            i0 = ecg_m->size[0];
            ecg_m->size[0] = ii->size[0];
            emxEnsureCapacity((emxArray__common *)ecg_m, i0, (int)sizeof(double));
            nm1d2 = ii->size[0];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                ecg_m->data[i0] = ii->data[i0];
            }

            /* 'pqrstDetection_22nov2018:320' b=b(1); */
            /* 'pqrstDetection_22nov2018:322' b=a(b); */
            /* 'pqrstDetection_22nov2018:324' Ploc(i,j)=b-12; */
            nois_i->data[nois_i->size[0] * j] = qrs_c->data[(int)ecg_m->data[0] - 1] -
                                                12.0;

            /* 'pqrstDetection_22nov2018:326' Pamp(i,j)=m; */
            SIGL_buf->data[SIGL_buf->size[0] * j] = ndbl;

            /* disp(Pamp); */
        }
        else
        {
            /* 'pqrstDetection_22nov2018:328' else */
            /* 'pqrstDetection_22nov2018:329' fprintf('Not able to detect p wave\n'); */
            cfprintf();

            /* check if first p wave value is less than 50 or last p wave value is zero */
        }

        /* The minima in the Window of   Rloc-100 to Rloc-10  is essentially the Q peak. */
        /*         %% Q  Detection */
        /* 'pqrstDetection_22nov2018:336' a=(Rloc(i,j)-30:Rloc(i,j)-5); */
        if (rtIsNaN(qrs_i->data[qrs_i->size[0] * j] - 30.0) || rtIsNaN(qrs_i->
                data[qrs_i->size[0] * j] - 5.0))
        {
            n = 1;
            anew = rtNaN;
            apnd = qrs_i->data[qrs_i->size[0] * j] - 5.0;
        }
        else if (qrs_i->data[qrs_i->size[0] * j] - 5.0 < qrs_i->data[qrs_i->size[0]
                                                                     * j] - 30.0)
        {
            n = 0;
            anew = qrs_i->data[qrs_i->size[0] * j] - 30.0;
            apnd = qrs_i->data[qrs_i->size[0] * j] - 5.0;
        }
        else if (rtIsInf(qrs_i->data[qrs_i->size[0] * j] - 30.0) || rtIsInf
                (qrs_i->data[qrs_i->size[0] * j] - 5.0))
        {
            n = 1;
            anew = rtNaN;
            apnd = qrs_i->data[qrs_i->size[0] * j] - 5.0;
        }
        else
        {
            anew = qrs_i->data[qrs_i->size[0] * j] - 30.0;
            ndbl = floor(((qrs_i->data[qrs_i->size[0] * j] - 5.0) - (qrs_i->data
                                                                     [qrs_i->size[0] * j] - 30.0)) + 0.5);
            apnd = (qrs_i->data[qrs_i->size[0] * j] - 30.0) + ndbl;
            cdiff = apnd - (qrs_i->data[qrs_i->size[0] * j] - 5.0);
            absa = fabs(qrs_i->data[qrs_i->size[0] * j] - 30.0);
            absb = fabs(qrs_i->data[qrs_i->size[0] * j] - 5.0);
            if ((absa >= absb) || rtIsNaN(absb))
            {
                absb = absa;
            }

            if (fabs(cdiff) < 4.4408920985006262E-16 * absb)
            {
                ndbl++;
                apnd = qrs_i->data[qrs_i->size[0] * j] - 5.0;
            }
            else if (cdiff > 0.0)
            {
                apnd = (qrs_i->data[qrs_i->size[0] * j] - 30.0) + (ndbl - 1.0);
            }
            else
            {
                ndbl++;
            }

            if (ndbl >= 0.0)
            {
                n = (int)ndbl;
            }
            else
            {
                n = 0;
            }
        }

        i0 = qrs_c->size[0] * qrs_c->size[1];
        qrs_c->size[0] = 1;
        qrs_c->size[1] = n;
        emxEnsureCapacity((emxArray__common *)qrs_c, i0, (int)sizeof(double));
        if (n > 0)
        {
            qrs_c->data[0] = anew;
            if (n > 1)
            {
                qrs_c->data[n - 1] = apnd;
                nm1d2 = (n - 1) / 2;
                for (nx = 1; nx < nm1d2; nx++)
                {
                    qrs_c->data[nx] = anew + (double)nx;
                    qrs_c->data[(n - nx) - 1] = apnd - (double)nx;
                }

                if (nm1d2 << 1 == n - 1)
                {
                    qrs_c->data[nm1d2] = (anew + apnd) / 2.0;
                }
                else
                {
                    qrs_c->data[nm1d2] = anew + (double)nm1d2;
                    qrs_c->data[nm1d2 + 1] = apnd - (double)nm1d2;
                }
            }
        }

        /* 'pqrstDetection_22nov2018:338' if all(a < 6250) */
        i0 = d_qrs_c->size[0] * d_qrs_c->size[1];
        d_qrs_c->size[0] = 1;
        d_qrs_c->size[1] = qrs_c->size[1];
        emxEnsureCapacity((emxArray__common *)d_qrs_c, i0, (int)sizeof(boolean_T));
        nm1d2 = qrs_c->size[0] * qrs_c->size[1];
        for (i0 = 0; i0 < nm1d2; i0++)
        {
            d_qrs_c->data[i0] = (qrs_c->data[i0] < 6250.0);
        }

        if (all(d_qrs_c))
        {
            /* 'pqrstDetection_22nov2018:339' m=min(y1(a)); */
            ixstart = 1;
            i0 = i_qrs_c->size[0];
            i_qrs_c->size[0] = qrs_c->size[1];
            emxEnsureCapacity((emxArray__common *)i_qrs_c, i0, (int)sizeof(int));
            nm1d2 = qrs_c->size[1];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                i_qrs_c->data[i0] = (int)qrs_c->data[qrs_c->size[0] * i0];
            }

            n = i_qrs_c->size[0];
            ndbl = ecg_bpass[(int)qrs_c->data[0] - 1];
            i0 = j_qrs_c->size[0];
            j_qrs_c->size[0] = qrs_c->size[1];
            emxEnsureCapacity((emxArray__common *)j_qrs_c, i0, (int)sizeof(int));
            nm1d2 = qrs_c->size[1];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                j_qrs_c->data[i0] = (int)qrs_c->data[qrs_c->size[0] * i0];
            }

            if (j_qrs_c->size[0] > 1)
            {
                if (rtIsNaN(ndbl))
                {
                    idx = 2;
                    exitg6 = false;
                    while ((!exitg6) && (idx <= n))
                    {
                        ixstart = idx;
                        if (!rtIsNaN(ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * (idx - 1)]
                                               - 1]))
                        {
                            ndbl = ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * (idx - 1)] - 1];
                            exitg6 = true;
                        }
                        else
                        {
                            idx++;
                        }
                    }
                }

                i0 = k_qrs_c->size[0];
                k_qrs_c->size[0] = qrs_c->size[1];
                emxEnsureCapacity((emxArray__common *)k_qrs_c, i0, (int)sizeof(int));
                nm1d2 = qrs_c->size[1];
                for (i0 = 0; i0 < nm1d2; i0++)
                {
                    k_qrs_c->data[i0] = (int)qrs_c->data[qrs_c->size[0] * i0];
                }

                if (ixstart < k_qrs_c->size[0])
                {
                    while (ixstart + 1 <= n)
                    {
                        if (ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * ixstart] - 1] < ndbl)
                        {
                            ndbl = ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * ixstart] - 1];
                        }

                        ixstart++;
                    }
                }
            }

            /* 'pqrstDetection_22nov2018:341' b=find(y1(a)==m); */
            i0 = x->size[0];
            x->size[0] = qrs_c->size[1];
            emxEnsureCapacity((emxArray__common *)x, i0, (int)sizeof(boolean_T));
            nm1d2 = qrs_c->size[1];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                x->data[i0] = (ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * i0] - 1] ==
                               ndbl);
            }

            nx = x->size[0];
            idx = 0;
            i0 = ii->size[0];
            ii->size[0] = x->size[0];
            emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
            nm1d2 = 1;
            exitg5 = false;
            while ((!exitg5) && (nm1d2 <= nx))
            {
                guard3 = false;
                if (x->data[nm1d2 - 1])
                {
                    idx++;
                    ii->data[idx - 1] = nm1d2;
                    if (idx >= nx)
                    {
                        exitg5 = true;
                    }
                    else
                    {
                        guard3 = true;
                    }
                }
                else
                {
                    guard3 = true;
                }

                if (guard3)
                {
                    nm1d2++;
                }
            }

            if (x->size[0] == 1)
            {
                if (idx == 0)
                {
                    i0 = ii->size[0];
                    ii->size[0] = 0;
                    emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
                }
            }
            else
            {
                i0 = ii->size[0];
                if (1 > idx)
                {
                    ii->size[0] = 0;
                }
                else
                {
                    ii->size[0] = idx;
                }

                emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
            }

            i0 = ecg_m->size[0];
            ecg_m->size[0] = ii->size[0];
            emxEnsureCapacity((emxArray__common *)ecg_m, i0, (int)sizeof(double));
            nm1d2 = ii->size[0];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                ecg_m->data[i0] = ii->data[i0];
            }

            /* 'pqrstDetection_22nov2018:343' b=b(1); */
            /* 'pqrstDetection_22nov2018:345' b=a(b); */
            /* 'pqrstDetection_22nov2018:347' Qloc(i,j)=b; */
            NOISL_buf->data[NOISL_buf->size[0] * j] = qrs_c->data[(int)ecg_m->data[0]
                                                                  - 1];

            /* 'pqrstDetection_22nov2018:349' Qamp(i,j)=m; */
            SIGL_buf1->data[SIGL_buf1->size[0] * j] = ndbl;

            /* disp(Qamp); */
        }
        else
        {
            /* 'pqrstDetection_22nov2018:352' else */
            /* 'pqrstDetection_22nov2018:353' fprintf('Not able to detect q wave\n'); */
            b_cfprintf();

            /* check if first q wave value is less than 50 or last q wave value is zero */
        }

        /* With similar logic you  can detect the S and T peaks. */
        /*         %% S  Detection */
        /* 'pqrstDetection_22nov2018:360' a=(Rloc(i,j)+5:Rloc(i,j)+30); */
        if (rtIsNaN(qrs_i->data[qrs_i->size[0] * j] + 5.0) || rtIsNaN(qrs_i->
                data[qrs_i->size[0] * j] + 30.0))
        {
            n = 1;
            anew = rtNaN;
            apnd = qrs_i->data[qrs_i->size[0] * j] + 30.0;
        }
        else if (qrs_i->data[qrs_i->size[0] * j] + 30.0 < qrs_i->data[qrs_i->size
                                                                      [0] * j] + 5.0)
        {
            n = 0;
            anew = qrs_i->data[qrs_i->size[0] * j] + 5.0;
            apnd = qrs_i->data[qrs_i->size[0] * j] + 30.0;
        }
        else if (rtIsInf(qrs_i->data[qrs_i->size[0] * j] + 5.0) || rtIsInf
                (qrs_i->data[qrs_i->size[0] * j] + 30.0))
        {
            n = 1;
            anew = rtNaN;
            apnd = qrs_i->data[qrs_i->size[0] * j] + 30.0;
        }
        else
        {
            anew = qrs_i->data[qrs_i->size[0] * j] + 5.0;
            ndbl = floor(((qrs_i->data[qrs_i->size[0] * j] + 30.0) - (qrs_i->
                    data[qrs_i->size[0] * j] + 5.0)) + 0.5);
            apnd = (qrs_i->data[qrs_i->size[0] * j] + 5.0) + ndbl;
            cdiff = apnd - (qrs_i->data[qrs_i->size[0] * j] + 30.0);
            absa = fabs(qrs_i->data[qrs_i->size[0] * j] + 5.0);
            absb = fabs(qrs_i->data[qrs_i->size[0] * j] + 30.0);
            if ((absa >= absb) || rtIsNaN(absb))
            {
                absb = absa;
            }

            if (fabs(cdiff) < 4.4408920985006262E-16 * absb)
            {
                ndbl++;
                apnd = qrs_i->data[qrs_i->size[0] * j] + 30.0;
            }
            else if (cdiff > 0.0)
            {
                apnd = (qrs_i->data[qrs_i->size[0] * j] + 5.0) + (ndbl - 1.0);
            }
            else
            {
                ndbl++;
            }

            if (ndbl >= 0.0)
            {
                n = (int)ndbl;
            }
            else
            {
                n = 0;
            }
        }

        i0 = qrs_c->size[0] * qrs_c->size[1];
        qrs_c->size[0] = 1;
        qrs_c->size[1] = n;
        emxEnsureCapacity((emxArray__common *)qrs_c, i0, (int)sizeof(double));
        if (n > 0)
        {
            qrs_c->data[0] = anew;
            if (n > 1)
            {
                qrs_c->data[n - 1] = apnd;
                nm1d2 = (n - 1) / 2;
                for (nx = 1; nx < nm1d2; nx++)
                {
                    qrs_c->data[nx] = anew + (double)nx;
                    qrs_c->data[(n - nx) - 1] = apnd - (double)nx;
                }

                if (nm1d2 << 1 == n - 1)
                {
                    qrs_c->data[nm1d2] = (anew + apnd) / 2.0;
                }
                else
                {
                    qrs_c->data[nm1d2] = anew + (double)nm1d2;
                    qrs_c->data[nm1d2 + 1] = apnd - (double)nm1d2;
                }
            }
        }

        /* 'pqrstDetection_22nov2018:362' if all(a < 6250) */
        i0 = c_qrs_c->size[0] * c_qrs_c->size[1];
        c_qrs_c->size[0] = 1;
        c_qrs_c->size[1] = qrs_c->size[1];
        emxEnsureCapacity((emxArray__common *)c_qrs_c, i0, (int)sizeof(boolean_T));
        nm1d2 = qrs_c->size[0] * qrs_c->size[1];
        for (i0 = 0; i0 < nm1d2; i0++)
        {
            c_qrs_c->data[i0] = (qrs_c->data[i0] < 6250.0);
        }

        if (all(c_qrs_c))
        {
            /* 'pqrstDetection_22nov2018:363' m=min(y1(a)); */
            ixstart = 1;
            i0 = l_qrs_c->size[0];
            l_qrs_c->size[0] = qrs_c->size[1];
            emxEnsureCapacity((emxArray__common *)l_qrs_c, i0, (int)sizeof(int));
            nm1d2 = qrs_c->size[1];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                l_qrs_c->data[i0] = (int)qrs_c->data[qrs_c->size[0] * i0];
            }

            n = l_qrs_c->size[0];
            ndbl = ecg_bpass[(int)qrs_c->data[0] - 1];
            i0 = m_qrs_c->size[0];
            m_qrs_c->size[0] = qrs_c->size[1];
            emxEnsureCapacity((emxArray__common *)m_qrs_c, i0, (int)sizeof(int));
            nm1d2 = qrs_c->size[1];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                m_qrs_c->data[i0] = (int)qrs_c->data[qrs_c->size[0] * i0];
            }

            if (m_qrs_c->size[0] > 1)
            {
                if (rtIsNaN(ndbl))
                {
                    idx = 2;
                    exitg4 = false;
                    while ((!exitg4) && (idx <= n))
                    {
                        ixstart = idx;
                        if (!rtIsNaN(ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * (idx - 1)]
                                               - 1]))
                        {
                            ndbl = ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * (idx - 1)] - 1];
                            exitg4 = true;
                        }
                        else
                        {
                            idx++;
                        }
                    }
                }

                i0 = n_qrs_c->size[0];
                n_qrs_c->size[0] = qrs_c->size[1];
                emxEnsureCapacity((emxArray__common *)n_qrs_c, i0, (int)sizeof(int));
                nm1d2 = qrs_c->size[1];
                for (i0 = 0; i0 < nm1d2; i0++)
                {
                    n_qrs_c->data[i0] = (int)qrs_c->data[qrs_c->size[0] * i0];
                }

                if (ixstart < n_qrs_c->size[0])
                {
                    while (ixstart + 1 <= n)
                    {
                        if (ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * ixstart] - 1] < ndbl)
                        {
                            ndbl = ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * ixstart] - 1];
                        }

                        ixstart++;
                    }
                }
            }

            /* 'pqrstDetection_22nov2018:365' b=find(y1(a)==m); */
            i0 = x->size[0];
            x->size[0] = qrs_c->size[1];
            emxEnsureCapacity((emxArray__common *)x, i0, (int)sizeof(boolean_T));
            nm1d2 = qrs_c->size[1];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                x->data[i0] = (ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * i0] - 1] ==
                               ndbl);
            }

            nx = x->size[0];
            idx = 0;
            i0 = ii->size[0];
            ii->size[0] = x->size[0];
            emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
            nm1d2 = 1;
            exitg3 = false;
            while ((!exitg3) && (nm1d2 <= nx))
            {
                guard2 = false;
                if (x->data[nm1d2 - 1])
                {
                    idx++;
                    ii->data[idx - 1] = nm1d2;
                    if (idx >= nx)
                    {
                        exitg3 = true;
                    }
                    else
                    {
                        guard2 = true;
                    }
                }
                else
                {
                    guard2 = true;
                }

                if (guard2)
                {
                    nm1d2++;
                }
            }

            if (x->size[0] == 1)
            {
                if (idx == 0)
                {
                    i0 = ii->size[0];
                    ii->size[0] = 0;
                    emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
                }
            }
            else
            {
                i0 = ii->size[0];
                if (1 > idx)
                {
                    ii->size[0] = 0;
                }
                else
                {
                    ii->size[0] = idx;
                }

                emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
            }

            i0 = ecg_m->size[0];
            ecg_m->size[0] = ii->size[0];
            emxEnsureCapacity((emxArray__common *)ecg_m, i0, (int)sizeof(double));
            nm1d2 = ii->size[0];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                ecg_m->data[i0] = ii->data[i0];
            }

            /* 'pqrstDetection_22nov2018:367' b=b(1); */
            /* 'pqrstDetection_22nov2018:369' b=a(b); */
            /* 'pqrstDetection_22nov2018:371' Sloc(i,j)=b; */
            NOISL_buf1->data[NOISL_buf1->size[0] * j] = qrs_c->data[(int)ecg_m->data[0]
                                                                    - 1];

            /* 'pqrstDetection_22nov2018:373' Samp(i,j)=m; */
            THRS_buf1->data[THRS_buf1->size[0] * j] = ndbl;

            /* disp(Samp); */
            scount++;
        }
        else
        {
            /* 'pqrstDetection_22nov2018:375' else */
            /* 'pqrstDetection_22nov2018:376' fprintf('Not able to detect s wave\n'); */
            c_cfprintf();

            /* check if f5irst s wave value is less than 50 or last s wave value is zero */
        }

        /*         %% T Peak */
        /* 'pqrstDetection_22nov2018:381' a=(Rloc(i,j)+25:Rloc(i,j)+70); */
        if (rtIsNaN(qrs_i->data[qrs_i->size[0] * j] + 25.0) || rtIsNaN(qrs_i->
                data[qrs_i->size[0] * j] + 70.0))
        {
            n = 1;
            anew = rtNaN;
            apnd = qrs_i->data[qrs_i->size[0] * j] + 70.0;
        }
        else if (qrs_i->data[qrs_i->size[0] * j] + 70.0 < qrs_i->data[qrs_i->size
                                                                      [0] * j] + 25.0)
        {
            n = 0;
            anew = qrs_i->data[qrs_i->size[0] * j] + 25.0;
            apnd = qrs_i->data[qrs_i->size[0] * j] + 70.0;
        }
        else if (rtIsInf(qrs_i->data[qrs_i->size[0] * j] + 25.0) || rtIsInf
                (qrs_i->data[qrs_i->size[0] * j] + 70.0))
        {
            n = 1;
            anew = rtNaN;
            apnd = qrs_i->data[qrs_i->size[0] * j] + 70.0;
        }
        else
        {
            anew = qrs_i->data[qrs_i->size[0] * j] + 25.0;
            ndbl = floor(((qrs_i->data[qrs_i->size[0] * j] + 70.0) - (qrs_i->
                    data[qrs_i->size[0] * j] + 25.0)) + 0.5);
            apnd = (qrs_i->data[qrs_i->size[0] * j] + 25.0) + ndbl;
            cdiff = apnd - (qrs_i->data[qrs_i->size[0] * j] + 70.0);
            absa = fabs(qrs_i->data[qrs_i->size[0] * j] + 25.0);
            absb = fabs(qrs_i->data[qrs_i->size[0] * j] + 70.0);
            if ((absa >= absb) || rtIsNaN(absb))
            {
                absb = absa;
            }

            if (fabs(cdiff) < 4.4408920985006262E-16 * absb)
            {
                ndbl++;
                apnd = qrs_i->data[qrs_i->size[0] * j] + 70.0;
            }
            else if (cdiff > 0.0)
            {
                apnd = (qrs_i->data[qrs_i->size[0] * j] + 25.0) + (ndbl - 1.0);
            }
            else
            {
                ndbl++;
            }

            if (ndbl >= 0.0)
            {
                n = (int)ndbl;
            }
            else
            {
                n = 0;
            }
        }

        i0 = qrs_c->size[0] * qrs_c->size[1];
        qrs_c->size[0] = 1;
        qrs_c->size[1] = n;
        emxEnsureCapacity((emxArray__common *)qrs_c, i0, (int)sizeof(double));
        if (n > 0)
        {
            qrs_c->data[0] = anew;
            if (n > 1)
            {
                qrs_c->data[n - 1] = apnd;
                nm1d2 = (n - 1) / 2;
                for (nx = 1; nx < nm1d2; nx++)
                {
                    qrs_c->data[nx] = anew + (double)nx;
                    qrs_c->data[(n - nx) - 1] = apnd - (double)nx;
                }

                if (nm1d2 << 1 == n - 1)
                {
                    qrs_c->data[nm1d2] = (anew + apnd) / 2.0;
                }
                else
                {
                    qrs_c->data[nm1d2] = anew + (double)nm1d2;
                    qrs_c->data[nm1d2 + 1] = apnd - (double)nm1d2;
                }
            }
        }

        /* 'pqrstDetection_22nov2018:383' if all(a < 6250) */
        i0 = b_qrs_c->size[0] * b_qrs_c->size[1];
        b_qrs_c->size[0] = 1;
        b_qrs_c->size[1] = qrs_c->size[1];
        emxEnsureCapacity((emxArray__common *)b_qrs_c, i0, (int)sizeof(boolean_T));
        nm1d2 = qrs_c->size[0] * qrs_c->size[1];
        for (i0 = 0; i0 < nm1d2; i0++)
        {
            b_qrs_c->data[i0] = (qrs_c->data[i0] < 6250.0);
        }

        if (all(b_qrs_c))
        {
            /* 'pqrstDetection_22nov2018:384' m=max(y1(a)); */
            ixstart = 1;
            i0 = o_qrs_c->size[0];
            o_qrs_c->size[0] = qrs_c->size[1];
            emxEnsureCapacity((emxArray__common *)o_qrs_c, i0, (int)sizeof(int));
            nm1d2 = qrs_c->size[1];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                o_qrs_c->data[i0] = (int)qrs_c->data[qrs_c->size[0] * i0];
            }

            n = o_qrs_c->size[0];
            ndbl = ecg_bpass[(int)qrs_c->data[0] - 1];
            i0 = p_qrs_c->size[0];
            p_qrs_c->size[0] = qrs_c->size[1];
            emxEnsureCapacity((emxArray__common *)p_qrs_c, i0, (int)sizeof(int));
            nm1d2 = qrs_c->size[1];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                p_qrs_c->data[i0] = (int)qrs_c->data[qrs_c->size[0] * i0];
            }

            if (p_qrs_c->size[0] > 1)
            {
                if (rtIsNaN(ndbl))
                {
                    idx = 2;
                    exitg2 = false;
                    while ((!exitg2) && (idx <= n))
                    {
                        ixstart = idx;
                        if (!rtIsNaN(ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * (idx - 1)]
                                               - 1]))
                        {
                            ndbl = ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * (idx - 1)] - 1];
                            exitg2 = true;
                        }
                        else
                        {
                            idx++;
                        }
                    }
                }

                i0 = q_qrs_c->size[0];
                q_qrs_c->size[0] = qrs_c->size[1];
                emxEnsureCapacity((emxArray__common *)q_qrs_c, i0, (int)sizeof(int));
                nm1d2 = qrs_c->size[1];
                for (i0 = 0; i0 < nm1d2; i0++)
                {
                    q_qrs_c->data[i0] = (int)qrs_c->data[qrs_c->size[0] * i0];
                }

                if (ixstart < q_qrs_c->size[0])
                {
                    while (ixstart + 1 <= n)
                    {
                        if (ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * ixstart] - 1] > ndbl)
                        {
                            ndbl = ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * ixstart] - 1];
                        }

                        ixstart++;
                    }
                }
            }

            /* 'pqrstDetection_22nov2018:386' b=find(y1(a)==m); */
            i0 = x->size[0];
            x->size[0] = qrs_c->size[1];
            emxEnsureCapacity((emxArray__common *)x, i0, (int)sizeof(boolean_T));
            nm1d2 = qrs_c->size[1];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                x->data[i0] = (ecg_bpass[(int)qrs_c->data[qrs_c->size[0] * i0] - 1] ==
                               ndbl);
            }

            nx = x->size[0];
            idx = 0;
            i0 = ii->size[0];
            ii->size[0] = x->size[0];
            emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
            nm1d2 = 1;
            exitg1 = false;
            while ((!exitg1) && (nm1d2 <= nx))
            {
                guard1 = false;
                if (x->data[nm1d2 - 1])
                {
                    idx++;
                    ii->data[idx - 1] = nm1d2;
                    if (idx >= nx)
                    {
                        exitg1 = true;
                    }
                    else
                    {
                        guard1 = true;
                    }
                }
                else
                {
                    guard1 = true;
                }

                if (guard1)
                {
                    nm1d2++;
                }
            }

            if (x->size[0] == 1)
            {
                if (idx == 0)
                {
                    i0 = ii->size[0];
                    ii->size[0] = 0;
                    emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
                }
            }
            else
            {
                i0 = ii->size[0];
                if (1 > idx)
                {
                    ii->size[0] = 0;
                }
                else
                {
                    ii->size[0] = idx;
                }

                emxEnsureCapacity((emxArray__common *)ii, i0, (int)sizeof(int));
            }

            i0 = ecg_m->size[0];
            ecg_m->size[0] = ii->size[0];
            emxEnsureCapacity((emxArray__common *)ecg_m, i0, (int)sizeof(double));
            nm1d2 = ii->size[0];
            for (i0 = 0; i0 < nm1d2; i0++)
            {
                ecg_m->data[i0] = ii->data[i0];
            }

            /* 'pqrstDetection_22nov2018:388' b=b(1); */
            /* 'pqrstDetection_22nov2018:390' b=a(b); */
            /* 'pqrstDetection_22nov2018:392' Tloc(i,j)=b+15; */
            THRS_buf->data[THRS_buf->size[0] * j] = qrs_c->data[(int)ecg_m->data[0] -
                                                                1] + 22.0;

            /* 'pqrstDetection_22nov2018:394' Tamp(i,j)=m; */
            Tamp->data[Tamp->size[0] * j] = ndbl;

            /* disp(Tamp); */
        }
        else
        {
            /* 'pqrstDetection_22nov2018:396' else */
            /* 'pqrstDetection_22nov2018:397' fprintf('Not able to detect t wave\n'); */
            d_cfprintf();

            /* check if last w wave value is greater than number of samples */
        }

        j++;
    }

    emxFree_int32_T(&q_qrs_c);
    emxFree_int32_T(&p_qrs_c);
    emxFree_int32_T(&o_qrs_c);
    emxFree_int32_T(&n_qrs_c);
    emxFree_int32_T(&m_qrs_c);
    emxFree_int32_T(&l_qrs_c);
    emxFree_int32_T(&k_qrs_c);
    emxFree_int32_T(&j_qrs_c);
    emxFree_int32_T(&i_qrs_c);
    emxFree_int32_T(&h_qrs_c);
    emxFree_int32_T(&g_qrs_c);
    emxFree_int32_T(&f_qrs_c);
    emxFree_boolean_T(&e_qrs_c);
    emxFree_boolean_T(&d_qrs_c);
    emxFree_boolean_T(&c_qrs_c);
    emxFree_boolean_T(&b_qrs_c);
    emxFree_int32_T(&ii);
    emxFree_boolean_T(&x);
    emxFree_real_T(&qrs_c);
    emxFree_real_T(&ecg_m);


    /* 'pqrstDetection_22nov2018:403' N=length(ecg); */
    /* 'pqrstDetection_22nov2018:404' duration_in_sec=N/fs; */
    /* 'pqrstDetection_22nov2018:405' duration_in_minutes = duration_in_sec/60; */
    /* 'pqrstDetection_22nov2018:406' heart=Beat_C/duration_in_minutes; */
    *heart = (double)Beat_C / (6250.0 / fs / 60.0);

    /* calculate heartrate */
    /* PR interval calculation */
    /* 'pqrstDetection_22nov2018:409' if((length(Ploc) >= 3)&&(Beat_C >=3)) */
    if ((nois_i->size[1] >= 3) && ((int)Beat_C >= 3))
    {
        /* 'pqrstDetection_22nov2018:411' ppeak1=Ploc(1); */
        /* 'pqrstDetection_22nov2018:412' ppeak2=Ploc(2); */
        /* 'pqrstDetection_22nov2018:413' ppeak3=Ploc(3); */
        /* 'pqrstDetection_22nov2018:415' qpeak1=Qloc(1); */
        /* 'pqrstDetection_22nov2018:416' qpeak2=Qloc(2); */
        /* 'pqrstDetection_22nov2018:417' qpeak3=Qloc(3); */
        /* 'pqrstDetection_22nov2018:420' pr_interval=(((qpeak1-ppeak1)+(qpeak2-ppeak2)+(qpeak3-ppeak3))/3); */
        *pr_interval = (((NOISL_buf->data[0] - nois_i->data[0]) + (NOISL_buf->data[1]
                                                                   - nois_i->data[1])) + (NOISL_buf->data[2] - nois_i->data[2])) / 3.0;

        /* +(rpeak3-ppeak3)+(rpeak4-ppeak4))/3); */
        /* 'pqrstDetection_22nov2018:422' fprintf('PR interval is :-> %f msec\n',((pr_interval/fs)*1000)); */
        if(change==true)
        {
            e_cfprintf(*pr_interval / fs * 1000.0);
        }

    }
    else
    {
        /* 'pqrstDetection_22nov2018:424' else */
        /* 'pqrstDetection_22nov2018:425' fprintf('Insufficient number of R and p peaks found!!,there number should be atleast 3\n'); */
        f_cfprintf();
    }

    /* QRS interval Calculation */
    /* 'pqrstDetection_22nov2018:431' if((length(Qloc) >= 3)&&(length(Sloc) >= 3)) */
    if ((NOISL_buf->size[1] >= 3) && (NOISL_buf1->size[1] >= 3))
    {
        /* 'pqrstDetection_22nov2018:433' qpeak1=Qloc(1); */
        /* 'pqrstDetection_22nov2018:434' qpeak2=Qloc(2); */
        /* 'pqrstDetection_22nov2018:435' qpeak3=Qloc(3); */
        /* 'pqrstDetection_22nov2018:437' speak1=Sloc(1); */
        /* 'pqrstDetection_22nov2018:438' speak2=Sloc(2); */
        /* 'pqrstDetection_22nov2018:439' speak3=Sloc(3); */
        /* 'pqrstDetection_22nov2018:441' qrs_interval=abs(((speak1-qpeak1)+(speak2-qpeak2)+(speak3-qpeak3))/3); */
        *qrs_interval = fabs((((NOISL_buf1->data[0] - NOISL_buf->data[0]) +
                               (NOISL_buf1->data[1] - NOISL_buf->data[1])) + (NOISL_buf1->data[2] -
                                                                              NOISL_buf->data[2])) / 3.0);

        /* 'pqrstDetection_22nov2018:442' fprintf('QRS interval is :-> %f msec\n',((qrs_interval/fs)*1000)); */
        if(change==true)
        {
            g_cfprintf(*qrs_interval / fs * 1000.0);
        }

    }
    else
    {
        /* 'pqrstDetection_22nov2018:444' else */
        /* 'pqrstDetection_22nov2018:445' fprintf('Insufficient number of S  and q peaks found!!,there number should be atleast 3\n'); */
        h_cfprintf();
    }

    /* QT interval calculation */
    /* 'pqrstDetection_22nov2018:449' if((length(Tloc) >= 3)&&(length(Qloc) >= 3)) */
    if ((THRS_buf->size[1] >= 3) && (NOISL_buf->size[1] >= 3))
    {
        /* 'pqrstDetection_22nov2018:451' qpeak1=Qloc(1); */
        /* 'pqrstDetection_22nov2018:452' qpeak2=Qloc(2); */
        /* 'pqrstDetection_22nov2018:453' qpeak3=Qloc(3); */
        /* 'pqrstDetection_22nov2018:455' tpeak1=Tloc(1); */
        /* 'pqrstDetection_22nov2018:456' tpeak2=Tloc(2); */
        /* 'pqrstDetection_22nov2018:457' tpeak3=Tloc(3); */
        /* 'pqrstDetection_22nov2018:459' qt_interval=abs(((tpeak1-qpeak1)+(tpeak2-qpeak2)+(tpeak3-qpeak3))/3); */
        *qt_interval = fabs((((THRS_buf->data[0] - NOISL_buf->data[0]) +
                              (THRS_buf->data[1] - NOISL_buf->data[1])) +
                             (THRS_buf->data[2] - NOISL_buf->data[2])) / 3.0);

        /* 'pqrstDetection_22nov2018:460' fprintf('QT interval is :-> %f msec\n',((qt_interval/fs)*1000)); */
        if(change==true)
        {
            i_cfprintf(*qt_interval / fs * 1000.0);
        }

    }
    else
    {
        /* 'pqrstDetection_22nov2018:461' else */
        /* 'pqrstDetection_22nov2018:462' fprintf('Insuffiecient number of Q and T peaks found ,there number should be atleast 3\n'); */
        j_cfprintf();
    }

    /* Heart rate calculation */
    /* 'pqrstDetection_22nov2018:466' if(Beat_C >=3) */
    if ((int)Beat_C >= 3)
    {
        /* 'pqrstDetection_22nov2018:468' fprintf('Heartrate:-> %f bpm\n',heart); */
        *count1=1;
        if(change==true)
        {
            k_cfprintf(*heart);
        }


        /* RR interval calculation */
        /* 'pqrstDetection_22nov2018:471' rpeak1=Rloc(1); */
        /* 'pqrstDetection_22nov2018:472' rpeak2=Rloc(2); */
        /* 'pqrstDetection_22nov2018:473' rpeak3=Rloc(3); */
        /* 'pqrstDetection_22nov2018:475' rr_interval=abs(((rpeak2-rpeak1)+(rpeak3-rpeak2))/2); */
        rr_interval = fabs(((qrs_i->data[1] - qrs_i->data[0]) + (qrs_i->data[2] -
                                                                 qrs_i->data[1])) / 2.0);

        /* 'pqrstDetection_22nov2018:477' fprintf('RR interval is :-> %f msec\n',((rr_interval/fs)*1000)); */
        if(change==true)
        {
            l_cfprintf(rr_interval / fs * 1000.0);
        }

    }
    else
    {
        /* 'pqrstDetection_22nov2018:479' else */
        /* 'pqrstDetection_22nov2018:480' fprintf('Insuffiecient number of R peaks found ,there number should be atleast 3\n'); */
        m_cfprintf();
        *heart=0;
        *count1=0;
    }

    /* QTc interval calculation */
    /* 'pqrstDetection_22nov2018:485' if((length(Rloc) >= 3)&&(length(Qloc) >=3)&&(length(Tloc) >=3)) */
    if ((qrs_i->size[1] >= 3) && (NOISL_buf->size[1] >= 3) && (THRS_buf->size[1] >=
                                                               3))
    {
        /* 'pqrstDetection_22nov2018:486' rr_Qtc=rr_interval/fs; */
        /* 'pqrstDetection_22nov2018:487' qt_Qtc=((qt_interval/fs)*1000); */
        /* qtc_interval=(qt_Qtc/(sqrt(rr_Qtc)));%Bazzett's formula */
        /* 'pqrstDetection_22nov2018:489' qtc_interval=(qt_Qtc/(rr_Qtc.^(1/3))); */
        *qtc_interval = *qt_interval / fs * 1000.0 / rt_powd_snf(rr_interval / fs,
                                                                 0.33333333333333331);

        /* Fridericia's formula */
        /* qtc_interval=1000*(qt_Qtc/1000+0.154*(1-rr_Qtc));%Sagie's formula */
        /* 'pqrstDetection_22nov2018:491' fprintf('QTc interval is :-> %f \n',qtc_interval); */
        if(change==true)
        {
            n_cfprintf(*qtc_interval);
        }

    }
    else
    {
        /* 'pqrstDetection_22nov2018:492' else */
        /* 'pqrstDetection_22nov2018:493' fprintf('Insuffiecient number of R ,Q and T peaks found ,there number should be atleast 3\n'); */
        o_cfprintf();
    }

    /* ST interval calculation */
    /* 'pqrstDetection_22nov2018:497' if((length(Sloc) >=3)&&(length(Tloc) >=3)) */
    if ((NOISL_buf1->size[1] >= 3) && (THRS_buf->size[1] >= 3))
    {
        /* 'pqrstDetection_22nov2018:499' speak1=Sloc(1); */
        /* 'pqrstDetection_22nov2018:500' speak2=Sloc(2); */
        /* 'pqrstDetection_22nov2018:501' speak3=Sloc(3); */
        /* 'pqrstDetection_22nov2018:503' tpeak1=Tloc(1); */
        /* 'pqrstDetection_22nov2018:504' tpeak2=Tloc(2); */
        /* 'pqrstDetection_22nov2018:505' tpeak3=Tloc(3); */
        /* 'pqrstDetection_22nov2018:507' st_interval=abs(((tpeak1-speak1)+(tpeak2-speak2)+(tpeak3-speak3))/3); */
        *st_interval = fabs((((THRS_buf->data[0] - NOISL_buf1->data[0]) +
                              (THRS_buf->data[1] - NOISL_buf1->data[1])) +
                             (THRS_buf->data[2] - NOISL_buf1->data[2])) / 3.0);

        /* 'pqrstDetection_22nov2018:508' fprintf('ST interval is :-> %f msec\n',((st_interval/fs)*1000)); */
        if(change==true)
        {
            p_cfprintf(*st_interval / fs * 1000.0);
        }

    }
    else
    {
        /* 'pqrstDetection_22nov2018:509' else */
        /* 'pqrstDetection_22nov2018:510' fprintf('Insuffiecient number of S and T peaks found ,there number should be atleast 3\n');*/

        q_cfprintf();


    }

    /****Calculating difference of R and S peak for lead1,lead2,lead3,aVR,aVL,aVF channels******/

    double RS[Beat_C];
    /* Amplitude detection */
    /* 'pqrstDetection_22nov2018:515' if((length(Ploc) >=3)&&(length(Qloc)>=3)&&(length(Rloc)>=3)&&(length(Sloc)>=3)&&(length(Tloc)>=3)) */

    if ((qrs_i->size[1] >=3) && (NOISL_buf1->size[1] >= 3) )
    {
        for(int i=0; i<(Beat_C-2); i++)
        {
            RS[i]=(double)(nois_c->data[i] + THRS_buf1->data[i]);
            //printf("RS difference[%d]:->%lf\n",i,RS[i]);
            //counter++;
        }

        double minimum=RS[0];
        int location=1;

        for(int i=0; i<Beat_C-2; i++)
        {
            if (RS[i] < minimum)
            {
                minimum = RS[i];
                location = i;
                //printf("Minimum:->%lf\n",minimum);
                *ramp=nois_c->data[i];
                *samp=THRS_buf1->data[i];
                //printf("R amplitude:->%lf\n",*ramp);
                //printf("S amplitude:->%lf\n",*samp);
            }
        }

        //printf("Minimum element is present at location %d and it's value is %f.\n", location, minimum);
        *minVal=minimum;
        /* 'pqrstDetection_22nov2018:516' pamp1=Pamp(1); */
        /* 'pqrstDetection_22nov2018:517' qamp1=Qamp(1); */
        /* 'pqrstDetection_22nov2018:518' ramp1=Ramp(1); */
        /* 'pqrstDetection_22nov2018:519' samp1=Samp(1); */
        /* 'pqrstDetection_22nov2018:520' tamp1=Tamp(1); */
        /* 'pqrstDetection_22nov2018:521' fprintf('Amplitude p wave is:-> %f mVolt\n',pamp1); */
        //r_cfprintf(SIGL_buf->data[0]);

        /* 'pqrstDetection_22nov2018:522' fprintf('Amplitude q wave is:-> %f mVolt\n',qamp1); */
        //s_cfprintf(SIGL_buf1->data[0]);

        /* 'pqrstDetection_22nov2018:523' fprintf('Amplitude r wave is:-> %f mVolt\n',ramp1); */
        //t_cfprintf(nois_c->data[0]);

        /* 'pqrstDetection_22nov2018:524' fprintf('Amplitude s wave is:-> %f mVolt\n',samp1); */
        //u_cfprintf(THRS_buf1->data[0]);

        /* 'pqrstDetection_22nov2018:525' fprintf('Amplitude t wave is:-> %f mVolt\n',tamp1); */
        //v_cfprintf(Tamp->data[0]);
    }
    else
    {
        /* 'pqrstDetection_22nov2018:527' else */
        /* 'pqrstDetection_22nov2018:528' fprintf('Insuffiecient number of R and S peaks found ,there number should be atleast 3\n'); */
        w_cfprintf();
    }

    //printf("Number of r waves:->%d\n",Beat_C);

    //int counter=0;

    /*
    double minimum=RS[0];
    int location=1;

    for(int i=0; i<Beat_C; i++)
    {
        if (RS[i] < minimum)
        {
            minimum = RS[i];
            location = i+1;
            *ramp=nois_c->data[i];
            *samp=THRS_buf1->data[i];
            //printf("R amplitude:->%lf\n",*ramp);
            //printf("S amplitude:->%lf\n",*samp);
        }
    }

    //printf("Minimum element is present at location %d and it's value is %f.\n", location, minimum);
    *minVal=minimum;
    */
    emxFree_real_T(&Tamp);
    emxFree_real_T(&THRS_buf);
    emxFree_real_T(&THRS_buf1);
    emxFree_real_T(&NOISL_buf1);
    emxFree_real_T(&SIGL_buf1);
    emxFree_real_T(&NOISL_buf);
    emxFree_real_T(&SIGL_buf);
    emxFree_real_T(&nois_i);
    emxFree_real_T(&nois_c);
    emxFree_real_T(&qrs_i);

    /* fprintf('not_nois:->%f\n',not_nois); */
    /* fprintf('delay is:->%f\n',delay); */
}

/*
 * File trailer for pqrstDetection_22nov2018.c
 *
 * [EOF]
 */

void b_filtfilt(const double x_in[6250], double y_out[6250])
{
    double d2;
    double d3;
    int i;
    double y[6286];
    double a[6];
    double b_y[6286];
    static const double b_a[6] = { -26.25, -77.5, -122.5, -127.5, -92.5, -36.25 };

    static const double dv3[7] = { 31.25, 51.25, 45.0, 5.0, -35.0, -56.25, -36.25
    };

    double c_y[6286];
    d2 = 2.0 * x_in[0];
    d3 = 2.0 * x_in[6249];
    for (i = 0; i < 18; i++) {
        y[i] = d2 - x_in[18 - i];
    }

    memcpy(&y[18], &x_in[0], 6250U * sizeof(double));
    for (i = 0; i < 18; i++) {
        y[i + 6268] = d3 - x_in[6248 - i];
    }

    for (i = 0; i < 6; i++) {
        a[i] = b_a[i] * y[0];
    }

    memcpy(&b_y[0], &y[0], 6286U * sizeof(double));
    b_filter(dv3, b_y, a, y);
    flipud(y);
    for (i = 0; i < 6; i++) {
        a[i] = b_a[i] * y[0];
    }

    memcpy(&c_y[0], &y[0], 6286U * sizeof(double));
    b_filter(dv3, c_y, a, y);
    flipud(y);
    memcpy(&y_out[0], &y[18], 6250U * sizeof(double));
}

/*
 * Arguments    : const double x_in[6250]
 *                double y_out[6250]
 * Return Type  : void
 */
void filtfilt(const double x_in[6250], double y_out[6250])
{
    double d0;
    double d1;
    int i;
    double y[6286];
    double a[6];
    double b_y[6286];
    static const double b_a[6] = { -0.0039439882316223574, -0.0039439881662410243,
                                   0.0078879763472313544, 0.007887976517609448, -0.0039439882538618879,
                                   -0.0039439882131160871 };

    static const double dv1[7] = { 0.00394398821922609, -0.0, -0.01183196465767827,
                                   -0.0, 0.01183196465767827, -0.0, -0.00394398821922609 };

    static const double dv2[7] = { 1.0, -5.2742756755571953, 11.633009985944149,
                                   -13.744306653246635, 9.17962351329148, -3.2869411038966216,
                                   0.49289054808090393 };

    double c_y[6286];
    d0 = 2.0 * x_in[0];
    d1 = 2.0 * x_in[6249];
    for (i = 0; i < 18; i++) {
        y[i] = d0 - x_in[18 - i];
    }

    memcpy(&y[18], &x_in[0], 6250U * sizeof(double));
    for (i = 0; i < 18; i++) {
        y[i + 6268] = d1 - x_in[6248 - i];
    }

    for (i = 0; i < 6; i++) {
        a[i] = b_a[i] * y[0];
    }

    memcpy(&b_y[0], &y[0], 6286U * sizeof(double));
    filter(dv1, dv2, b_y, a, y);
    flipud(y);
    for (i = 0; i < 6; i++) {
        a[i] = b_a[i] * y[0];
    }

    memcpy(&c_y[0], &y[0], 6286U * sizeof(double));
    filter(dv1, dv2, c_y, a, y);
    flipud(y);
    memcpy(&y_out[0], &y[18], 6250U * sizeof(double));
}

void flipud(double x[6286])
{
    int i;
    double xtmp;
    for (i = 0; i < 3143; i++) {
        xtmp = x[i];
        x[i] = x[6285 - i];
        x[6285 - i] = xtmp;
    }
}

void b_abs(const double x[6250], double y[6250])
{
    int k;
    for (k = 0; k < 6250; k++) {
        y[k] = fabs(x[k]);
    }
}


real_T rtInf;
real_T rtMinusInf;
real_T rtNaN;
real32_T rtInfF;
real32_T rtMinusInfF;
real32_T rtNaNF;

/* Function: rt_InitInfAndNaN ==================================================
 * Abstract:
 * Initialize the rtInf, rtMinusInf, and rtNaN needed by the
 * generated code. NaN is initialized as non-signaling. Assumes IEEE.
 */
void rt_InitInfAndNaN(size_t realSize)
{
    (void) (realSize);
    rtNaN = rtGetNaN();
    rtNaNF = rtGetNaNF();
    rtInf = rtGetInf();
    rtInfF = rtGetInfF();
    rtMinusInf = rtGetMinusInf();
    rtMinusInfF = rtGetMinusInfF();
}

/* Function: rtIsInf ==================================================
 * Abstract:
 * Test if value is infinite
 */
boolean_T rtIsInf(real_T value)
{
    return ((value==rtInf || value==rtMinusInf) ? 1U : 0U);
}

/* Function: rtIsInfF =================================================
 * Abstract:
 * Test if single-precision value is infinite
 */
boolean_T rtIsInfF(real32_T value)
{
    return(((value)==rtInfF || (value)==rtMinusInfF) ? 1U : 0U);
}

/* Function: rtIsNaN ==================================================
 * Abstract:
 * Test if value is not a number
 */
boolean_T rtIsNaN(real_T value)
{

#if defined(_MSC_VER) && (_MSC_VER <= 1200)

    return _isnan(value)? TRUE:FALSE;

#else

    return (value!=value)? 1U:0U;

#endif

}

/* Function: rtIsNaNF =================================================
 * Abstract:
 * Test if single-precision value is not a number
 */
boolean_T rtIsNaNF(real32_T value)
{

#if defined(_MSC_VER) && (_MSC_VER <= 1200)

    return _isnan((real_T)value)? true:false;

#else

    return (value!=value)? 1U:0U;

#endif

}

/*
 * File trailer for rt_nonfinite.c
 *
 * [EOF]
 */


/* Function Definitions */

/*
 * Arguments    : emxArray__common *emxArray
 *                int oldNumel
 *                int elementSize
 * Return Type  : void
 */
void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int elementSize)
{
    int newNumel;
    int i;
    void *newData;
    newNumel = 1;
    for (i = 0; i < emxArray->numDimensions; i++) {
        newNumel *= emxArray->size[i];
    }

    if (newNumel > emxArray->allocatedSize) {
        i = emxArray->allocatedSize;
        if (i < 16) {
            i = 16;
        }

        while (i < newNumel) {
            if (i > 1073741823) {
                i = MAX_int32_T;
            } else {
                i <<= 1;
            }
        }

        newData = calloc((unsigned int)i, (unsigned int)elementSize);
        if (emxArray->data != NULL) {
            memcpy(newData, emxArray->data, (unsigned int)(elementSize * oldNumel));
            if (emxArray->canFreeData) {
                free(emxArray->data);
            }
        }

        emxArray->data = newData;
        emxArray->allocatedSize = i;
        emxArray->canFreeData = true;
    }
}

/*
 * Arguments    : emxArray_boolean_T **pEmxArray
 * Return Type  : void
 */
void emxFree_boolean_T(emxArray_boolean_T **pEmxArray)
{
    if (*pEmxArray != (emxArray_boolean_T *)NULL) {
        if (((*pEmxArray)->data != (boolean_T *)NULL) && (*pEmxArray)->canFreeData)
        {
            free((void *)(*pEmxArray)->data);
        }

        free((void *)(*pEmxArray)->size);
        free((void *)*pEmxArray);
        *pEmxArray = (emxArray_boolean_T *)NULL;
    }
}

/*
 * Arguments    : emxArray_int32_T **pEmxArray
 * Return Type  : void
 */
void emxFree_int32_T(emxArray_int32_T **pEmxArray)
{
    if (*pEmxArray != (emxArray_int32_T *)NULL) {
        if (((*pEmxArray)->data != (int *)NULL) && (*pEmxArray)->canFreeData) {
            free((void *)(*pEmxArray)->data);
        }

        free((void *)(*pEmxArray)->size);
        free((void *)*pEmxArray);
        *pEmxArray = (emxArray_int32_T *)NULL;
    }
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 * Return Type  : void
 */
void emxFree_real_T(emxArray_real_T **pEmxArray)
{
    if (*pEmxArray != (emxArray_real_T *)NULL) {
        if (((*pEmxArray)->data != (double *)NULL) && (*pEmxArray)->canFreeData) {
            free((void *)(*pEmxArray)->data);
        }

        free((void *)(*pEmxArray)->size);
        free((void *)*pEmxArray);
        *pEmxArray = (emxArray_real_T *)NULL;
    }
}

/*
 * Arguments    : emxArray_boolean_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions)
{
    emxArray_boolean_T *emxArray;
    int i;
    *pEmxArray = (emxArray_boolean_T *)malloc(sizeof(emxArray_boolean_T));
    emxArray = *pEmxArray;
    emxArray->data = (boolean_T *)NULL;
    emxArray->numDimensions = numDimensions;
    emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
    emxArray->allocatedSize = 0;
    emxArray->canFreeData = true;
    for (i = 0; i < numDimensions; i++) {
        emxArray->size[i] = 0;
    }
}

/*
 * Arguments    : emxArray_boolean_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInit_boolean_T1(emxArray_boolean_T **pEmxArray, int numDimensions)
{
    emxArray_boolean_T *emxArray;
    int i;
    *pEmxArray = (emxArray_boolean_T *)malloc(sizeof(emxArray_boolean_T));
    emxArray = *pEmxArray;
    emxArray->data = (boolean_T *)NULL;
    emxArray->numDimensions = numDimensions;
    emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
    emxArray->allocatedSize = 0;
    emxArray->canFreeData = true;
    for (i = 0; i < numDimensions; i++) {
        emxArray->size[i] = 0;
    }
}

/*
 * Arguments    : emxArray_int32_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
{
    emxArray_int32_T *emxArray;
    int i;
    *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
    emxArray = *pEmxArray;
    emxArray->data = (int *)NULL;
    emxArray->numDimensions = numDimensions;
    emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
    emxArray->allocatedSize = 0;
    emxArray->canFreeData = true;
    for (i = 0; i < numDimensions; i++) {
        emxArray->size[i] = 0;
    }
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
    emxArray_real_T *emxArray;
    int i;
    *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
    emxArray = *pEmxArray;
    emxArray->data = (double *)NULL;
    emxArray->numDimensions = numDimensions;
    emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
    emxArray->allocatedSize = 0;
    emxArray->canFreeData = true;
    for (i = 0; i < numDimensions; i++) {
        emxArray->size[i] = 0;
    }
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 *                int numDimensions
 * Return Type  : void
 */
void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions)
{
    emxArray_real_T *emxArray;
    int i;
    *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
    emxArray = *pEmxArray;
    emxArray->data = (double *)NULL;
    emxArray->numDimensions = numDimensions;
    emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
    emxArray->allocatedSize = 0;
    emxArray->canFreeData = true;
    for (i = 0; i < numDimensions; i++) {
        emxArray->size[i] = 0;
    }
}

/*
 * File trailer for pqrstDetection_22nov2018_emxutil.c
 *
 * [EOF]
 */

/* Function Definitions */

/*
 * Arguments    : const double a[6250]
 *                double y[6250]
 * Return Type  : void
 */
void power(const double a[6250], double y[6250])
{
    int k;
    for (k = 0; k < 6250; k++) {
        y[k] = a[k] * a[k];
    }
}

/*
 * File trailer for power.c
 *
 * [EOF]
 */

/* Function Definitions */

/*
 * Arguments    : const double x[8]
 * Return Type  : double
 */
double b_mean(const double x[8])
{
    double y;
    int k;
    y = x[0];
    for (k = 0; k < 7; k++) {
        y += x[k + 1];
    }

    y /= 8.0;
    return y;
}

/*
 * Arguments    : const emxArray_real_T *x
 * Return Type  : double
 */
double mean(const emxArray_real_T *x)
{
    double y;
    int k;
    if (x->size[0] == 0) {
        y = 0.0;
    } else {
        y = x->data[0];
        for (k = 2; k <= x->size[0]; k++) {
            y += x->data[k - 1];
        }
    }

    y /= (double)x->size[0];
    return y;
}

/*
 * File trailer for mean.c
 *
 * [EOF]
 */
/* Function Definitions */

/*
 * Arguments    : const double A[6250]
 *                const emxArray_real_T *B
 *                emxArray_real_T *C
 * Return Type  : void
 */
void conv(const double A[6250], const emxArray_real_T *B, emxArray_real_T *C)
{
    int nB;
    int nC;
    int jA2;
    int jC;
    double s;
    int k;
    nB = B->size[1];
    if (B->size[1] == 0) {
        nC = 6250;
    } else {
        nC = B->size[1] + 6249;
    }

    jA2 = C->size[0];
    C->size[0] = nC;
    emxEnsureCapacity((emxArray__common *)C, jA2, (int)sizeof(double));
    if (B->size[1] == 0) {
        jC = C->size[0];
        jA2 = C->size[0];
        C->size[0] = jC;
        emxEnsureCapacity((emxArray__common *)C, jA2, (int)sizeof(double));
        for (jA2 = 0; jA2 < jC; jA2++) {
            C->data[jA2] = 0.0;
        }
    } else {
        for (jC = 1; jC <= nC; jC++) {
            if (6250 <= jC) {
                jA2 = 6250;
            } else {
                jA2 = jC;
            }

            s = 0.0;
            if (nB < jC + 1) {
                k = jC - nB;
            } else {
                k = 0;
            }

            while (k + 1 <= jA2) {
                s += A[k] * B->data[(jC - k) - 1];
                k++;
            }

            C->data[jC - 1] = s;
        }
    }
}

/*
 * File trailer for conv.c
 *
 * [EOF]
 */


/* Function Declarations */
static void findLocalMaxima(emxArray_real_T *yTemp, emxArray_real_T *iPk,
                            emxArray_real_T *iInflect);
static void parse_inputs(const emxArray_real_T *Yin, double varargin_2,
                         emxArray_real_T *y, boolean_T *yIsRow, emxArray_real_T *x, boolean_T *xIsRow,
                         double *Pd, double *NpOut);

/* Function Definitions */

/*
 * Arguments    : emxArray_real_T *yTemp
 *                emxArray_real_T *iPk
 *                emxArray_real_T *iInflect
 * Return Type  : void
 */
static void findLocalMaxima(emxArray_real_T *yTemp, emxArray_real_T *iPk,
                            emxArray_real_T *iInflect)
{
    emxArray_real_T *r3;
    int i3;
    int cdiff;
    int ndbl;
    int apnd;
    int absb;
    emxArray_real_T *y;
    emxArray_real_T *iTemp;
    emxArray_boolean_T *yFinite;
    emxArray_boolean_T *x;
    emxArray_int32_T *ii;
    boolean_T exitg3;
    boolean_T guard3 = false;
    emxArray_int32_T *r4;
    emxArray_real_T *b_iTemp;
    emxArray_real_T *b_yTemp;
    emxArray_real_T *s;
    emxArray_real_T *r5;
    double b_x;
    boolean_T exitg2;
    boolean_T guard2 = false;
    emxArray_int32_T *b_ii;
    boolean_T exitg1;
    boolean_T guard1 = false;
    emxInit_real_T(&r3, 1);
    i3 = r3->size[0];
    r3->size[0] = 2 + yTemp->size[0];
    emxEnsureCapacity((emxArray__common *)r3, i3, (int)sizeof(double));
    r3->data[0] = rtNaN;
    cdiff = yTemp->size[0];
    for (i3 = 0; i3 < cdiff; i3++) {
        r3->data[i3 + 1] = yTemp->data[i3];
    }

    r3->data[1 + yTemp->size[0]] = rtNaN;
    i3 = yTemp->size[0];
    yTemp->size[0] = r3->size[0];
    emxEnsureCapacity((emxArray__common *)yTemp, i3, (int)sizeof(double));
    cdiff = r3->size[0];
    for (i3 = 0; i3 < cdiff; i3++) {
        yTemp->data[i3] = r3->data[i3];
    }

    emxFree_real_T(&r3);
    ndbl = (int)floor(((double)yTemp->size[0] - 1.0) + 0.5);
    apnd = ndbl + 1;
    cdiff = (ndbl - yTemp->size[0]) + 1;
    absb = yTemp->size[0];
    if (abs(cdiff) < 4.4408920985006262E-16 * (double)absb) {
        ndbl++;
        apnd = yTemp->size[0];
    } else if (cdiff > 0) {
        apnd = ndbl;
    } else {
        ndbl++;
    }

    emxInit_real_T1(&y, 2);
    i3 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = ndbl;
    emxEnsureCapacity((emxArray__common *)y, i3, (int)sizeof(double));
    y->data[0] = 1.0;
    if (ndbl > 1) {
        y->data[ndbl - 1] = apnd;
        cdiff = (ndbl - 1) / 2;
        for (absb = 1; absb < cdiff; absb++) {
            y->data[absb] = 1.0 + (double)absb;
            y->data[(ndbl - absb) - 1] = apnd - absb;
        }

        if (cdiff << 1 == ndbl - 1) {
            y->data[cdiff] = (1.0 + (double)apnd) / 2.0;
        } else {
            y->data[cdiff] = 1.0 + (double)cdiff;
            y->data[cdiff + 1] = apnd - cdiff;
        }
    }

    emxInit_real_T(&iTemp, 1);
    i3 = iTemp->size[0];
    iTemp->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)iTemp, i3, (int)sizeof(double));
    cdiff = y->size[1];
    for (i3 = 0; i3 < cdiff; i3++) {
        iTemp->data[i3] = y->data[y->size[0] * i3];
    }

    emxFree_real_T(&y);
    emxInit_boolean_T(&yFinite, 1);
    i3 = yFinite->size[0];
    yFinite->size[0] = yTemp->size[0];
    emxEnsureCapacity((emxArray__common *)yFinite, i3, (int)sizeof(boolean_T));
    cdiff = yTemp->size[0];
    for (i3 = 0; i3 < cdiff; i3++) {
        yFinite->data[i3] = rtIsNaN(yTemp->data[i3]);
    }

    i3 = yFinite->size[0];
    emxEnsureCapacity((emxArray__common *)yFinite, i3, (int)sizeof(boolean_T));
    cdiff = yFinite->size[0];
    for (i3 = 0; i3 < cdiff; i3++) {
        yFinite->data[i3] = !yFinite->data[i3];
    }

    emxInit_boolean_T(&x, 1);
    cdiff = yTemp->size[0] - 2;
    i3 = x->size[0];
    x->size[0] = cdiff + 1;
    emxEnsureCapacity((emxArray__common *)x, i3, (int)sizeof(boolean_T));
    for (i3 = 0; i3 <= cdiff; i3++) {
        x->data[i3] = ((yTemp->data[i3] != yTemp->data[1 + i3]) && (yFinite->data[i3]
                                                                    || yFinite->data[1 + i3]));
    }

    emxFree_boolean_T(&yFinite);
    emxInit_int32_T(&ii, 1);
    ndbl = x->size[0];
    absb = 0;
    i3 = ii->size[0];
    ii->size[0] = x->size[0];
    emxEnsureCapacity((emxArray__common *)ii, i3, (int)sizeof(int));
    cdiff = 1;
    exitg3 = false;
    while ((!exitg3) && (cdiff <= ndbl)) {
        guard3 = false;
        if (x->data[cdiff - 1]) {
            absb++;
            ii->data[absb - 1] = cdiff;
            if (absb >= ndbl) {
                exitg3 = true;
            } else {
                guard3 = true;
            }
        } else {
            guard3 = true;
        }

        if (guard3) {
            cdiff++;
        }
    }

    if (x->size[0] == 1) {
        if (absb == 0) {
            i3 = ii->size[0];
            ii->size[0] = 0;
            emxEnsureCapacity((emxArray__common *)ii, i3, (int)sizeof(int));
        }
    } else {
        i3 = ii->size[0];
        if (1 > absb) {
            ii->size[0] = 0;
        } else {
            ii->size[0] = absb;
        }

        emxEnsureCapacity((emxArray__common *)ii, i3, (int)sizeof(int));
    }

    emxInit_int32_T(&r4, 1);
    i3 = r4->size[0];
    r4->size[0] = 1 + ii->size[0];
    emxEnsureCapacity((emxArray__common *)r4, i3, (int)sizeof(int));
    r4->data[0] = 1;
    cdiff = ii->size[0];
    for (i3 = 0; i3 < cdiff; i3++) {
        r4->data[i3 + 1] = ii->data[i3] + 1;
    }

    emxInit_real_T(&b_iTemp, 1);
    i3 = b_iTemp->size[0];
    b_iTemp->size[0] = r4->size[0];
    emxEnsureCapacity((emxArray__common *)b_iTemp, i3, (int)sizeof(double));
    cdiff = r4->size[0];
    for (i3 = 0; i3 < cdiff; i3++) {
        b_iTemp->data[i3] = iTemp->data[r4->data[i3] - 1];
    }

    emxFree_int32_T(&r4);
    i3 = iTemp->size[0];
    iTemp->size[0] = b_iTemp->size[0];
    emxEnsureCapacity((emxArray__common *)iTemp, i3, (int)sizeof(double));
    cdiff = b_iTemp->size[0];
    for (i3 = 0; i3 < cdiff; i3++) {
        iTemp->data[i3] = b_iTemp->data[i3];
    }

    emxFree_real_T(&b_iTemp);
    emxInit_real_T(&b_yTemp, 1);
    i3 = b_yTemp->size[0];
    b_yTemp->size[0] = iTemp->size[0];
    emxEnsureCapacity((emxArray__common *)b_yTemp, i3, (int)sizeof(double));
    cdiff = iTemp->size[0];
    for (i3 = 0; i3 < cdiff; i3++) {
        b_yTemp->data[i3] = yTemp->data[(int)iTemp->data[i3] - 1];
    }

    emxInit_real_T(&s, 1);
    diff(b_yTemp, s);
    ndbl = s->size[0];
    absb = 0;
    emxFree_real_T(&b_yTemp);
    while (absb + 1 <= ndbl) {
        if (s->data[absb] < 0.0) {
            b_x = -1.0;
        } else if (s->data[absb] > 0.0) {
            b_x = 1.0;
        } else if (s->data[absb] == 0.0) {
            b_x = 0.0;
        } else {
            b_x = s->data[absb];
        }

        s->data[absb] = b_x;
        absb++;
    }

    emxInit_real_T(&r5, 1);
    diff(s, r5);
    i3 = x->size[0];
    x->size[0] = r5->size[0];
    emxEnsureCapacity((emxArray__common *)x, i3, (int)sizeof(boolean_T));
    cdiff = r5->size[0];
    for (i3 = 0; i3 < cdiff; i3++) {
        x->data[i3] = (r5->data[i3] < 0.0);
    }

    emxFree_real_T(&r5);
    ndbl = x->size[0];
    absb = 0;
    i3 = ii->size[0];
    ii->size[0] = x->size[0];
    emxEnsureCapacity((emxArray__common *)ii, i3, (int)sizeof(int));
    cdiff = 1;
    exitg2 = false;
    while ((!exitg2) && (cdiff <= ndbl)) {
        guard2 = false;
        if (x->data[cdiff - 1]) {
            absb++;
            ii->data[absb - 1] = cdiff;
            if (absb >= ndbl) {
                exitg2 = true;
            } else {
                guard2 = true;
            }
        } else {
            guard2 = true;
        }

        if (guard2) {
            cdiff++;
        }
    }

    if (x->size[0] == 1) {
        if (absb == 0) {
            i3 = ii->size[0];
            ii->size[0] = 0;
            emxEnsureCapacity((emxArray__common *)ii, i3, (int)sizeof(int));
        }
    } else {
        i3 = ii->size[0];
        if (1 > absb) {
            ii->size[0] = 0;
        } else {
            ii->size[0] = absb;
        }

        emxEnsureCapacity((emxArray__common *)ii, i3, (int)sizeof(int));
    }

    if (1 > s->size[0] - 1) {
        cdiff = 0;
    } else {
        cdiff = s->size[0] - 1;
    }

    if (2 > s->size[0]) {
        i3 = 0;
    } else {
        i3 = 1;
    }

    absb = x->size[0];
    x->size[0] = cdiff;
    emxEnsureCapacity((emxArray__common *)x, absb, (int)sizeof(boolean_T));
    for (absb = 0; absb < cdiff; absb++) {
        x->data[absb] = (s->data[absb] != s->data[i3 + absb]);
    }

    emxFree_real_T(&s);
    emxInit_int32_T(&b_ii, 1);
    ndbl = x->size[0];
    absb = 0;
    i3 = b_ii->size[0];
    b_ii->size[0] = x->size[0];
    emxEnsureCapacity((emxArray__common *)b_ii, i3, (int)sizeof(int));
    cdiff = 1;
    exitg1 = false;
    while ((!exitg1) && (cdiff <= ndbl)) {
        guard1 = false;
        if (x->data[cdiff - 1]) {
            absb++;
            b_ii->data[absb - 1] = cdiff;
            if (absb >= ndbl) {
                exitg1 = true;
            } else {
                guard1 = true;
            }
        } else {
            guard1 = true;
        }

        if (guard1) {
            cdiff++;
        }
    }

    if (x->size[0] == 1) {
        if (absb == 0) {
            i3 = b_ii->size[0];
            b_ii->size[0] = 0;
            emxEnsureCapacity((emxArray__common *)b_ii, i3, (int)sizeof(int));
        }
    } else {
        i3 = b_ii->size[0];
        if (1 > absb) {
            b_ii->size[0] = 0;
        } else {
            b_ii->size[0] = absb;
        }

        emxEnsureCapacity((emxArray__common *)b_ii, i3, (int)sizeof(int));
    }

    emxFree_boolean_T(&x);
    i3 = iInflect->size[0];
    iInflect->size[0] = b_ii->size[0];
    emxEnsureCapacity((emxArray__common *)iInflect, i3, (int)sizeof(double));
    cdiff = b_ii->size[0];
    for (i3 = 0; i3 < cdiff; i3++) {
        iInflect->data[i3] = iTemp->data[b_ii->data[i3]] - 1.0;
    }

    emxFree_int32_T(&b_ii);
    i3 = iPk->size[0];
    iPk->size[0] = ii->size[0];
    emxEnsureCapacity((emxArray__common *)iPk, i3, (int)sizeof(double));
    cdiff = ii->size[0];
    for (i3 = 0; i3 < cdiff; i3++) {
        iPk->data[i3] = iTemp->data[ii->data[i3]] - 1.0;
    }

    emxFree_int32_T(&ii);
    emxFree_real_T(&iTemp);
}

/*
 * Arguments    : const emxArray_real_T *Yin
 *                double varargin_2
 *                emxArray_real_T *y
 *                boolean_T *yIsRow
 *                emxArray_real_T *x
 *                boolean_T *xIsRow
 *                double *Pd
 *                double *NpOut
 * Return Type  : void
 */
static void parse_inputs(const emxArray_real_T *Yin, double varargin_2,
                         emxArray_real_T *y, boolean_T *yIsRow, emxArray_real_T *x, boolean_T *xIsRow,
                         double *Pd, double *NpOut)
{
    int cdiff;
    int nm1d2;
    int ndbl;
    int apnd;
    emxArray_real_T *b_y;
    *yIsRow = false;
    cdiff = y->size[0];
    y->size[0] = Yin->size[0];
    emxEnsureCapacity((emxArray__common *)y, cdiff, (int)sizeof(double));
    nm1d2 = Yin->size[0];
    for (cdiff = 0; cdiff < nm1d2; cdiff++) {
        y->data[cdiff] = Yin->data[cdiff];
    }

    *xIsRow = false;
    nm1d2 = Yin->size[0];
    if (nm1d2 < 1) {
        ndbl = 0;
        apnd = 0;
    } else {
        nm1d2 = Yin->size[0];
        ndbl = (int)floor(((double)nm1d2 - 1.0) + 0.5);
        apnd = ndbl + 1;
        nm1d2 = Yin->size[0];
        cdiff = (ndbl - nm1d2) + 1;
        nm1d2 = Yin->size[0];
        if (1 >= nm1d2) {
            nm1d2 = 1;
        }

        if (abs(cdiff) < 4.4408920985006262E-16 * (double)nm1d2) {
            ndbl++;
            apnd = Yin->size[0];
        } else if (cdiff > 0) {
            apnd = ndbl;
        } else {
            ndbl++;
        }
    }

    emxInit_real_T1(&b_y, 2);
    cdiff = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = ndbl;
    emxEnsureCapacity((emxArray__common *)b_y, cdiff, (int)sizeof(double));
    if (ndbl > 0) {
        b_y->data[0] = 1.0;
        if (ndbl > 1) {
            b_y->data[ndbl - 1] = apnd;
            nm1d2 = (ndbl - 1) / 2;
            for (cdiff = 1; cdiff < nm1d2; cdiff++) {
                b_y->data[cdiff] = 1.0 + (double)cdiff;
                b_y->data[(ndbl - cdiff) - 1] = apnd - cdiff;
            }

            if (nm1d2 << 1 == ndbl - 1) {
                b_y->data[nm1d2] = (1.0 + (double)apnd) / 2.0;
            } else {
                b_y->data[nm1d2] = 1.0 + (double)nm1d2;
                b_y->data[nm1d2 + 1] = apnd - nm1d2;
            }
        }
    }

    cdiff = x->size[0];
    x->size[0] = b_y->size[1];
    emxEnsureCapacity((emxArray__common *)x, cdiff, (int)sizeof(double));
    nm1d2 = b_y->size[1];
    for (cdiff = 0; cdiff < nm1d2; cdiff++) {
        x->data[cdiff] = b_y->data[b_y->size[0] * cdiff];
    }

    emxFree_real_T(&b_y);
    *Pd = varargin_2;
    nm1d2 = Yin->size[0];
    *NpOut = nm1d2;
}

/*
 * Arguments    : const emxArray_real_T *Yin
 *                double varargin_2
 *                emxArray_real_T *Ypk
 *                emxArray_real_T *Xpk
 * Return Type  : void
 */
void findpeaks(const emxArray_real_T *Yin, double varargin_2, emxArray_real_T
*Ypk, emxArray_real_T *Xpk)
{
    emxArray_real_T *y;
    emxArray_real_T *x;
    emxArray_boolean_T *idelete;
    boolean_T yIsRow;
    boolean_T xIsRow;
    double minD;
    double maxN;
    int i2;
    int cdiff;
    emxArray_int32_T *ii;
    int nx;
    int idx;
    boolean_T exitg1;
    boolean_T guard1 = false;
    emxArray_real_T *iInfite;
    emxArray_real_T *yTemp;
    emxArray_real_T *iPk;
    emxArray_real_T *locs;
    int ndbl;
    double b_locs;
    emxArray_real_T *b_iPk;
    emxArray_int32_T *ib;
    emxArray_boolean_T *r2;
    emxArray_real_T *b_y;
    int apnd;
    unsigned int unnamed_idx_0;
    double c_locs;
    emxInit_real_T(&y, 1);
    emxInit_real_T(&x, 1);
    emxInit_boolean_T(&idelete, 1);
    parse_inputs(Yin, varargin_2, y, &yIsRow, x, &xIsRow, &minD, &maxN);
    i2 = idelete->size[0];
    idelete->size[0] = y->size[0];
    emxEnsureCapacity((emxArray__common *)idelete, i2, (int)sizeof(boolean_T));
    cdiff = y->size[0];
    for (i2 = 0; i2 < cdiff; i2++) {
        idelete->data[i2] = rtIsInf(y->data[i2]);
    }

    i2 = idelete->size[0];
    emxEnsureCapacity((emxArray__common *)idelete, i2, (int)sizeof(boolean_T));
    cdiff = idelete->size[0];
    for (i2 = 0; i2 < cdiff; i2++) {
        idelete->data[i2] = (idelete->data[i2] && (y->data[i2] > 0.0));
    }

    emxInit_int32_T(&ii, 1);
    nx = idelete->size[0];
    idx = 0;
    i2 = ii->size[0];
    ii->size[0] = idelete->size[0];
    emxEnsureCapacity((emxArray__common *)ii, i2, (int)sizeof(int));
    cdiff = 1;
    exitg1 = false;
    while ((!exitg1) && (cdiff <= nx)) {
        guard1 = false;
        if (idelete->data[cdiff - 1]) {
            idx++;
            ii->data[idx - 1] = cdiff;
            if (idx >= nx) {
                exitg1 = true;
            } else {
                guard1 = true;
            }
        } else {
            guard1 = true;
        }

        if (guard1) {
            cdiff++;
        }
    }

    if (idelete->size[0] == 1) {
        if (idx == 0) {
            i2 = ii->size[0];
            ii->size[0] = 0;
            emxEnsureCapacity((emxArray__common *)ii, i2, (int)sizeof(int));
        }
    } else {
        i2 = ii->size[0];
        if (1 > idx) {
            ii->size[0] = 0;
        } else {
            ii->size[0] = idx;
        }

        emxEnsureCapacity((emxArray__common *)ii, i2, (int)sizeof(int));
    }

    emxInit_real_T(&iInfite, 1);
    i2 = iInfite->size[0];
    iInfite->size[0] = ii->size[0];
    emxEnsureCapacity((emxArray__common *)iInfite, i2, (int)sizeof(double));
    cdiff = ii->size[0];
    for (i2 = 0; i2 < cdiff; i2++) {
        iInfite->data[i2] = ii->data[i2];
    }

    emxInit_real_T(&yTemp, 1);
    i2 = yTemp->size[0];
    yTemp->size[0] = y->size[0];
    emxEnsureCapacity((emxArray__common *)yTemp, i2, (int)sizeof(double));
    cdiff = y->size[0];
    for (i2 = 0; i2 < cdiff; i2++) {
        yTemp->data[i2] = y->data[i2];
    }

    i2 = ii->size[0];
    ii->size[0] = iInfite->size[0];
    emxEnsureCapacity((emxArray__common *)ii, i2, (int)sizeof(int));
    cdiff = iInfite->size[0];
    for (i2 = 0; i2 < cdiff; i2++) {
        ii->data[i2] = (int)iInfite->data[i2];
    }

    cdiff = ii->size[0];
    for (i2 = 0; i2 < cdiff; i2++) {
        yTemp->data[ii->data[i2] - 1] = rtNaN;
    }

    emxInit_real_T(&iPk, 1);
    emxInit_real_T(&locs, 1);
    findLocalMaxima(yTemp, iPk, locs);
    if (!(iPk->size[0] == 0)) {
        nx = iPk->size[0] - 1;
        idx = 0;
        for (ndbl = 0; ndbl <= nx; ndbl++) {
            if (y->data[(int)iPk->data[ndbl] - 1] > rtMinusInf) {
                idx++;
            }
        }

        cdiff = 0;
        for (ndbl = 0; ndbl <= nx; ndbl++) {
            if (y->data[(int)iPk->data[ndbl] - 1] > rtMinusInf) {
                iPk->data[cdiff] = iPk->data[ndbl];
                cdiff++;
            }
        }

        i2 = iPk->size[0];
        iPk->size[0] = idx;
        emxEnsureCapacity((emxArray__common *)iPk, i2, (int)sizeof(double));
    }

    cdiff = iPk->size[0];
    i2 = yTemp->size[0];
    yTemp->size[0] = cdiff;
    emxEnsureCapacity((emxArray__common *)yTemp, i2, (int)sizeof(double));
    for (idx = 0; idx + 1 <= cdiff; idx++) {
        if ((y->data[(int)(iPk->data[idx] - 1.0) - 1] >= y->data[(int)(iPk->data[idx]
                                                                       + 1.0) - 1]) || rtIsNaN(y->data[(int)(iPk->data[idx] + 1.0) - 1])) {
            b_locs = y->data[(int)(iPk->data[idx] - 1.0) - 1];
        } else {
            b_locs = y->data[(int)(iPk->data[idx] + 1.0) - 1];
        }

        yTemp->data[idx] = b_locs;
    }

    nx = iPk->size[0] - 1;
    idx = 0;
    for (ndbl = 0; ndbl <= nx; ndbl++) {
        if (y->data[(int)iPk->data[ndbl] - 1] - yTemp->data[ndbl] >= 0.0) {
            idx++;
        }
    }

    cdiff = 0;
    for (ndbl = 0; ndbl <= nx; ndbl++) {
        if (y->data[(int)iPk->data[ndbl] - 1] - yTemp->data[ndbl] >= 0.0) {
            iPk->data[cdiff] = iPk->data[ndbl];
            cdiff++;
        }
    }

    emxInit_real_T(&b_iPk, 1);
    emxInit_int32_T(&ib, 1);
    i2 = iPk->size[0];
    iPk->size[0] = idx;
    emxEnsureCapacity((emxArray__common *)iPk, i2, (int)sizeof(double));
    do_vectors(iPk, iInfite, b_iPk, ii, ib);
    emxFree_int32_T(&ib);
    emxInit_boolean_T(&r2, 1);
    emxInit_real_T1(&b_y, 2);
    if ((b_iPk->size[0] == 0) || (minD == 0.0)) {
        if (b_iPk->size[0] < 1) {
            ndbl = 0;
            apnd = 0;
        } else {
            ndbl = (int)floor(((double)b_iPk->size[0] - 1.0) + 0.5);
            apnd = ndbl + 1;
            cdiff = (ndbl - b_iPk->size[0]) + 1;
            nx = b_iPk->size[0];
            if (abs(cdiff) < 4.4408920985006262E-16 * (double)nx) {
                ndbl++;
                apnd = b_iPk->size[0];
            } else if (cdiff > 0) {
                apnd = ndbl;
            } else {
                ndbl++;
            }
        }

        i2 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = 1;
        b_y->size[1] = ndbl;
        emxEnsureCapacity((emxArray__common *)b_y, i2, (int)sizeof(double));
        if (ndbl > 0) {
            b_y->data[0] = 1.0;
            if (ndbl > 1) {
                b_y->data[ndbl - 1] = apnd;
                nx = (ndbl - 1) / 2;
                for (idx = 1; idx < nx; idx++) {
                    b_y->data[idx] = 1.0 + (double)idx;
                    b_y->data[(ndbl - idx) - 1] = apnd - idx;
                }

                if (nx << 1 == ndbl - 1) {
                    b_y->data[nx] = (1.0 + (double)apnd) / 2.0;
                } else {
                    b_y->data[nx] = 1.0 + (double)nx;
                    b_y->data[nx + 1] = apnd - nx;
                }
            }
        }

        i2 = yTemp->size[0];
        yTemp->size[0] = b_y->size[1];
        emxEnsureCapacity((emxArray__common *)yTemp, i2, (int)sizeof(double));
        cdiff = b_y->size[1];
        for (i2 = 0; i2 < cdiff; i2++) {
            yTemp->data[i2] = b_y->data[b_y->size[0] * i2];
        }
    } else {
        i2 = locs->size[0];
        locs->size[0] = b_iPk->size[0];
        emxEnsureCapacity((emxArray__common *)locs, i2, (int)sizeof(double));
        cdiff = b_iPk->size[0];
        for (i2 = 0; i2 < cdiff; i2++) {
            locs->data[i2] = x->data[(int)b_iPk->data[i2] - 1];
        }

        i2 = yTemp->size[0];
        yTemp->size[0] = b_iPk->size[0];
        emxEnsureCapacity((emxArray__common *)yTemp, i2, (int)sizeof(double));
        cdiff = b_iPk->size[0];
        for (i2 = 0; i2 < cdiff; i2++) {
            yTemp->data[i2] = y->data[(int)b_iPk->data[i2] - 1];
        }

        sort(yTemp, ii);
        i2 = iInfite->size[0];
        iInfite->size[0] = ii->size[0];
        emxEnsureCapacity((emxArray__common *)iInfite, i2, (int)sizeof(double));
        cdiff = ii->size[0];
        for (i2 = 0; i2 < cdiff; i2++) {
            iInfite->data[i2] = ii->data[i2];
        }

        i2 = yTemp->size[0];
        yTemp->size[0] = iInfite->size[0];
        emxEnsureCapacity((emxArray__common *)yTemp, i2, (int)sizeof(double));
        cdiff = iInfite->size[0];
        for (i2 = 0; i2 < cdiff; i2++) {
            yTemp->data[i2] = locs->data[(int)iInfite->data[i2] - 1];
        }

        unnamed_idx_0 = (unsigned int)iInfite->size[0];
        i2 = idelete->size[0];
        idelete->size[0] = (int)unnamed_idx_0;
        emxEnsureCapacity((emxArray__common *)idelete, i2, (int)sizeof(boolean_T));
        cdiff = (int)unnamed_idx_0;
        for (i2 = 0; i2 < cdiff; i2++) {
            idelete->data[i2] = false;
        }

        for (ndbl = 0; ndbl < iInfite->size[0]; ndbl++) {
            if (!idelete->data[ndbl]) {
                b_locs = locs->data[(int)iInfite->data[ndbl] - 1] - minD;
                c_locs = locs->data[(int)iInfite->data[ndbl] - 1] + minD;
                i2 = r2->size[0];
                r2->size[0] = yTemp->size[0];
                emxEnsureCapacity((emxArray__common *)r2, i2, (int)sizeof(boolean_T));
                cdiff = yTemp->size[0];
                for (i2 = 0; i2 < cdiff; i2++) {
                    r2->data[i2] = ((yTemp->data[i2] >= b_locs) && (yTemp->data[i2] <=
                                                                    c_locs));
                }

                i2 = idelete->size[0];
                emxEnsureCapacity((emxArray__common *)idelete, i2, (int)sizeof(boolean_T));
                cdiff = idelete->size[0];
                for (i2 = 0; i2 < cdiff; i2++) {
                    idelete->data[i2] = (idelete->data[i2] || r2->data[i2]);
                }

                idelete->data[ndbl] = false;
            }
        }

        nx = idelete->size[0] - 1;
        idx = 0;
        for (ndbl = 0; ndbl <= nx; ndbl++) {
            if (!idelete->data[ndbl]) {
                idx++;
            }
        }

        i2 = ii->size[0];
        ii->size[0] = idx;
        emxEnsureCapacity((emxArray__common *)ii, i2, (int)sizeof(int));
        cdiff = 0;
        for (ndbl = 0; ndbl <= nx; ndbl++) {
            if (!idelete->data[ndbl]) {
                ii->data[cdiff] = ndbl + 1;
                cdiff++;
            }
        }

        i2 = yTemp->size[0];
        yTemp->size[0] = ii->size[0];
        emxEnsureCapacity((emxArray__common *)yTemp, i2, (int)sizeof(double));
        cdiff = ii->size[0];
        for (i2 = 0; i2 < cdiff; i2++) {
            yTemp->data[i2] = iInfite->data[ii->data[i2] - 1];
        }

        c_sort(yTemp);
    }

    emxFree_real_T(&b_y);
    emxFree_boolean_T(&r2);
    emxFree_boolean_T(&idelete);
    emxFree_real_T(&locs);
    emxFree_int32_T(&ii);
    emxFree_real_T(&iInfite);
    if (yTemp->size[0] > maxN) {
        i2 = yTemp->size[0];
        yTemp->size[0] = (int)maxN;
        emxEnsureCapacity((emxArray__common *)yTemp, i2, (int)sizeof(double));
    }

    i2 = iPk->size[0];
    iPk->size[0] = yTemp->size[0];
    emxEnsureCapacity((emxArray__common *)iPk, i2, (int)sizeof(double));
    cdiff = yTemp->size[0];
    for (i2 = 0; i2 < cdiff; i2++) {
        iPk->data[i2] = b_iPk->data[(int)yTemp->data[i2] - 1];
    }

    emxFree_real_T(&b_iPk);
    i2 = Ypk->size[0] * Ypk->size[1];
    Ypk->size[0] = yTemp->size[0];
    Ypk->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)Ypk, i2, (int)sizeof(double));
    cdiff = yTemp->size[0];
    for (i2 = 0; i2 < cdiff; i2++) {
        Ypk->data[i2] = y->data[(int)iPk->data[i2] - 1];
    }

    emxFree_real_T(&y);
    i2 = Xpk->size[0] * Xpk->size[1];
    Xpk->size[0] = yTemp->size[0];
    Xpk->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)Xpk, i2, (int)sizeof(double));
    cdiff = yTemp->size[0];
    emxFree_real_T(&yTemp);
    for (i2 = 0; i2 < cdiff; i2++) {
        Xpk->data[i2] = x->data[(int)iPk->data[i2] - 1];
    }

    emxFree_real_T(&x);
    emxFree_real_T(&iPk);
}

/*
 * File trailer for findpeaks.c
 *
 * [EOF]
 */

/* Function Definitions */

/*
 * Arguments    : const double x[9]
 *                double y[8]
 * Return Type  : void
 */
void b_diff(const double x[9], double y[8])
{
    int ixLead;
    int iyLead;
    double work;
    int m;
    double tmp2;
    ixLead = 1;
    iyLead = 0;
    work = x[0];
    for (m = 0; m < 8; m++) {
        tmp2 = work;
        work = x[ixLead];
        tmp2 = x[ixLead] - tmp2;
        ixLead++;
        y[iyLead] = tmp2;
        iyLead++;
    }
}

/*
 * Arguments    : const emxArray_real_T *x
 *                emxArray_real_T *y
 * Return Type  : void
 */
void diff(const emxArray_real_T *x, emxArray_real_T *y)
{
    int iyLead;
    int orderForDim;
    double work_data_idx_0;
    int m;
    double tmp1;
    double tmp2;
    if (x->size[0] == 0) {
        iyLead = y->size[0];
        y->size[0] = 0;
        emxEnsureCapacity((emxArray__common *)y, iyLead, (int)sizeof(double));
    } else {
        if (x->size[0] - 1 <= 1) {
            orderForDim = x->size[0] - 1;
        } else {
            orderForDim = 1;
        }

        if (orderForDim < 1) {
            iyLead = y->size[0];
            y->size[0] = 0;
            emxEnsureCapacity((emxArray__common *)y, iyLead, (int)sizeof(double));
        } else {
            orderForDim = x->size[0] - 1;
            iyLead = y->size[0];
            y->size[0] = orderForDim;
            emxEnsureCapacity((emxArray__common *)y, iyLead, (int)sizeof(double));
            if (!(y->size[0] == 0)) {
                orderForDim = 1;
                iyLead = 0;
                work_data_idx_0 = x->data[0];
                for (m = 2; m <= x->size[0]; m++) {
                    tmp1 = x->data[orderForDim];
                    tmp2 = work_data_idx_0;
                    work_data_idx_0 = tmp1;
                    tmp1 -= tmp2;
                    orderForDim++;
                    y->data[iyLead] = tmp1;
                    iyLead++;
                }
            }
        }
    }
}

/*
 * File trailer for diff.c
 *
 * [EOF]
 */
/* Function Definitions */

/*
 * Arguments    : const emxArray_boolean_T *x
 * Return Type  : boolean_T
 */
boolean_T all(const emxArray_boolean_T *x)
{
    boolean_T y;
    int ix;
    boolean_T exitg1;
    y = true;
    ix = 1;
    exitg1 = false;
    while ((!exitg1) && (ix <= x->size[1])) {
        if (!x->data[ix - 1]) {
            y = false;
            exitg1 = true;
        } else {
            ix++;
        }
    }

    return y;
}

/*
 * File trailer for all.c
 *
 * [EOF]
 */

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : int
 */
int b_cfprintf(void)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[27] = { 'N', 'o', 't', ' ', 'a', 'b', 'l', 'e', ' ',
                                   't', 'o', ' ', 'd', 'e', 't', 'e', 'c', 't', ' ', 'q', ' ', 'w', 'a', 'v',
                                   'e', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    b_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, "%s",cfmt);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : void
 * Return Type  : int
 */
int c_cfprintf(void)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[27] = { 'N', 'o', 't', ' ', 'a', 'b', 'l', 'e', ' ',
                                   't', 'o', ' ', 'd', 'e', 't', 'e', 'c', 't', ' ', 's', ' ', 'w', 'a', 'v',
                                   'e', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    c_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar,"%s",cfmt);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : void
 * Return Type  : int
 */
int cfprintf(void)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[27] = { 'N', 'o', 't', ' ', 'a', 'b', 'l', 'e', ' ',
                                   't', 'o', ' ', 'd', 'e', 't', 'e', 'c', 't', ' ', 'p', ' ', 'w', 'a', 'v',
                                   'e', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, "%s",cfmt);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : void
 * Return Type  : int
 */
int d_cfprintf(void)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[27] = { 'N', 'o', 't', ' ', 'a', 'b', 'l', 'e', ' ',
                                   't', 'o', ' ', 'd', 'e', 't', 'e', 'c', 't', ' ', 't', ' ', 'w', 'a', 'v',
                                   'e', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    d_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, "%s",cfmt);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : double varargin_1
 * Return Type  : int
 */
int e_cfprintf(double varargin_1)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[28] = { 'P', 'R', ' ', 'i', 'n', 't', 'e', 'r', 'v',
                                   'a', 'l', ' ', 'i', 's', ' ', ':', '-', '>', ' ', '%', 'f', ' ', 'm', 's',
                                   'e', 'c', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    e_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, cfmt, varargin_1);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : void
 * Return Type  : int
 */
int f_cfprintf(void)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[79] = { 'I', 'n', 's', 'u', 'f', 'f', 'i', 'c', 'i',
                                   'e', 'n', 't', ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ', 'o', 'f', ' ', 'R',
                                   ' ', 'a', 'n', 'd', ' ', 'p', ' ', 'p', 'e', 'a', 'k', 's', ' ', 'f', 'o',
                                   'u', 'n', 'd', '!', '!', ',', 't', 'h', 'e', 'r', 'e', ' ', 'n', 'u', 'm',
                                   'b', 'e', 'r', ' ', 's', 'h', 'o', 'u', 'l', 'd', ' ', 'b', 'e', ' ', 'a',
                                   't', 'l', 'e', 'a', 's', 't', ' ', '3', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    f_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, "%s",cfmt);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : double varargin_1
 * Return Type  : int
 */
int g_cfprintf(double varargin_1)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[29] = { 'Q', 'R', 'S', ' ', 'i', 'n', 't', 'e', 'r',
                                   'v', 'a', 'l', ' ', 'i', 's', ' ', ':', '-', '>', ' ', '%', 'f', ' ', 'm',
                                   's', 'e', 'c', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    g_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, cfmt, varargin_1);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : void
 * Return Type  : int
 */
int h_cfprintf(void)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[80] = { 'I', 'n', 's', 'u', 'f', 'f', 'i', 'c', 'i',
                                   'e', 'n', 't', ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ', 'o', 'f', ' ', 'S',
                                   ' ', ' ', 'a', 'n', 'd', ' ', 'q', ' ', 'p', 'e', 'a', 'k', 's', ' ', 'f',
                                   'o', 'u', 'n', 'd', '!', '!', ',', 't', 'h', 'e', 'r', 'e', ' ', 'n', 'u',
                                   'm', 'b', 'e', 'r', ' ', 's', 'h', 'o', 'u', 'l', 'd', ' ', 'b', 'e', ' ',
                                   'a', 't', 'l', 'e', 'a', 's', 't', ' ', '3', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    h_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, "%s",cfmt);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : double varargin_1
 * Return Type  : int
 */
int i_cfprintf(double varargin_1)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[28] = { 'Q', 'T', ' ', 'i', 'n', 't', 'e', 'r', 'v',
                                   'a', 'l', ' ', 'i', 's', ' ', ':', '-', '>', ' ', '%', 'f', ' ', 'm', 's',
                                   'e', 'c', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    i_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, cfmt, varargin_1);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : void
 * Return Type  : int
 */
int j_cfprintf(void)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[79] = { 'I', 'n', 's', 'u', 'f', 'f', 'i', 'e', 'c',
                                   'i', 'e', 'n', 't', ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ', 'o', 'f', ' ',
                                   'Q', ' ', 'a', 'n', 'd', ' ', 'T', ' ', 'p', 'e', 'a', 'k', 's', ' ', 'f',
                                   'o', 'u', 'n', 'd', ' ', ',', 't', 'h', 'e', 'r', 'e', ' ', 'n', 'u', 'm',
                                   'b', 'e', 'r', ' ', 's', 'h', 'o', 'u', 'l', 'd', ' ', 'b', 'e', ' ', 'a',
                                   't', 'l', 'e', 'a', 's', 't', ' ', '3', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    j_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, "%s",cfmt);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : double varargin_1
 * Return Type  : int
 */
int k_cfprintf(double varargin_1)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[21] = { 'H', 'e', 'a', 'r', 't', 'r', 'a', 't', 'e',
                                   ':', '-', '>', ' ', '%', 'f', ' ', 'b', 'p', 'm', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    k_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, cfmt, varargin_1);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : double varargin_1
 * Return Type  : int
 */
int l_cfprintf(double varargin_1)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[28] = { 'R', 'R', ' ', 'i', 'n', 't', 'e', 'r', 'v',
                                   'a', 'l', ' ', 'i', 's', ' ', ':', '-', '>', ' ', '%', 'f', ' ', 'm', 's',
                                   'e', 'c', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    l_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, cfmt, varargin_1);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : void
 * Return Type  : int
 */
int m_cfprintf(void)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[73] = { 'I', 'n', 's', 'u', 'f', 'f', 'i', 'e', 'c',
                                   'i', 'e', 'n', 't', ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ', 'o', 'f', ' ',
                                   'R', ' ', 'p', 'e', 'a', 'k', 's', ' ', 'f', 'o', 'u', 'n', 'd', ' ', ',',
                                   't', 'h', 'e', 'r', 'e', ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ', 's', 'h',
                                   'o', 'u', 'l', 'd', ' ', 'b', 'e', ' ', 'a', 't', 'l', 'e', 'a', 's', 't',
                                   ' ', '3', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    m_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar,"%s", cfmt);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : double varargin_1
 * Return Type  : int
 */
int n_cfprintf(double varargin_1)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[25] = { 'Q', 'T', 'c', ' ', 'i', 'n', 't', 'e', 'r',
                                   'v', 'a', 'l', ' ', 'i', 's', ' ', ':', '-', '>', ' ', '%', 'f', ' ', '\x0a',
                                   '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    n_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, cfmt, varargin_1);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : void
 * Return Type  : int
 */
int o_cfprintf(void)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[82] = { 'I', 'n', 's', 'u', 'f', 'f', 'i', 'e', 'c',
                                   'i', 'e', 'n', 't', ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ', 'o', 'f', ' ',
                                   'R', ' ', ',', 'Q', ' ', 'a', 'n', 'd', ' ', 'T', ' ', 'p', 'e', 'a', 'k',
                                   's', ' ', 'f', 'o', 'u', 'n', 'd', ' ', ',', 't', 'h', 'e', 'r', 'e', ' ',
                                   'n', 'u', 'm', 'b', 'e', 'r', ' ', 's', 'h', 'o', 'u', 'l', 'd', ' ', 'b',
                                   'e', ' ', 'a', 't', 'l', 'e', 'a', 's', 't', ' ', '3', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    o_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, "%s",cfmt);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : double varargin_1
 * Return Type  : int
 */
int p_cfprintf(double varargin_1)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[28] = { 'S', 'T', ' ', 'i', 'n', 't', 'e', 'r', 'v',
                                   'a', 'l', ' ', 'i', 's', ' ', ':', '-', '>', ' ', '%', 'f', ' ', 'm', 's',
                                   'e', 'c', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    p_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, cfmt, varargin_1);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : void
 * Return Type  : int
 */
int q_cfprintf(void)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[79] = { 'I', 'n', 's', 'u', 'f', 'f', 'i', 'e', 'c',
                                   'i', 'e', 'n', 't', ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ', 'o', 'f', ' ',
                                   'S', ' ', 'a', 'n', 'd', ' ', 'T', ' ', 'p', 'e', 'a', 'k', 's', ' ', 'f',
                                   'o', 'u', 'n', 'd', ' ', ',', 't', 'h', 'e', 'r', 'e', ' ', 'n', 'u', 'm',
                                   'b', 'e', 'r', ' ', 's', 'h', 'o', 'u', 'l', 'd', ' ', 'b', 'e', ' ', 'a',
                                   't', 'l', 'e', 'a', 's', 't', ' ', '3', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    q_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, "%s",cfmt);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : double varargin_1
 * Return Type  : int
 */
int r_cfprintf(double varargin_1)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[33] = { 'A', 'm', 'p', 'l', 'i', 't', 'u', 'd', 'e',
                                   ' ', 'p', ' ', 'w', 'a', 'v', 'e', ' ', 'i', 's', ':', '-', '>', ' ', '%',
                                   'f', ' ', 'm', 'V', 'o', 'l', 't', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    r_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, cfmt, varargin_1);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : double varargin_1
 * Return Type  : int
 */
int s_cfprintf(double varargin_1)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[33] = { 'A', 'm', 'p', 'l', 'i', 't', 'u', 'd', 'e',
                                   ' ', 'q', ' ', 'w', 'a', 'v', 'e', ' ', 'i', 's', ':', '-', '>', ' ', '%',
                                   'f', ' ', 'm', 'V', 'o', 'l', 't', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    s_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, cfmt, varargin_1);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : double varargin_1
 * Return Type  : int
 */
int t_cfprintf(double varargin_1)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[33] = { 'A', 'm', 'p', 'l', 'i', 't', 'u', 'd', 'e',
                                   ' ', 'r', ' ', 'w', 'a', 'v', 'e', ' ', 'i', 's', ':', '-', '>', ' ', '%',
                                   'f', ' ', 'm', 'V', 'o', 'l', 't', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    t_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, cfmt, varargin_1);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : double varargin_1
 * Return Type  : int
 */
int u_cfprintf(double varargin_1)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[33] = { 'A', 'm', 'p', 'l', 'i', 't', 'u', 'd', 'e',
                                   ' ', 's', ' ', 'w', 'a', 'v', 'e', ' ', 'i', 's', ':', '-', '>', ' ', '%',
                                   'f', ' ', 'm', 'V', 'o', 'l', 't', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    u_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, cfmt, varargin_1);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : double varargin_1
 * Return Type  : int
 */
int v_cfprintf(double varargin_1)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[33] = { 'A', 'm', 'p', 'l', 'i', 't', 'u', 'd', 'e',
                                   ' ', 't', ' ', 'w', 'a', 'v', 'e', ' ', 'i', 's', ':', '-', '>', ' ', '%',
                                   'f', ' ', 'm', 'V', 'o', 'l', 't', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    v_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, cfmt, varargin_1);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * Arguments    : void
 * Return Type  : int
 */
int w_cfprintf(void)
{
    int nbytesint;
    FILE * b_NULL;
    FILE * filestar;
    boolean_T autoflush;
    static const char cfmt[79] = { 'I', 'n', 's', 'u', 'f', 'f', 'i', 'e', 'c',
                                   'i', 'e', 'n', 't', ' ', 'n', 'u', 'm', 'b', 'e', 'r', ' ', 'o', 'f', ' ',
                                   'R', ' ', 'a', 'n', 'd', ' ', 'S', ' ', 'p', 'e', 'a', 'k', 's', ' ', 'f',
                                   'o', 'u', 'n', 'd', ' ', ',', 't', 'h', 'e', 'r', 'e', ' ', 'n', 'u', 'm',
                                   'b', 'e', 'r', ' ', 's', 'h', 'o', 'u', 'l', 'd', ' ', 'b', 'e', ' ', 'a',
                                   't', 'l', 'e', 'a', 's', 't', ' ', '3', '\x0a', '\x00' };

    b_NULL = NULL;
    nbytesint = 0;
    w_fileManager(&filestar, &autoflush);
    if (filestar == b_NULL) {
    } else {
        nbytesint = fprintf(filestar, "%s",cfmt);
        fflush(filestar);
    }

    return nbytesint;
}

/*
 * File trailer for fprintf.c
 *
 * [EOF]
 */

/* Function Definitions */

/*
 * Arguments    : const double b[7]
 *                const double x[6286]
 *                const double zi[6]
 *                double y[6286]
 * Return Type  : void
 */
void b_filter(const double b[7], const double x[6286], const double zi[6],
              double y[6286])
{
    double dbuffer[7];
    int k;
    int j;
    for (k = 0; k < 6; k++) {
        dbuffer[k + 1] = zi[k];
    }

    for (j = 0; j < 6286; j++) {
        for (k = 0; k < 6; k++) {
            dbuffer[k] = dbuffer[k + 1];
        }

        dbuffer[6] = 0.0;
        for (k = 0; k < 7; k++) {
            dbuffer[k] += x[j] * b[k];
        }

        y[j] = dbuffer[0];
    }
}

/*
 * Arguments    : const double b[7]
 *                const double a[7]
 *                const double x[6286]
 *                const double zi[6]
 *                double y[6286]
 * Return Type  : void
 */
void filter(const double b[7], const double a[7], const double x[6286], const
double zi[6], double y[6286])
{
    double dbuffer[7];
    int k;
    int j;
    for (k = 0; k < 6; k++) {
        dbuffer[k + 1] = zi[k];
    }

    for (j = 0; j < 6286; j++) {
        for (k = 0; k < 6; k++) {
            dbuffer[k] = dbuffer[k + 1];
        }

        dbuffer[6] = 0.0;
        for (k = 0; k < 7; k++) {
            dbuffer[k] += x[j] * b[k];
        }

        for (k = 0; k < 6; k++) {
            dbuffer[k + 1] -= dbuffer[0] * a[k + 1];
        }

        y[j] = dbuffer[0];
    }
}

/*
 * File trailer for filter.c
 *
 * [EOF]
 */

/* Function Definitions */

/*
 * Arguments    : emxArray_real_T *x
 *                int dim
 *                emxArray_int32_T *idx
 * Return Type  : void
 */
static void b_sort(emxArray_real_T *x, int dim, emxArray_int32_T *idx)
{
    int i4;
    emxArray_real_T *vwork;
    int vstride;
    int x_idx_0;
    int j;
    emxArray_int32_T *iidx;
    if (dim <= 1) {
        i4 = x->size[0];
    } else {
        i4 = 1;
    }

    emxInit_real_T(&vwork, 1);
    vstride = vwork->size[0];
    vwork->size[0] = i4;
    emxEnsureCapacity((emxArray__common *)vwork, vstride, (int)sizeof(double));
    x_idx_0 = x->size[0];
    vstride = idx->size[0];
    idx->size[0] = x_idx_0;
    emxEnsureCapacity((emxArray__common *)idx, vstride, (int)sizeof(int));
    vstride = 1;
    x_idx_0 = 1;
    while (x_idx_0 <= dim - 1) {
        vstride *= x->size[0];
        x_idx_0 = 2;
    }

    j = 0;
    emxInit_int32_T(&iidx, 1);
    while (j + 1 <= vstride) {
        for (x_idx_0 = 0; x_idx_0 + 1 <= i4; x_idx_0++) {
            vwork->data[x_idx_0] = x->data[j + x_idx_0 * vstride];
        }

        sortIdx(vwork, iidx);
        for (x_idx_0 = 0; x_idx_0 + 1 <= i4; x_idx_0++) {
            x->data[j + x_idx_0 * vstride] = vwork->data[x_idx_0];
            idx->data[j + x_idx_0 * vstride] = iidx->data[x_idx_0];
        }

        j++;
    }

    emxFree_int32_T(&iidx);
    emxFree_real_T(&vwork);
}

/*
 * Arguments    : emxArray_real_T *x
 * Return Type  : void
 */
void c_sort(emxArray_real_T *x)
{
    int dim;
    int i5;
    emxArray_real_T *vwork;
    int j;
    int vstride;
    int k;
    emxArray_int32_T *b_vwork;
    dim = 2;
    if (x->size[0] != 1) {
        dim = 1;
    }

    if (dim <= 1) {
        i5 = x->size[0];
    } else {
        i5 = 1;
    }

    emxInit_real_T(&vwork, 1);
    j = vwork->size[0];
    vwork->size[0] = i5;
    emxEnsureCapacity((emxArray__common *)vwork, j, (int)sizeof(double));
    vstride = 1;
    k = 1;
    while (k <= dim - 1) {
        vstride *= x->size[0];
        k = 2;
    }

    j = 0;
    emxInit_int32_T(&b_vwork, 1);
    while (j + 1 <= vstride) {
        for (k = 0; k + 1 <= i5; k++) {
            vwork->data[k] = x->data[j + k * vstride];
        }

        b_sortIdx(vwork, b_vwork);
        for (k = 0; k + 1 <= i5; k++) {
            x->data[j + k * vstride] = vwork->data[k];
        }

        j++;
    }

    emxFree_int32_T(&b_vwork);
    emxFree_real_T(&vwork);
}

/*
 * Arguments    : emxArray_real_T *x
 *                emxArray_int32_T *idx
 * Return Type  : void
 */
void sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
    int dim;
    dim = 2;
    if (x->size[0] != 1) {
        dim = 1;
    }

    b_sort(x, dim, idx);
}

/*
 * File trailer for sort1.c
 *
 * [EOF]
 */


/* Function Definitions */

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void b_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void c_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void d_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void e_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void f_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void g_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void h_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void i_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void j_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void k_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void l_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void m_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void n_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void o_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void p_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void q_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void r_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void s_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void t_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void u_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void v_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * Arguments    : FILE * *f
 *                boolean_T *a
 * Return Type  : void
 */
void w_fileManager(FILE * *f, boolean_T *a)
{
    *f = stdout;
    *a = true;
}

/*
 * File trailer for fileManager.c
 *
 * [EOF]
 */


/* Function Declarations */
static double skip_to_last_equal_value(int *k, const emxArray_real_T *x);

/* Function Definitions */

/*
 * Arguments    : int *k
 *                const emxArray_real_T *x
 * Return Type  : double
 */
static double skip_to_last_equal_value(int *k, const emxArray_real_T *x)
{
    double xk;
    boolean_T exitg1;
    double absxk;
    int exponent;
    boolean_T p;
    xk = x->data[*k - 1];
    exitg1 = false;
    while ((!exitg1) && (*k < x->size[0])) {
        absxk = fabs(xk / 2.0);
        if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
            if (absxk <= 2.2250738585072014E-308) {
                absxk = 4.94065645841247E-324;
            } else {
                frexp(absxk, &exponent);
                absxk = ldexp(1.0, exponent - 53);
            }
        } else {
            absxk = rtNaN;
        }

        if ((fabs(xk - x->data[*k]) < absxk) || (rtIsInf(x->data[*k]) && rtIsInf(xk)))
        {
            p = true;
        } else {
            p = false;
        }

        if (p) {
            (*k)++;
        } else {
            exitg1 = true;
        }
    }

    return xk;
}

/*
 * Arguments    : const emxArray_real_T *a
 *                const emxArray_real_T *b
 *                emxArray_real_T *c
 *                emxArray_int32_T *ia
 *                emxArray_int32_T *ib
 * Return Type  : void
 */
void do_vectors(const emxArray_real_T *a, const emxArray_real_T *b,
                emxArray_real_T *c, emxArray_int32_T *ia, emxArray_int32_T *ib)
{
    int na;
    int nb;
    int ncmax;
    int ibfirst;
    int nc;
    int nia;
    int nib;
    int iafirst;
    int ialast;
    int iblast;
    int b_ialast;
    double ak;
    int b_iblast;
    double bk;
    double absxk;
    emxArray_int32_T *b_ia;
    int exponent;
    emxArray_int32_T *b_ib;
    boolean_T p;
    emxArray_real_T *b_c;
    na = a->size[0];
    nb = b->size[0];
    ncmax = a->size[0] + b->size[0];
    ibfirst = c->size[0];
    c->size[0] = ncmax;
    emxEnsureCapacity((emxArray__common *)c, ibfirst, (int)sizeof(double));
    ibfirst = ia->size[0];
    ia->size[0] = a->size[0];
    emxEnsureCapacity((emxArray__common *)ia, ibfirst, (int)sizeof(int));
    ibfirst = ib->size[0];
    ib->size[0] = b->size[0];
    emxEnsureCapacity((emxArray__common *)ib, ibfirst, (int)sizeof(int));
    nc = -1;
    nia = -1;
    nib = 0;
    iafirst = 1;
    ialast = 1;
    ibfirst = 0;
    iblast = 1;
    while ((ialast <= na) && (iblast <= nb)) {
        b_ialast = ialast;
        ak = skip_to_last_equal_value(&b_ialast, a);
        ialast = b_ialast;
        b_iblast = iblast;
        bk = skip_to_last_equal_value(&b_iblast, b);
        iblast = b_iblast;
        absxk = fabs(bk / 2.0);
        if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
            if (absxk <= 2.2250738585072014E-308) {
                absxk = 4.94065645841247E-324;
            } else {
                frexp(absxk, &exponent);
                absxk = ldexp(1.0, exponent - 53);
            }
        } else {
            absxk = rtNaN;
        }

        if ((fabs(bk - ak) < absxk) || (rtIsInf(ak) && rtIsInf(bk))) {
            p = true;
        } else {
            p = false;
        }

        if (p) {
            nc++;
            c->data[nc] = ak;
            nia++;
            ia->data[nia] = iafirst;
            ialast = b_ialast + 1;
            iafirst = b_ialast + 1;
            iblast = b_iblast + 1;
            ibfirst = b_iblast;
        } else {
            if ((ak < bk) || rtIsNaN(bk)) {
                p = true;
            } else {
                p = false;
            }

            if (p) {
                nc++;
                nia++;
                c->data[nc] = ak;
                ia->data[nia] = iafirst;
                ialast = b_ialast + 1;
                iafirst = b_ialast + 1;
            } else {
                nc++;
                nib++;
                c->data[nc] = bk;
                ib->data[nib - 1] = ibfirst + 1;
                iblast = b_iblast + 1;
                ibfirst = b_iblast;
            }
        }
    }

    while (ialast <= na) {
        iafirst = ialast;
        ak = skip_to_last_equal_value(&iafirst, a);
        nc++;
        nia++;
        c->data[nc] = ak;
        ia->data[nia] = ialast;
        ialast = iafirst + 1;
    }

    while (iblast <= nb) {
        iafirst = iblast;
        bk = skip_to_last_equal_value(&iafirst, b);
        nc++;
        nib++;
        c->data[nc] = bk;
        ib->data[nib - 1] = iblast;
        iblast = iafirst + 1;
    }

    if (a->size[0] > 0) {
        if (1 > nia + 1) {
            iafirst = -1;
        } else {
            iafirst = nia;
        }

        emxInit_int32_T(&b_ia, 1);
        ibfirst = b_ia->size[0];
        b_ia->size[0] = iafirst + 1;
        emxEnsureCapacity((emxArray__common *)b_ia, ibfirst, (int)sizeof(int));
        for (ibfirst = 0; ibfirst <= iafirst; ibfirst++) {
            b_ia->data[ibfirst] = ia->data[ibfirst];
        }

        ibfirst = ia->size[0];
        ia->size[0] = b_ia->size[0];
        emxEnsureCapacity((emxArray__common *)ia, ibfirst, (int)sizeof(int));
        iafirst = b_ia->size[0];
        for (ibfirst = 0; ibfirst < iafirst; ibfirst++) {
            ia->data[ibfirst] = b_ia->data[ibfirst];
        }

        emxFree_int32_T(&b_ia);
    }

    if (b->size[0] > 0) {
        if (1 > nib) {
            iafirst = 0;
        } else {
            iafirst = nib;
        }

        emxInit_int32_T(&b_ib, 1);
        ibfirst = b_ib->size[0];
        b_ib->size[0] = iafirst;
        emxEnsureCapacity((emxArray__common *)b_ib, ibfirst, (int)sizeof(int));
        for (ibfirst = 0; ibfirst < iafirst; ibfirst++) {
            b_ib->data[ibfirst] = ib->data[ibfirst];
        }

        ibfirst = ib->size[0];
        ib->size[0] = b_ib->size[0];
        emxEnsureCapacity((emxArray__common *)ib, ibfirst, (int)sizeof(int));
        iafirst = b_ib->size[0];
        for (ibfirst = 0; ibfirst < iafirst; ibfirst++) {
            ib->data[ibfirst] = b_ib->data[ibfirst];
        }

        emxFree_int32_T(&b_ib);
    }

    if (ncmax > 0) {
        if (1 > nc + 1) {
            iafirst = -1;
        } else {
            iafirst = nc;
        }

        emxInit_real_T(&b_c, 1);
        ibfirst = b_c->size[0];
        b_c->size[0] = iafirst + 1;
        emxEnsureCapacity((emxArray__common *)b_c, ibfirst, (int)sizeof(double));
        for (ibfirst = 0; ibfirst <= iafirst; ibfirst++) {
            b_c->data[ibfirst] = c->data[ibfirst];
        }

        ibfirst = c->size[0];
        c->size[0] = b_c->size[0];
        emxEnsureCapacity((emxArray__common *)c, ibfirst, (int)sizeof(double));
        iafirst = b_c->size[0];
        for (ibfirst = 0; ibfirst < iafirst; ibfirst++) {
            c->data[ibfirst] = b_c->data[ibfirst];
        }

        emxFree_real_T(&b_c);
    }
}

/*
 * File trailer for eml_setop.c
 *
 * [EOF]
 */
#define NumBitsPerChar                 8U

/* Function: rtGetNaN ==================================================
 * Abstract:
 * Initialize rtNaN needed by the generated code.
 * NaN is initialized as non-signaling. Assumes IEEE.
 */
real_T rtGetNaN(void)
{
    size_t bitsPerReal = sizeof(real_T) * (NumBitsPerChar);
    real_T nan = 0.0;
    if (bitsPerReal == 32U) {
        nan = rtGetNaNF();
    } else {
        uint16_T one = 1U;
        enum {
            LittleEndian,
            BigEndian
        } machByteOrder = (*((uint8_T *) &one) == 1U) ? LittleEndian : BigEndian;
        switch (machByteOrder) {
            case LittleEndian:
            {
                union {
                    LittleEndianIEEEDouble bitVal;
                    real_T fltVal;
                } tmpVal;

                tmpVal.bitVal.words.wordH = 0xFFF80000U;
                tmpVal.bitVal.words.wordL = 0x00000000U;
                nan = tmpVal.fltVal;
                break;
            }

            case BigEndian:
            {
                union {
                    BigEndianIEEEDouble bitVal;
                    real_T fltVal;
                } tmpVal;

                tmpVal.bitVal.words.wordH = 0x7FFFFFFFU;
                tmpVal.bitVal.words.wordL = 0xFFFFFFFFU;
                nan = tmpVal.fltVal;
                break;
            }
        }
    }

    return nan;
}

/* Function: rtGetNaNF ==================================================
 * Abstract:
 * Initialize rtNaNF needed by the generated code.
 * NaN is initialized as non-signaling. Assumes IEEE.
 */
real32_T rtGetNaNF(void)
{
    IEEESingle nanF = { { 0 } };

    uint16_T one = 1U;
    enum {
        LittleEndian,
        BigEndian
    } machByteOrder = (*((uint8_T *) &one) == 1U) ? LittleEndian : BigEndian;
    switch (machByteOrder) {
        case LittleEndian:
        {
            nanF.wordL.wordLuint = 0xFFC00000U;
            break;
        }

        case BigEndian:
        {
            nanF.wordL.wordLuint = 0x7FFFFFFFU;
            break;
        }
    }

    return nanF.wordL.wordLreal;
}

/*
 * File trailer for rtGetNaN.c
 *
 * [EOF]
 */


/* Function: rtGetInf ==================================================
 * Abstract:
 * Initialize rtInf needed by the generated code.
 * Inf is initialized as non-signaling. Assumes IEEE.
 */
real_T rtGetInf(void)
{
    size_t bitsPerReal = sizeof(real_T) * (NumBitsPerChar);
    real_T inf = 0.0;
    if (bitsPerReal == 32U) {
        inf = rtGetInfF();
    } else {
        uint16_T one = 1U;
        enum {
            LittleEndian,
            BigEndian
        } machByteOrder = (*((uint8_T *) &one) == 1U) ? LittleEndian : BigEndian;
        switch (machByteOrder) {
            case LittleEndian:
            {
                union {
                    LittleEndianIEEEDouble bitVal;
                    real_T fltVal;
                } tmpVal;

                tmpVal.bitVal.words.wordH = 0x7FF00000U;
                tmpVal.bitVal.words.wordL = 0x00000000U;
                inf = tmpVal.fltVal;
                break;
            }

            case BigEndian:
            {
                union {
                    BigEndianIEEEDouble bitVal;
                    real_T fltVal;
                } tmpVal;

                tmpVal.bitVal.words.wordH = 0x7FF00000U;
                tmpVal.bitVal.words.wordL = 0x00000000U;
                inf = tmpVal.fltVal;
                break;
            }
        }
    }

    return inf;
}

/* Function: rtGetInfF ==================================================
 * Abstract:
 * Initialize rtInfF needed by the generated code.
 * Inf is initialized as non-signaling. Assumes IEEE.
 */
real32_T rtGetInfF(void)
{
    IEEESingle infF;
    infF.wordL.wordLuint = 0x7F800000U;
    return infF.wordL.wordLreal;
}

/* Function: rtGetMinusInf ==================================================
 * Abstract:
 * Initialize rtMinusInf needed by the generated code.
 * Inf is initialized as non-signaling. Assumes IEEE.
 */
real_T rtGetMinusInf(void)
{
    size_t bitsPerReal = sizeof(real_T) * (NumBitsPerChar);
    real_T minf = 0.0;
    if (bitsPerReal == 32U) {
        minf = rtGetMinusInfF();
    } else {
        uint16_T one = 1U;
        enum {
            LittleEndian,
            BigEndian
        } machByteOrder = (*((uint8_T *) &one) == 1U) ? LittleEndian : BigEndian;
        switch (machByteOrder) {
            case LittleEndian:
            {
                union {
                    LittleEndianIEEEDouble bitVal;
                    real_T fltVal;
                } tmpVal;

                tmpVal.bitVal.words.wordH = 0xFFF00000U;
                tmpVal.bitVal.words.wordL = 0x00000000U;
                minf = tmpVal.fltVal;
                break;
            }

            case BigEndian:
            {
                union {
                    BigEndianIEEEDouble bitVal;
                    real_T fltVal;
                } tmpVal;

                tmpVal.bitVal.words.wordH = 0xFFF00000U;
                tmpVal.bitVal.words.wordL = 0x00000000U;
                minf = tmpVal.fltVal;
                break;
            }
        }
    }

    return minf;
}

/* Function: rtGetMinusInfF ==================================================
 * Abstract:
 * Initialize rtMinusInfF needed by the generated code.
 * Inf is initialized as non-signaling. Assumes IEEE.
 */
real32_T rtGetMinusInfF(void)
{
    IEEESingle minfF;
    minfF.wordL.wordLuint = 0xFF800000U;
    return minfF.wordL.wordLreal;
}

/*
 * File trailer for rtGetInf.c
 *
 * [EOF]
 */

/* Function Declarations */
static void b_merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int
np, int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void b_merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
                          int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
                        int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork);

/* Function Definitions */

/*
 * Arguments    : emxArray_int32_T *idx
 *                emxArray_real_T *x
 *                int offset
 *                int np
 *                int nq
 *                emxArray_int32_T *iwork
 *                emxArray_real_T *xwork
 * Return Type  : void
 */
static void b_merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int
np, int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
    int n;
    int qend;
    int p;
    int iout;
    int exitg1;
    if (nq == 0) {
    } else {
        n = np + nq;
        for (qend = 0; qend + 1 <= n; qend++) {
            iwork->data[qend] = idx->data[offset + qend];
            xwork->data[qend] = x->data[offset + qend];
        }

        p = 0;
        n = np;
        qend = np + nq;
        iout = offset - 1;
        do {
            exitg1 = 0;
            iout++;
            if (xwork->data[p] <= xwork->data[n]) {
                idx->data[iout] = iwork->data[p];
                x->data[iout] = xwork->data[p];
                if (p + 1 < np) {
                    p++;
                } else {
                    exitg1 = 1;
                }
            } else {
                idx->data[iout] = iwork->data[n];
                x->data[iout] = xwork->data[n];
                if (n + 1 < qend) {
                    n++;
                } else {
                    n = (iout - p) + 1;
                    while (p + 1 <= np) {
                        idx->data[n + p] = iwork->data[p];
                        x->data[n + p] = xwork->data[p];
                        p++;
                    }

                    exitg1 = 1;
                }
            }
        } while (exitg1 == 0);
    }
}

/*
 * Arguments    : emxArray_int32_T *idx
 *                emxArray_real_T *x
 *                int offset
 *                int n
 *                int preSortLevel
 *                emxArray_int32_T *iwork
 *                emxArray_real_T *xwork
 * Return Type  : void
 */
static void b_merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
                          int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
    int nPairs;
    int bLen;
    int tailOffset;
    int nTail;
    nPairs = n >> preSortLevel;
    bLen = 1 << preSortLevel;
    while (nPairs > 1) {
        if ((nPairs & 1) != 0) {
            nPairs--;
            tailOffset = bLen * nPairs;
            nTail = n - tailOffset;
            if (nTail > bLen) {
                b_merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
            }
        }

        tailOffset = bLen << 1;
        nPairs >>= 1;
        for (nTail = 1; nTail <= nPairs; nTail++) {
            b_merge(idx, x, offset + (nTail - 1) * tailOffset, bLen, bLen, iwork,
                    xwork);
        }

        bLen = tailOffset;
    }

    if (n > bLen) {
        b_merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
    }
}

/*
 * Arguments    : emxArray_int32_T *idx
 *                emxArray_real_T *x
 *                int offset
 *                int np
 *                int nq
 *                emxArray_int32_T *iwork
 *                emxArray_real_T *xwork
 * Return Type  : void
 */
static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
    int n;
    int qend;
    int p;
    int iout;
    int exitg1;
    if (nq == 0) {
    } else {
        n = np + nq;
        for (qend = 0; qend + 1 <= n; qend++) {
            iwork->data[qend] = idx->data[offset + qend];
            xwork->data[qend] = x->data[offset + qend];
        }

        p = 0;
        n = np;
        qend = np + nq;
        iout = offset - 1;
        do {
            exitg1 = 0;
            iout++;
            if (xwork->data[p] >= xwork->data[n]) {
                idx->data[iout] = iwork->data[p];
                x->data[iout] = xwork->data[p];
                if (p + 1 < np) {
                    p++;
                } else {
                    exitg1 = 1;
                }
            } else {
                idx->data[iout] = iwork->data[n];
                x->data[iout] = xwork->data[n];
                if (n + 1 < qend) {
                    n++;
                } else {
                    n = (iout - p) + 1;
                    while (p + 1 <= np) {
                        idx->data[n + p] = iwork->data[p];
                        x->data[n + p] = xwork->data[p];
                        p++;
                    }

                    exitg1 = 1;
                }
            }
        } while (exitg1 == 0);
    }
}

/*
 * Arguments    : emxArray_int32_T *idx
 *                emxArray_real_T *x
 *                int offset
 *                int n
 *                int preSortLevel
 *                emxArray_int32_T *iwork
 *                emxArray_real_T *xwork
 * Return Type  : void
 */
static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
                        int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
    int nPairs;
    int bLen;
    int tailOffset;
    int nTail;
    nPairs = n >> preSortLevel;
    bLen = 1 << preSortLevel;
    while (nPairs > 1) {
        if ((nPairs & 1) != 0) {
            nPairs--;
            tailOffset = bLen * nPairs;
            nTail = n - tailOffset;
            if (nTail > bLen) {
                merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
            }
        }

        tailOffset = bLen << 1;
        nPairs >>= 1;
        for (nTail = 1; nTail <= nPairs; nTail++) {
            merge(idx, x, offset + (nTail - 1) * tailOffset, bLen, bLen, iwork, xwork);
        }

        bLen = tailOffset;
    }

    if (n > bLen) {
        merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
    }
}

/*
 * Arguments    : emxArray_real_T *x
 *                emxArray_int32_T *idx
 * Return Type  : void
 */
void b_sortIdx(emxArray_real_T *x, emxArray_int32_T *idx)
{
    emxArray_real_T *b_x;
    int ib;
    int wOffset;
    int m;
    int n;
    double x4[4];
    int idx4[4];
    emxArray_int32_T *iwork;
    emxArray_real_T *xwork;
    int nNaNs;
    int k;
    signed char perm[4];
    int nNonNaN;
    int i3;
    int i4;
    int nBlocks;
    int b_iwork[256];
    double b_xwork[256];
    int bLen2;
    int nPairs;
    int exitg1;
    emxInit_real_T(&b_x, 1);
    ib = x->size[0];
    wOffset = b_x->size[0];
    b_x->size[0] = x->size[0];
    emxEnsureCapacity((emxArray__common *)b_x, wOffset, (int)sizeof(double));
    m = x->size[0];
    for (wOffset = 0; wOffset < m; wOffset++) {
        b_x->data[wOffset] = x->data[wOffset];
    }

    wOffset = idx->size[0];
    idx->size[0] = ib;
    emxEnsureCapacity((emxArray__common *)idx, wOffset, (int)sizeof(int));
    for (wOffset = 0; wOffset < ib; wOffset++) {
        idx->data[wOffset] = 0;
    }

    n = x->size[0];
    for (m = 0; m < 4; m++) {
        x4[m] = 0.0;
        idx4[m] = 0;
    }

    emxInit_int32_T(&iwork, 1);
    wOffset = iwork->size[0];
    iwork->size[0] = ib;
    emxEnsureCapacity((emxArray__common *)iwork, wOffset, (int)sizeof(int));
    m = iwork->size[0];
    wOffset = iwork->size[0];
    iwork->size[0] = m;
    emxEnsureCapacity((emxArray__common *)iwork, wOffset, (int)sizeof(int));
    for (wOffset = 0; wOffset < m; wOffset++) {
        iwork->data[wOffset] = 0;
    }

    emxInit_real_T(&xwork, 1);
    m = x->size[0];
    wOffset = xwork->size[0];
    xwork->size[0] = m;
    emxEnsureCapacity((emxArray__common *)xwork, wOffset, (int)sizeof(double));
    m = xwork->size[0];
    wOffset = xwork->size[0];
    xwork->size[0] = m;
    emxEnsureCapacity((emxArray__common *)xwork, wOffset, (int)sizeof(double));
    for (wOffset = 0; wOffset < m; wOffset++) {
        xwork->data[wOffset] = 0.0;
    }

    nNaNs = 1;
    ib = 0;
    for (k = 0; k + 1 <= n; k++) {
        if (rtIsNaN(b_x->data[k])) {
            idx->data[n - nNaNs] = k + 1;
            xwork->data[n - nNaNs] = b_x->data[k];
            nNaNs++;
        } else {
            ib++;
            idx4[ib - 1] = k + 1;
            x4[ib - 1] = b_x->data[k];
            if (ib == 4) {
                ib = k - nNaNs;
                if (x4[0] <= x4[1]) {
                    m = 1;
                    wOffset = 2;
                } else {
                    m = 2;
                    wOffset = 1;
                }

                if (x4[2] <= x4[3]) {
                    i3 = 3;
                    i4 = 4;
                } else {
                    i3 = 4;
                    i4 = 3;
                }

                if (x4[m - 1] <= x4[i3 - 1]) {
                    if (x4[wOffset - 1] <= x4[i3 - 1]) {
                        perm[0] = (signed char)m;
                        perm[1] = (signed char)wOffset;
                        perm[2] = (signed char)i3;
                        perm[3] = (signed char)i4;
                    } else if (x4[wOffset - 1] <= x4[i4 - 1]) {
                        perm[0] = (signed char)m;
                        perm[1] = (signed char)i3;
                        perm[2] = (signed char)wOffset;
                        perm[3] = (signed char)i4;
                    } else {
                        perm[0] = (signed char)m;
                        perm[1] = (signed char)i3;
                        perm[2] = (signed char)i4;
                        perm[3] = (signed char)wOffset;
                    }
                } else if (x4[m - 1] <= x4[i4 - 1]) {
                    if (x4[wOffset - 1] <= x4[i4 - 1]) {
                        perm[0] = (signed char)i3;
                        perm[1] = (signed char)m;
                        perm[2] = (signed char)wOffset;
                        perm[3] = (signed char)i4;
                    } else {
                        perm[0] = (signed char)i3;
                        perm[1] = (signed char)m;
                        perm[2] = (signed char)i4;
                        perm[3] = (signed char)wOffset;
                    }
                } else {
                    perm[0] = (signed char)i3;
                    perm[1] = (signed char)i4;
                    perm[2] = (signed char)m;
                    perm[3] = (signed char)wOffset;
                }

                idx->data[ib - 2] = idx4[perm[0] - 1];
                idx->data[ib - 1] = idx4[perm[1] - 1];
                idx->data[ib] = idx4[perm[2] - 1];
                idx->data[ib + 1] = idx4[perm[3] - 1];
                b_x->data[ib - 2] = x4[perm[0] - 1];
                b_x->data[ib - 1] = x4[perm[1] - 1];
                b_x->data[ib] = x4[perm[2] - 1];
                b_x->data[ib + 1] = x4[perm[3] - 1];
                ib = 0;
            }
        }
    }

    wOffset = x->size[0] - nNaNs;
    if (ib > 0) {
        for (m = 0; m < 4; m++) {
            perm[m] = 0;
        }

        if (ib == 1) {
            perm[0] = 1;
        } else if (ib == 2) {
            if (x4[0] <= x4[1]) {
                perm[0] = 1;
                perm[1] = 2;
            } else {
                perm[0] = 2;
                perm[1] = 1;
            }
        } else if (x4[0] <= x4[1]) {
            if (x4[1] <= x4[2]) {
                perm[0] = 1;
                perm[1] = 2;
                perm[2] = 3;
            } else if (x4[0] <= x4[2]) {
                perm[0] = 1;
                perm[1] = 3;
                perm[2] = 2;
            } else {
                perm[0] = 3;
                perm[1] = 1;
                perm[2] = 2;
            }
        } else if (x4[0] <= x4[2]) {
            perm[0] = 2;
            perm[1] = 1;
            perm[2] = 3;
        } else if (x4[1] <= x4[2]) {
            perm[0] = 2;
            perm[1] = 3;
            perm[2] = 1;
        } else {
            perm[0] = 3;
            perm[1] = 2;
            perm[2] = 1;
        }

        for (k = 1; k <= ib; k++) {
            idx->data[(wOffset - ib) + k] = idx4[perm[k - 1] - 1];
            b_x->data[(wOffset - ib) + k] = x4[perm[k - 1] - 1];
        }
    }

    m = (nNaNs - 1) >> 1;
    for (k = 1; k <= m; k++) {
        ib = idx->data[wOffset + k];
        idx->data[wOffset + k] = idx->data[n - k];
        idx->data[n - k] = ib;
        b_x->data[wOffset + k] = xwork->data[n - k];
        b_x->data[n - k] = xwork->data[wOffset + k];
    }

    if (((nNaNs - 1) & 1) != 0) {
        b_x->data[(wOffset + m) + 1] = xwork->data[(wOffset + m) + 1];
    }

    nNonNaN = (x->size[0] - nNaNs) + 1;
    m = 2;
    if (nNonNaN > 1) {
        if (x->size[0] >= 256) {
            nBlocks = nNonNaN >> 8;
            if (nBlocks > 0) {
                for (i3 = 1; i3 <= nBlocks; i3++) {
                    i4 = ((i3 - 1) << 8) - 1;
                    for (nNaNs = 0; nNaNs < 6; nNaNs++) {
                        n = 1 << (nNaNs + 2);
                        bLen2 = n << 1;
                        nPairs = 256 >> (nNaNs + 3);
                        for (k = 1; k <= nPairs; k++) {
                            m = i4 + (k - 1) * bLen2;
                            for (ib = 1; ib <= bLen2; ib++) {
                                b_iwork[ib - 1] = idx->data[m + ib];
                                b_xwork[ib - 1] = b_x->data[m + ib];
                            }

                            wOffset = 0;
                            ib = n;
                            do {
                                exitg1 = 0;
                                m++;
                                if (b_xwork[wOffset] <= b_xwork[ib]) {
                                    idx->data[m] = b_iwork[wOffset];
                                    b_x->data[m] = b_xwork[wOffset];
                                    if (wOffset + 1 < n) {
                                        wOffset++;
                                    } else {
                                        exitg1 = 1;
                                    }
                                } else {
                                    idx->data[m] = b_iwork[ib];
                                    b_x->data[m] = b_xwork[ib];
                                    if (ib + 1 < bLen2) {
                                        ib++;
                                    } else {
                                        ib = m - wOffset;
                                        while (wOffset + 1 <= n) {
                                            idx->data[(ib + wOffset) + 1] = b_iwork[wOffset];
                                            b_x->data[(ib + wOffset) + 1] = b_xwork[wOffset];
                                            wOffset++;
                                        }

                                        exitg1 = 1;
                                    }
                                }
                            } while (exitg1 == 0);
                        }
                    }
                }

                m = nBlocks << 8;
                ib = nNonNaN - m;
                if (ib > 0) {
                    b_merge_block(idx, b_x, m, ib, 2, iwork, xwork);
                }

                m = 8;
            }
        }

        b_merge_block(idx, b_x, 0, nNonNaN, m, iwork, xwork);
    }

    emxFree_real_T(&xwork);
    emxFree_int32_T(&iwork);
    wOffset = x->size[0];
    x->size[0] = b_x->size[0];
    emxEnsureCapacity((emxArray__common *)x, wOffset, (int)sizeof(double));
    m = b_x->size[0];
    for (wOffset = 0; wOffset < m; wOffset++) {
        x->data[wOffset] = b_x->data[wOffset];
    }

    emxFree_real_T(&b_x);
}

/*
 * Arguments    : emxArray_real_T *x
 *                emxArray_int32_T *idx
 * Return Type  : void
 */
void sortIdx(emxArray_real_T *x, emxArray_int32_T *idx)
{
    emxArray_real_T *b_x;
    int ib;
    int wOffset;
    int m;
    int n;
    double x4[4];
    int idx4[4];
    emxArray_int32_T *iwork;
    emxArray_real_T *xwork;
    int nNaNs;
    int k;
    signed char perm[4];
    int nNonNaN;
    int i3;
    int i4;
    int nBlocks;
    int b_iwork[256];
    double b_xwork[256];
    int bLen;
    int bLen2;
    int nPairs;
    int exitg1;
    emxInit_real_T(&b_x, 1);
    ib = x->size[0];
    wOffset = b_x->size[0];
    b_x->size[0] = x->size[0];
    emxEnsureCapacity((emxArray__common *)b_x, wOffset, (int)sizeof(double));
    m = x->size[0];
    for (wOffset = 0; wOffset < m; wOffset++) {
        b_x->data[wOffset] = x->data[wOffset];
    }

    wOffset = idx->size[0];
    idx->size[0] = ib;
    emxEnsureCapacity((emxArray__common *)idx, wOffset, (int)sizeof(int));
    for (wOffset = 0; wOffset < ib; wOffset++) {
        idx->data[wOffset] = 0;
    }

    n = x->size[0];
    for (m = 0; m < 4; m++) {
        x4[m] = 0.0;
        idx4[m] = 0;
    }

    emxInit_int32_T(&iwork, 1);
    wOffset = iwork->size[0];
    iwork->size[0] = ib;
    emxEnsureCapacity((emxArray__common *)iwork, wOffset, (int)sizeof(int));
    m = iwork->size[0];
    wOffset = iwork->size[0];
    iwork->size[0] = m;
    emxEnsureCapacity((emxArray__common *)iwork, wOffset, (int)sizeof(int));
    for (wOffset = 0; wOffset < m; wOffset++) {
        iwork->data[wOffset] = 0;
    }

    emxInit_real_T(&xwork, 1);
    m = x->size[0];
    wOffset = xwork->size[0];
    xwork->size[0] = m;
    emxEnsureCapacity((emxArray__common *)xwork, wOffset, (int)sizeof(double));
    m = xwork->size[0];
    wOffset = xwork->size[0];
    xwork->size[0] = m;
    emxEnsureCapacity((emxArray__common *)xwork, wOffset, (int)sizeof(double));
    for (wOffset = 0; wOffset < m; wOffset++) {
        xwork->data[wOffset] = 0.0;
    }

    nNaNs = 0;
    ib = 0;
    for (k = 0; k + 1 <= n; k++) {
        if (rtIsNaN(b_x->data[k])) {
            idx->data[(n - nNaNs) - 1] = k + 1;
            xwork->data[(n - nNaNs) - 1] = b_x->data[k];
            nNaNs++;
        } else {
            ib++;
            idx4[ib - 1] = k + 1;
            x4[ib - 1] = b_x->data[k];
            if (ib == 4) {
                ib = k - nNaNs;
                if (x4[0] >= x4[1]) {
                    m = 1;
                    wOffset = 2;
                } else {
                    m = 2;
                    wOffset = 1;
                }

                if (x4[2] >= x4[3]) {
                    i3 = 3;
                    i4 = 4;
                } else {
                    i3 = 4;
                    i4 = 3;
                }

                if (x4[m - 1] >= x4[i3 - 1]) {
                    if (x4[wOffset - 1] >= x4[i3 - 1]) {
                        perm[0] = (signed char)m;
                        perm[1] = (signed char)wOffset;
                        perm[2] = (signed char)i3;
                        perm[3] = (signed char)i4;
                    } else if (x4[wOffset - 1] >= x4[i4 - 1]) {
                        perm[0] = (signed char)m;
                        perm[1] = (signed char)i3;
                        perm[2] = (signed char)wOffset;
                        perm[3] = (signed char)i4;
                    } else {
                        perm[0] = (signed char)m;
                        perm[1] = (signed char)i3;
                        perm[2] = (signed char)i4;
                        perm[3] = (signed char)wOffset;
                    }
                } else if (x4[m - 1] >= x4[i4 - 1]) {
                    if (x4[wOffset - 1] >= x4[i4 - 1]) {
                        perm[0] = (signed char)i3;
                        perm[1] = (signed char)m;
                        perm[2] = (signed char)wOffset;
                        perm[3] = (signed char)i4;
                    } else {
                        perm[0] = (signed char)i3;
                        perm[1] = (signed char)m;
                        perm[2] = (signed char)i4;
                        perm[3] = (signed char)wOffset;
                    }
                } else {
                    perm[0] = (signed char)i3;
                    perm[1] = (signed char)i4;
                    perm[2] = (signed char)m;
                    perm[3] = (signed char)wOffset;
                }

                idx->data[ib - 3] = idx4[perm[0] - 1];
                idx->data[ib - 2] = idx4[perm[1] - 1];
                idx->data[ib - 1] = idx4[perm[2] - 1];
                idx->data[ib] = idx4[perm[3] - 1];
                b_x->data[ib - 3] = x4[perm[0] - 1];
                b_x->data[ib - 2] = x4[perm[1] - 1];
                b_x->data[ib - 1] = x4[perm[2] - 1];
                b_x->data[ib] = x4[perm[3] - 1];
                ib = 0;
            }
        }
    }

    wOffset = (x->size[0] - nNaNs) - 1;
    if (ib > 0) {
        for (m = 0; m < 4; m++) {
            perm[m] = 0;
        }

        if (ib == 1) {
            perm[0] = 1;
        } else if (ib == 2) {
            if (x4[0] >= x4[1]) {
                perm[0] = 1;
                perm[1] = 2;
            } else {
                perm[0] = 2;
                perm[1] = 1;
            }
        } else if (x4[0] >= x4[1]) {
            if (x4[1] >= x4[2]) {
                perm[0] = 1;
                perm[1] = 2;
                perm[2] = 3;
            } else if (x4[0] >= x4[2]) {
                perm[0] = 1;
                perm[1] = 3;
                perm[2] = 2;
            } else {
                perm[0] = 3;
                perm[1] = 1;
                perm[2] = 2;
            }
        } else if (x4[0] >= x4[2]) {
            perm[0] = 2;
            perm[1] = 1;
            perm[2] = 3;
        } else if (x4[1] >= x4[2]) {
            perm[0] = 2;
            perm[1] = 3;
            perm[2] = 1;
        } else {
            perm[0] = 3;
            perm[1] = 2;
            perm[2] = 1;
        }

        for (k = 1; k <= ib; k++) {
            idx->data[(wOffset - ib) + k] = idx4[perm[k - 1] - 1];
            b_x->data[(wOffset - ib) + k] = x4[perm[k - 1] - 1];
        }
    }

    m = nNaNs >> 1;
    for (k = 1; k <= m; k++) {
        ib = idx->data[wOffset + k];
        idx->data[wOffset + k] = idx->data[n - k];
        idx->data[n - k] = ib;
        b_x->data[wOffset + k] = xwork->data[n - k];
        b_x->data[n - k] = xwork->data[wOffset + k];
    }

    if ((nNaNs & 1) != 0) {
        b_x->data[(wOffset + m) + 1] = xwork->data[(wOffset + m) + 1];
    }

    nNonNaN = x->size[0] - nNaNs;
    m = 2;
    if (nNonNaN > 1) {
        if (x->size[0] >= 256) {
            nBlocks = nNonNaN >> 8;
            if (nBlocks > 0) {
                for (i3 = 1; i3 <= nBlocks; i3++) {
                    i4 = ((i3 - 1) << 8) - 1;
                    for (n = 0; n < 6; n++) {
                        bLen = 1 << (n + 2);
                        bLen2 = bLen << 1;
                        nPairs = 256 >> (n + 3);
                        for (k = 1; k <= nPairs; k++) {
                            m = i4 + (k - 1) * bLen2;
                            for (ib = 1; ib <= bLen2; ib++) {
                                b_iwork[ib - 1] = idx->data[m + ib];
                                b_xwork[ib - 1] = b_x->data[m + ib];
                            }

                            wOffset = 0;
                            ib = bLen;
                            do {
                                exitg1 = 0;
                                m++;
                                if (b_xwork[wOffset] >= b_xwork[ib]) {
                                    idx->data[m] = b_iwork[wOffset];
                                    b_x->data[m] = b_xwork[wOffset];
                                    if (wOffset + 1 < bLen) {
                                        wOffset++;
                                    } else {
                                        exitg1 = 1;
                                    }
                                } else {
                                    idx->data[m] = b_iwork[ib];
                                    b_x->data[m] = b_xwork[ib];
                                    if (ib + 1 < bLen2) {
                                        ib++;
                                    } else {
                                        ib = m - wOffset;
                                        while (wOffset + 1 <= bLen) {
                                            idx->data[(ib + wOffset) + 1] = b_iwork[wOffset];
                                            b_x->data[(ib + wOffset) + 1] = b_xwork[wOffset];
                                            wOffset++;
                                        }

                                        exitg1 = 1;
                                    }
                                }
                            } while (exitg1 == 0);
                        }
                    }
                }

                m = nBlocks << 8;
                ib = nNonNaN - m;
                if (ib > 0) {
                    merge_block(idx, b_x, m, ib, 2, iwork, xwork);
                }

                m = 8;
            }
        }

        merge_block(idx, b_x, 0, nNonNaN, m, iwork, xwork);
    }

    if ((nNaNs > 0) && (nNonNaN > 0)) {
        for (k = 0; k + 1 <= nNaNs; k++) {
            xwork->data[k] = b_x->data[nNonNaN + k];
            iwork->data[k] = idx->data[nNonNaN + k];
        }

        for (k = nNonNaN - 1; k + 1 > 0; k--) {
            b_x->data[nNaNs + k] = b_x->data[k];
            idx->data[nNaNs + k] = idx->data[k];
        }

        for (k = 0; k + 1 <= nNaNs; k++) {
            b_x->data[k] = xwork->data[k];
            idx->data[k] = iwork->data[k];
        }
    }

    emxFree_real_T(&xwork);
    emxFree_int32_T(&iwork);
    wOffset = x->size[0];
    x->size[0] = b_x->size[0];
    emxEnsureCapacity((emxArray__common *)x, wOffset, (int)sizeof(double));
    m = b_x->size[0];
    for (wOffset = 0; wOffset < m; wOffset++) {
        x->data[wOffset] = b_x->data[wOffset];
    }

    emxFree_real_T(&b_x);
}

/*
 * File trailer for sortIdx.c
 *
 * [EOF]
 */



