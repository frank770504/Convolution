#include "stdio.h"
#include "stdlib.h"

#define MAXC(a,b) ((a>b)?a:b)
#define MINC(a,b) ((a<b)?a:b)
#define ABSC(a) ((a>=0)?a:a*-1)

#define SIZE_CVBUFF       25
#define GET_WINDOW_SIZE(wid) sizeof(wid)/2

static int CV_Buf[SIZE_CVBUFF]; //convolution buffer

    //-Filter Coefficient-//
const short Ma_walk_coeff[]   = {5, 9, 14, 19, 23, 28, 33, 37, 42, 47}; // verse [10:1]./55
const short Ma_uni10_coeff[]   = {25, 25, 25, 25, 25, 25, 25, 25, 25, 26}; // 1/10
const short Ma_uni5_coeff[]   = {51, 51, 51, 51, 52}; // 1/10
const short Ma_uni20_coeff[]   = {12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
                                 13, 13, 13, 13, 14, 14, 14, 14, 14, 14}; // 1/10
const short Ma_run_coeff[]   = {0, 0, 0, 0, 0, 17, 34, 51, 64, 71}; // verse [10:1]./55
//const SINT16 Ma_run_coeff[]   = {0, 0, 0, 0, 0, 0, 0, 42, 85, 109}; // verse [10:1]./55

const short Diff_run_coeff[] = {-2, 2, 1, 0, -1}; // >>[2 1 0 1 2] (j)
const short Diff_walk_coeff[] = {-4, 4, 2, 1, 0, -1, -2}; // >>[2 1 0 1 2] (j)
const short Diff_uni_coeff[] = {-73, -37, -18, 0, 18, 37, 73}; // >>[2 1 0 1 2] (j)
const short Diff_uni2_coeff[] = {-68, -34, -17, -9, 0, 9, 17, 34, 68}; // >>[2 1 0 1 2] (j)
const short Diff_uni3_coeff[] = {-66, -33, -17, -8, -4, 0, 4, 8, 17, 33, 66}; // >>[2 1 0 1 2] (j)

typedef int (*routine_fp)(int, uint32_t*, const short*, int);

void reset_cv_buf(int val)
// fill the convolution buffer by using a value
// CV_Buf is global variable
{
    int i;
    for(i=0;i<SIZE_CVBUFF;i++)
        CV_Buf[i] = val;

    return; //void
}

int avg_routine(int st_pt, uint32_t* ptrData, const short* coeff, int Round)
{
    int i;
    int j;
    int Avg_result;

    Avg_result = 0;
    j = 0;

    for(i = 0; i < Round; i++)
    {
            j = (i + st_pt) % Round; //check
            Avg_result = Avg_result + coeff[i] * *(ptrData + j);
    }
    Avg_result = Avg_result>>8;
    return Avg_result;
}

int subst_routine(int st_pt, uint32_t* ptrData, const short* coeff, int Round)
{
    int i;
    int j, head;
    int Avg_result;

    Avg_result = 0;
    j = 0;

    for(i = 0; i < Round; i++)
    {
        j = (i + st_pt) % Round; //check
        Avg_result = Avg_result + coeff[i] * *(ptrData + j);
        if (i==Round-1) { head = *(ptrData + j); }
    }
    Avg_result = Avg_result>>8;
    return head - (Avg_result>>1);
}

void convolution(int* souc, int* dest, int len, const short* coeff, int Round, routine_fp routine)
{
    int* ptrData32;
    int   i, j, ind;

    reset_cv_buf(0);
    ptrData32 = CV_Buf;
    j = 0;
    for(i=0;i<len;i++)
    {
        *(ptrData32 + j%Round) =  *(souc + i);
//        printf("%d..%d\n", i, (j+1)%Round);
        *(dest + i) = routine(j+1, CV_Buf, coeff, Round);
        j++;
    }
}

#define SIZE_X 2000

void gvtless(int* acc, int len)
{
    convolution(acc, acc, SIZE_X, Ma_uni10_coeff, GET_WINDOW_SIZE(Ma_uni10_coeff), subst_routine);
}

int sqrt_q8(int a)
// Newton method
// x1 = ((a/x0)+x0)/2
{
    int i, x = 128, temp, last = 0, err = 1;

    while(ABSC(last-x)>1)
    {
        last = x;
        temp = (a<<8)/x;
        temp = temp + x;
        x = temp>>1;
    }
    return x;
}

int main()
{
    printf("start...\n");
    int x, y, z;
    int i = 0;
    int accx[SIZE_X], accy[SIZE_X], accz[SIZE_X];
    int ma[SIZE_X], diff[SIZE_X];

    FILE* f;
    f = fopen("./PedoAcc.txt","r");
    if (f==NULL) { printf("file read failed \n"); }
    while (fscanf(f,"%d %d %d\n", &x, &y, &z)!=-1)
    {
        accx[i] = x;
        accy[i] = y;
        accz[i] = z;
        i++;
    }
    fclose(f);

    convolution(accy, ma, SIZE_X, Ma_uni10_coeff, GET_WINDOW_SIZE(Ma_uni10_coeff), avg_routine);
    convolution(ma, diff, SIZE_X, Diff_uni3_coeff, GET_WINDOW_SIZE(Diff_uni3_coeff), avg_routine);

    return 0;
}

