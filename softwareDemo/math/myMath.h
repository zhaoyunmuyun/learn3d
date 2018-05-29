//
// Created by chen on 18-5-29.
//

#ifndef SOFTWAREDEMO_MYMATH_H
#define SOFTWAREDEMO_MYMATH_H


#include <cstring>

// defines for small numbers
#define EPSILON_E3 (float)(1E-3)
#define EPSILON_E4 (float)(1E-4)
#define EPSILON_E5 (float)(1E-5)
#define EPSILON_E6 (float)(1E-6)

//一、数据结构和类型
//1、向量和点
//不包含w分量的2D向量和点
typedef struct VECTOR2D_TYP
{
    union
    {
        float M[2];
        struct
        {
            float x,y;
        };
    };
}VECTOR2D,POINT2D,*VECTOR2D_PTR,*POINT2D_PTR;
//不包含w分量的3D向量和点
typedef struct VECTOR3D_TYP
{
    union
    {
        float M[3];
        struct
        {
            float x,y,z;
        };
    };
}VECTOR3D,POINT3D,*VECTOR3D_PTR,*POINT3D_PTR;
//包含w分量的4D齐次向量和点
typedef struct VECTOR4D_TYP
{
    union
    {
        float M[4];
        struct
        {
            float x,y,z,w;
        };
    };
}VECTOR4D,POINT4D,*VECTOR4D_PTR,*POINT4D_PTR;
//2D顶点
typedef struct VERTEX2DI_TYP
{
    int x,y;
}VERTEX2DI,*VERTEX2DI_PTR;
//2D顶点
typedef struct VERTEX2DF_TYP
{
    float x,y;
}VERTEX2DF,*VERTEX2DF_PTR;


//2、参数化直线
//2d参数化直线
typedef struct PARMLINE2D_TYP
{
    POINT2D p0;
    POINT2D p1;
    VECTOR2D v;
}PARMLINE2D,*PARMLINE2D_PTR;
//3d参数化直线
typedef struct PARMLINE3D_TYP
{
    POINT3D p0;
    POINT3D p1;
    VECTOR3D v;
}PARMLINE3D,*PARMLINE3D_PTR;


//3、3d平面
typedef struct PLANE3D_TYP
{
    POINT3D p0;
    VECTOR3D n;
}PLANE3D,*PLANE3D_PTR;


//4、矩阵
//4x4矩阵
typedef struct MATRIX4X4_TYP
{
    union
    {
        float M[4][4];
        struct
        {
            float M00,M01,M02,M03;
            float M10,M11,M12,M13;
            float M20,M21,M22,M23;
            float M30,M31,M32,M33;
        };
    };
}MATRIX4X4,*MATRIX4X4_PTR;
//4x3矩阵
typedef struct MATRIX4X3_TYP
{
    union
    {
        float M[4][3];
        struct
        {
            float M00,M01,M02;
            float M10,M11,M12;
            float M20,M21,M22;
            float M30,M31,M32;
        };
    };
}MATRIX4X3,*MATRIX4X3_PTR;
//1x4矩阵
typedef struct MATRIX1X4_TYP
{
    union
    {
        float M[4];
        struct
        {
            float M00,M01,M02,M03;
        };
    };
}MATRIX1X4,*MATRIX1X4_PTR;
//3x3矩阵
typedef struct MATRIX3X3_TYP
{
    union
    {
        float M[3][3];
        struct
        {
            float M00,M01,M02;
            float M10,M11,M12;
            float M20,M21,M22;
        };
    };
}MATRIX3X3,*MATRIX3X3_PTR;
//1x3矩阵
typedef struct MATRIX1X3_TYP
{
    union
    {
        float M[3];
        struct
        {
            float M00,M01,M02;
        };
    };
}MATRIX1X3,*MATRIX1X3_PTR;
//3x2矩阵
typedef struct MATRIX3X2_TYP
{
    union
    {
        float M[3][2];
        struct
        {
            float M00,M01;
            float M10,M11;
            float M20,M21;
        };
    };
}MATRIX3X2,*MATRIX3X2_PTR;
//2x2矩阵
typedef struct MATRIX2X2_TYP
{
    union
    {
        float M[2][2];
        struct
        {
            float M00,M01;
            float M10,M11;
        };
    };
}MATRIX2X2,*MATRIX2X2_PTR;
//1x2矩阵
typedef struct MATRIX1X2_TYP
{
    union
    {
        float M[2];
        struct
        {
            float M00,M01;
        };
    };
}MATRIX1X2,*MATRIX1X2_PTR;


//5、四元数
typedef struct QUAT_TYP
{
    union
    {
        float M[4];
        struct
        {
            float q0;
            VECTOR3D qv;
        };
        struct
        {
            float w,x,y,z;
        };
    };
}QUAT,*QUAT_PTR;

//6、角坐标系

//7、2d极坐标
typedef struct POLAR2D_TYP
{
    float r;
    float theta;
}POLAR2D,*POLAR2D_PTR;

//8、3d柱面坐标
typedef struct CYLINDRICAL3D_TYP
{
    float r;
    float theta;
    float z;
}CYLINDRICAL3D,*CYLINDRICAL3D_PTR;
//9、3d球面坐标
typedef struct SPHERICAL3D_TYP
{
    float p;
    float theta;
    float phi;
}SPHERICAL3D,*SPHERICAL3D_PTR;

//10、定点数
typedef int FIXP16;
typedef int *FIXP16_PTR;

//二、数学常量
//1、与PI相关的常量
#define PI          ((float)3.141592654f)
#define PI2         ((float)6.283185307f)
#define PI_DIV_2    ((float)1.570796327f)
#define PI_DIV_4    ((float)0.785398163f)
#define PI_INV      ((float)0.318309886f)
//2、与定点数运算相关的常量
#define FIXP16_SHIFT    16
#define FIXP16_MAG      65536
#define FIXP16_DP_MASK  0x0000ffff
#define FIXP16_WP_MASK  0xffff0000
#define FIXP16_ROUND_UP 0x00008000
//3、用于参数化直线交点的常量
#define PARM_LINE_NO_INTERSECT          0
#define PARM_LINE_INTERSECT_IN_SEGMENT  1
#define PARM_LINE_INTERSECT_OUT_SEGMENT 2
#define PARM_LINE_INTERSECT_EVERYWHERE  3
//4、单位矩阵
//4x4单位矩阵
const MATRIX4X4 IMAT_4X4=
        {
                1,0,0,0,
                0,1,0,0,
                0,0,1,0,
                0,0,0,1,
        };
//4x3单位矩阵
const MATRIX4X3 IMAT_4X3=
        {
                1,0,0,
                0,1,0,
                0,0,1,
                0,0,0,
        };
//3x3单位矩阵
const MATRIX3X3 IMAT_3X3=
        {
                1,0,0,
                0,1,0,
                0,0,1,
        };
//2x2单位矩阵
const MATRIX2X2 IMAT_2X2=
        {
                1,0,
                0,1,
        };

//三、宏和内联函数
//矩阵复制宏
#define MAT_COPY_2X2(src_mat, dest_mat) {memcpy((void *)(dest_mat), (void *)(src_mat), sizeof(MATRIX2X2) ); }
#define MAT_COPY_3X3(src_mat, dest_mat) {memcpy((void *)(dest_mat), (void *)(src_mat), sizeof(MATRIX3X3) ); }
#define MAT_COPY_4X4(src_mat, dest_mat) {memcpy((void *)(dest_mat), (void *)(src_mat), sizeof(MATRIX4X4) ); }
#define MAT_COPY_4X3(src_mat, dest_mat) {memcpy((void *)(dest_mat), (void *)(src_mat), sizeof(MATRIX4X3) ); }

//1、通用宏
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define SWAP(a,b,t) {t=a;a=b;b=t;}
#define DEG_TO_RAD(ang) ((ang)*PI/180.0)
#define RAD_TO_DEG(rads) ((rads)*180.0/PI)
#define RAND_RANGE(x,y) ((x)+((rand)()%((y)-(x)+1)))
//2、点和向量宏
//向量归零宏
inline void VECTOR2D_ZERO(VECTOR2D_PTR v)
{v->x=v->y=0.0;}
inline void VECTOR3D_ZERO(VECTOR3D_PTR v)
{v->x=v->y=v->z=0.0;}
inline void VECTOR4D_ZERO(VECTOR4D_PTR v)
{v->x=v->y=v->z=0.0;v->w=1.0;}
//使用分量初始化向量的宏
inline void VECTOR2D_INITXY(VECTOR2D_PTR v,float x,float y)
{v->x=x;v->y=y;}
inline void VECTOR3D_INITXY(VECTOR3D_PTR v,float x,float y, float z)
{v->x=x;v->y=y;v->z=z;}
inline void VECTOR4D_INITXY(VECTOR4D_PTR v,float x,float y, float z)
{v->x=x;v->y=y;v->z=z;v->w=1.0;}
//使用另一个向量来初始化向量的宏
inline void VECTOR2D_INIT(VECTOR2D_PTR vdst,VECTOR2D_PTR vsrc)
{vdst->x=vsrc->x;vdst->y=vsrc->y;}
inline void VECTOR3D_INIT(VECTOR3D_PTR vdst,VECTOR3D_PTR vsrc)
{vdst->x=vsrc->x;vdst->y=vsrc->y;vdst->z=vsrc->z;}
inline void VECTOR4D_INIT(VECTOR4D_PTR vdst,VECTOR4D_PTR vsrc)
{vdst->x=vsrc->x;vdst->y=vsrc->y;vdst->z=vsrc->z;vdst->w=vsrc->w;}
//复制向量的宏
inline void VECTOR2D_COPY(VECTOR2D_PTR vdst,VECTOR2D_PTR vsrc)
{vdst->x=vsrc->x;vdst->y=vsrc->y;}
inline void VECTOR3D_COPY(VECTOR3D_PTR vdst,VECTOR3D_PTR vsrc)
{vdst->x=vsrc->x;vdst->y=vsrc->y;vdst->z=vsrc->z;}
inline void VECTOR4D_COPY(VECTOR4D_PTR vdst,VECTOR4D_PTR vsrc)
{vdst->x=vsrc->x;vdst->y=vsrc->y;vdst->z=vsrc->z;vdst->w=vsrc->w;}
//初始化点的宏
inline void POINT2D_INIT(POINT2D_PTR vdst,POINT2D_PTR vsrc)
{vdst->x=vsrc->x;vdst->y=vsrc->y;}
inline void POINT3D_INI(POINT3D_PTR vdst,POINT3D_PTR vsrc)
{vdst->x=vsrc->x;vdst->y=vsrc->y;vdst->z=vsrc->z;}
inline void POINT4D_INI(POINT4D_PTR vdst,POINT4D_PTR vsrc)
{vdst->x=vsrc->x;vdst->y=vsrc->y;vdst->z=vsrc->z;vdst->w=vsrc->w;}
//复制点的宏
inline void POINT2D_COPY(POINT2D_PTR vdst,POINT2D_PTR vsrc)
{vdst->x=vsrc->x;vdst->y=vsrc->y;}
inline void POINT3D_COPY(POINT3D_PTR vdst,POINT3D_PTR vsrc)
{vdst->x=vsrc->x;vdst->y=vsrc->y;vdst->z=vsrc->z;}
inline void POINT4D_COPY(POINT4D_PTR vdst,POINT4D_PTR vsrc)
{vdst->x=vsrc->x;vdst->y=vsrc->y;vdst->z=vsrc->z;vdst->w=vsrc->w;}
//3、矩阵宏
//清空矩阵的宏
#define MAT_ZERO_2x2(m) {memset((void*)(m),0,sizeof(MATRIX2x2));}
#define MAT_ZERO_3x3(m) {memset((void*)(m),0,sizeof(MATRIX3x3));}
#define MAT_ZERO_4x4(m) {memset((void*)(m),0,sizeof(MATRIX4x4));}
#define MAT_ZERO_4x3(m) {memset((void*)(m),0,sizeof(MATRIX4x3));}
//设置单位矩阵的宏
#define MAT_IDENTITY_2x2(m) {memcpy((void *)(m),(void *)&IMAT_2x2,sizeof(MATRIX2x2));}
#define MAT_IDENTITY_3x3(m) {memcpy((void *)(m),(void *)&IMAT_3x3,sizeof(MATRIX3x3));}
#define MAT_IDENTITY_4x4(m) {memcpy((void *)(m),(void *)&IMAT_4x4,sizeof(MATRIX4x4));}
#define MAT_IDENTITY_4x3(m) {memcpy((void *)(m),(void *)&IMAT_4x3,sizeof(MATRIX4x3));}
//复制矩阵的宏
#define MAT_COPY_2x2(src_mat,dest_mat)  {memcpy((void *)(dest_mat),(void *)(src_mat),sizeof(MATRIX2x2));}
#define MAT_COPY_3x3(src_mat,dest_mat)  {memcpy((void *)(dest_mat),(void *)(src_mat),sizeof(MATRIX3x3));}
#define MAT_COPY_4x4(src_mat,dest_mat)  {memcpy((void *)(dest_mat),(void *)(src_mat),sizeof(MATRIX4x4));}
#define MAT_COPY_4x3(src_mat,dest_mat)  {memcpy((void *)(dest_mat),(void *)(src_mat),sizeof(MATRIX4x3));}
//对矩阵进行转置的宏
inline void MAT_TRANSPOSE_3X3(MATRIX3X3_PTR m)
{
    MATRIX3X3 mt;
    mt.M00=m->M00;mt.M01=m->M10;mt.M02=m->M20;
    mt.M10=m->M01;mt.M11=m->M11;mt.M12=m->M21;
    mt.M20=m->M02;mt.M21=m->M12;mt.M22=m->M22;
    memcpy((void *)m,(void *)&mt, sizeof(MATRIX3X3));
}
inline void MAT_TRANSPOSE_4X4(MATRIX4X4_PTR m)
{
    MATRIX4X4 mt;
    mt.M00=m->M00;mt.M01=m->M10;mt.M02=m->M20;mt.M03=m->M30;
    mt.M10=m->M01;mt.M11=m->M11;mt.M12=m->M21;mt.M13=m->M31;
    mt.M20=m->M02;mt.M21=m->M12;mt.M22=m->M22;mt.M23=m->M32;
    mt.M30=m->M03;mt.M31=m->M13;mt.M32=m->M23;mt.M33=m->M33;
    memcpy((void *)m,(void *)&mt, sizeof(MATRIX4X4));
}
inline void MAT_TRANSPOSE_3X3(MATRIX3X3_PTR m,MATRIX3X3_PTR mt)
{
    mt->M00=m->M00;mt->M01=m->M10;mt->M02=m->M20;
    mt->M10=m->M01;mt->M11=m->M11;mt->M12=m->M21;
    mt->M20=m->M02;mt->M21=m->M12;mt->M22=m->M22;
}
inline void MAT_TRANSPOSE_4X4(MATRIX4X4_PTR m,MATRIX4X4_PTR mt)
{
    mt->M00=m->M00;mt->M01=m->M10;mt->M02=m->M20;mt->M03=m->M30;
    mt->M10=m->M01;mt->M11=m->M11;mt->M12=m->M21;mt->M13=m->M31;
    mt->M20=m->M02;mt->M21=m->M12;mt->M22=m->M22;mt->M23=m->M32;
    mt->M30=m->M03;mt->M31=m->M13;mt->M32=m->M23;mt->M33=m->M33;
}
//矩阵和向量列互换宏
inline void MAT_COLUMN_SWAP_2X2(MATRIX2X2_PTR m,int c,MATRIX1X2_PTR v)
{m->M[0][c]=v->M[0];m->M[1][c]=v->M[1];}
inline void MAT_COLUMN_SWAP_3X3(MATRIX3X3_PTR m,int c,MATRIX1X3_PTR v)
{m->M[0][c]=v->M[0];m->M[1][c]=v->M[1];m->M[2][c]=v->M[2];}
inline void MAT_COLUMN_SWAP_4X4(MATRIX4X4_PTR m,int c,MATRIX1X4_PTR v)
{m->M[0][c]=v->M[0];m->M[1][c]=v->M[1];m->M[2][c]=v->M[2];m->M[3][c]=v->M[3];}
inline void MAT_COLUMN_SWAP_4X3(MATRIX4X3_PTR m,int c,MATRIX1X4_PTR v)
{m->M[0][c]=v->M[0];m->M[1][c]=v->M[1];m->M[2][c]=v->M[2];m->M[3][c]=v->M[3];}
//4、四元数
//四元数宏
inline void QUAT_ZERO(QUAT_PTR q)
{q->x=q->y=q->z=q->w=0.0;}
inline void QUAT_INITWXYZ(QUAT_PTR q,float w,float x,float y, float z)
{q->x=x;q->y=y;q->z=z;q->w=w;}
inline void QUAT_INIT_VECTOR3D(QUAT_PTR q,VECTOR3D_PTR v)
{q->x=v->x;q->y=v->y;q->z=v->z;q->w=0;}
inline void QUAT_INIT(QUAT_PTR qdst,QUAT_PTR qsrc)
{qdst->w=qsrc->w;qdst->x=qsrc->x;qdst->y=qsrc->y;qdst->z=qsrc->z;}
inline void QUAT_COPY(QUAT_PTR qdst,QUAT_PTR qsrc)
{qdst->w=qsrc->w;qdst->x=qsrc->x;qdst->y=qsrc->y;qdst->z=qsrc->z;}
//5、定点数宏
#define FIXP16_WP(fp) ((fp)>>FIXP16_SHIFT)
#define FIXP16_DP(fp) ((fp)&&FIXP16_DP_MASK)

#define INT_TO_FIXP16(i) ((i)>>FIXP16_SHIFT)
#define FLOAT_TO_FIXP16(f) ((float)(f)*(float)FIXP16_MAG+0.5)
//四函数原型
#define FIXP16_TO_FLOAT(fp) (((float)(fp))/FIXP16_MAG)


//五、函数原型
//通用三角函数
float Fast_Sin(float theta);
float Fast_Cos(float theta);

//距离函数
int Fast_Distance_2D(int x,int y);
float  Fast_Distance_3D(float x, float y, float z);

//极坐标、柱面坐标和球面坐标函数
//将一个用r和theta表示的2d极坐标点转换为一个(x,y)点，并将它存储在一个POINT2D结构中
void POLAR2D_To_POINT2D(POLAR2D_PTR polar, POINT2D_PTR rect);
//将一个2d极坐标转换为x，y坐标
void POLAR2D_To_RectXY(POLAR2D_PTR polar, float *x, float *y);
//将一个用直角坐标表示的点转换为2d极坐标格式
void POINT2D_To_POLAR2D(POINT2D_PTR rect, POLAR2D_PTR polar);
//将一个2d极坐标点转换为r和theta
void POINT2D_To_PolarRTh(POINT2D_PTR rect, float *r, float *theta);
//将一个3d柱面坐标转换为一个3d直角坐标系
void CYLINDRICAL3D_To_POINT3D(CYLINDRICAL3D_PTR cyl, POINT3D_PTR rect);
//将一个3d柱面坐标点转换为x、y、z坐标
void CYLINDRICAL3D_To_RectXYZ(CYLINDRICAL3D_PTR cyl, float *x, float *y, float *z);
//将3d直角坐标点转换为3d柱面坐标点
void POINT3D_To_CYLINDRICAL3D(POINT3D_PTR rect, CYLINDRICAL3D_PTR cyl);
//将3d点转换为柱面坐标r、theta、z
void POINT3D_To_CylindricalRThZ(POINT3D_PTR rect, float *r, float *theta, float *z);
//将一个3d球面坐标转换为3d直角坐标点
void SPHERICAL3D_To_POINT3D(SPHERICAL3D_PTR sph, POINT3D_PTR rect);
//将一个3d球面坐标点转换为x、y、z坐标
void SPHERICAL3D_To_RectXYZ(SPHERICAL3D_PTR sph, float *x, float *y, float *z);
//将一个3d点转换为球面坐标点
void POINT3D_To_SPHERICAL3D(POINT3D_PTR rect, SPHERICAL3D_PTR sph);
//将一个3d点转换为球面坐标rho、theta和phi
void POINT3D_To_SphericalPThPh(POINT3D_PTR rect, float *p, float *theta, float *phi);

//2d向量函数
//指定向量va和vb相加，然后将结果存储在vsum中
void VECTOR2D_Add(VECTOR2D_PTR va, VECTOR2D_PTR vb, VECTOR2D_PTR vsum);
//向量相加函数的堆栈版本，将结果返回到堆栈中
VECTOR2D VECTOR2D_Add(VECTOR2D_PTR va, VECTOR2D_PTR vb);
//将两个向量相减，并将结果存储到vdiff中
void VECTOR2D_Sub(VECTOR2D_PTR va, VECTOR2D_PTR vb, VECTOR2D_PTR vdiff);
//向量相减函数的堆栈版本，将结果返回到堆栈
VECTOR2D VECTOR2D_Sub(VECTOR2D_PTR va, VECTOR2D_PTR vb);
//使用缩放因子k对向量va进行缩放，然后将结果(k*va)存储在va中
void VECTOR2D_Scale(float k, VECTOR2D_PTR va);
//使用缩放因子k对向量va进行缩放，然后将结果(k*va)存储在vscaled中
void VECTOR2D_Scale(float k, VECTOR2D_PTR va, VECTOR2D_PTR vscaled);
//计算点积，并返回标量结果
float VECTOR2D_Dot(VECTOR2D_PTR va, VECTOR2D_PTR vb);
//使用标准的“平方和的平方根”算法计算指定向量的长度
float VECTOR2D_Length(VECTOR2D_PTR va);
//使用泰勒技术近似计算制定向量的长度
float VECTOR2D_Length_Fast(VECTOR2D_PTR va);
//将向量va归一化，即将向量va的各个分量除以其长度，使之成为一个单位向量
void VECTOR2D_Normalize(VECTOR2D_PTR va);
//将向量va归一化，并将结果存储在vn中
void VECTOR2D_Normalize(VECTOR2D_PTR va, VECTOR2D_PTR vn);
//创建一个从init到term的向量，并将它存储到result中
void VECTOR2D_Build(VECTOR2D_PTR init, VECTOR2D_PTR term, VECTOR2D_PTR result);
//计算两个向量va和vb之间的夹角的余弦
float VECTOR2D_CosTh(VECTOR2D_PTR va, VECTOR2D_PTR vb);

//3d向量函数
void VECTOR3D_Add(VECTOR3D_PTR va, VECTOR3D_PTR vb, VECTOR3D_PTR vsum);
VECTOR3D VECTOR3D_Add(VECTOR3D_PTR va, VECTOR3D_PTR vb);
void VECTOR3D_Sub(VECTOR3D_PTR va, VECTOR3D_PTR vb, VECTOR3D_PTR vdiff);
VECTOR3D VECTOR3D_Sub(VECTOR3D_PTR va, VECTOR3D_PTR vb);
void VECTOR3D_Scale(float k, VECTOR3D_PTR va);
void VECTOR3D_Scale(float k, VECTOR3D_PTR va, VECTOR3D_PTR vscaled);
float VECTOR3D_Dot(VECTOR3D_PTR va, VECTOR3D_PTR vb);
//计算向量叉积，即与向量和vb都垂直的向量，并将结果范湖掉堆栈中
void VECTOR3D_Cross(VECTOR3D_PTR va,VECTOR3D_PTR vb,VECTOR3D_PTR vn);
VECTOR3D VECTOR3D_Cross(VECTOR3D_PTR va, VECTOR3D_PTR vb);
float VECTOR3D_Length(VECTOR3D_PTR va);
float VECTOR3D_Length_Fast(VECTOR3D_PTR va);
void VECTOR3D_Normalize(VECTOR3D_PTR va);
void VECTOR3D_Normalize(VECTOR3D_PTR va, VECTOR3D_PTR vn);
void VECTOR3D_Build(VECTOR3D_PTR init, VECTOR3D_PTR term, VECTOR3D_PTR result);
float VECTOR3D_CosTh(VECTOR3D_PTR va, VECTOR3D_PTR vb);

//4d向量函数
void VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vsum);
VECTOR4D VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb);
void VECTOR4D_Sub(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vdiff);
VECTOR4D VECTOR4D_Sub(VECTOR4D_PTR va, VECTOR4D_PTR vb);
void VECTOR4D_Scale(float k, VECTOR4D_PTR va);
void VECTOR4D_Scale(float k, VECTOR4D_PTR va, VECTOR4D_PTR vscaled);
float VECTOR4D_Dot(VECTOR4D_PTR va, VECTOR4D_PTR vb);
void VECTOR4D_Cross(VECTOR4D_PTR va,VECTOR4D_PTR vb,VECTOR4D_PTR vn);
VECTOR4D VECTOR4D_Cross(VECTOR4D_PTR va, VECTOR4D_PTR vb);
float VECTOR4D_Length(VECTOR4D_PTR va);
float VECTOR4D_Length_Fast(VECTOR4D_PTR va);
void VECTOR4D_Normalize(VECTOR4D_PTR va);
void VECTOR4D_Normalize(VECTOR4D_PTR va, VECTOR4D_PTR vn);
void VECTOR4D_Build(VECTOR4D_PTR init, VECTOR4D_PTR term, VECTOR4D_PTR result);
float VECTOR4D_CosTh(VECTOR4D_PTR va, VECTOR4D_PTR vb);

//2x2矩阵函数
//使用制定的浮点数按先行后列的顺序初始化矩阵ma
void Mat_Init_2X2(MATRIX2X2_PTR ma, float m00, float m01, float m10, float m11);
//计算矩阵m的行列式，并将结果返回到堆栈
float Mat_Det_2X2(MATRIX2X2_PTR m);
//将两个矩阵相加，并将结果存储到msum中
void Mat_Add_2X2(MATRIX2X2_PTR ma, MATRIX2X2_PTR mb, MATRIX2X2_PTR msum);
//将两个矩阵相乘，并将结果存储到mprod中
void Mat_Mul_2X2(MATRIX2X2_PTR ma, MATRIX2X2_PTR mb, MATRIX2X2_PTR mprod);
//计算矩阵m的逆矩阵，并将结果存储到mi中，如果逆矩阵存在，该函数返回1，否则返回0，且mi未定义
int Mat_Inverse_2X2(MATRIX2X2_PTR m, MATRIX2X2_PTR mi);
//
int Solve_2X2_System(MATRIX2X2_PTR A, MATRIX1X2_PTR X, MATRIX1X2_PTR B);

//3x3矩阵函数
int Mat_Mul_1x2_3x2(MATRIX1X2_PTR ma,MATRIX3X2_PTR mb,MATRIX1X2_PTR mprod);
//将一个1x3矩阵与一个3x3矩阵相乘
int Mat_Mul_1x3_3x3(MATRIX1X3_PTR ma,MATRIX3X3_PTR mb,MATRIX1X3_PTR mprod);
//将两个3x3矩阵相乘,并将结果存储到mprod中
int Mat_Mul_3x3(MATRIX3X3_PTR ma,MATRIX3X3_PTR mb,MATRIX3X3_PTR mprod);
//使用传入的浮点数以先行后列的顺序初始化矩阵ma
int Mat_Init_3x2(MATRIX3X2_PTR ma,
                 float m00,float m01,
                 float m10,float m11,
                 float m20,float m21);
//将两个矩阵相加，并将结果存储到msum中
void Mat_Add_3X3(MATRIX3X3_PTR ma, MATRIX3X3_PTR mb, MATRIX3X3_PTR msum);
//将1x3行向量va与3x3矩阵mb相乘，并将结果存储到1x3的行向量vprod中
void Mat_Mul_VECTOR3D_3X3(VECTOR3D_PTR  va, MATRIX3X3_PTR mb,VECTOR3D_PTR  vprod);
int Mat_Inverse_3X3(MATRIX3X3_PTR m, MATRIX3X3_PTR mi);
//使用传入的浮点数以先行后列的顺序初始化矩阵ma
void Mat_Init_3X3(MATRIX3X3_PTR ma,
                  float m00, float m01, float m02,
                  float m10, float m11, float m12,
                  float m20, float m21, float m22);
//计算m的行列式，并将结果返回到堆栈
float Mat_Det_3X3(MATRIX3X3_PTR m);
int Solve_3X3_System(MATRIX3X3_PTR A, MATRIX1X3_PTR X, MATRIX1X3_PTR B);

//4x4矩阵函数
//将两个矩阵相加，并将结果存储到msum中
void Mat_Add_4X4(MATRIX4X4_PTR ma, MATRIX4X4_PTR mb, MATRIX4X4_PTR msum);
//将两个矩阵相乘，并将结果存储到mprod中
void Mat_Mul_4X4(MATRIX4X4_PTR ma, MATRIX4X4_PTR mb, MATRIX4X4_PTR mprod);
//将一个1x4的行向量与一个4x4矩阵相乘，并将结果存储到mprod中
void Mat_Mul_1X4_4X4(MATRIX1X4_PTR ma, MATRIX4X4_PTR mb, MATRIX1X4_PTR mprod);
//将一个3d向量与一个4x4矩阵相乘
void Mat_Mul_VECTOR3D_4X4(VECTOR3D_PTR  va, MATRIX4X4_PTR mb, VECTOR3D_PTR  vprod);
//将一个3d向量与一个4x3矩阵相乘
void Mat_Mul_VECTOR3D_4X3(VECTOR3D_PTR  va, MATRIX4X3_PTR mb, VECTOR3D_PTR  vprod);
//将一个1x4行向量va与4x4矩阵相乘，并将结果存储到1x4行向量vprod中
void Mat_Mul_VECTOR4D_4X4(VECTOR4D_PTR  va, MATRIX4X4_PTR mb, VECTOR4D_PTR  vprod);
//将一个4d向量与一个4x3矩阵相乘
void Mat_Mul_VECTOR4D_4X3(VECTOR4D_PTR  va, MATRIX4X4_PTR mb, VECTOR4D_PTR  vprod);
//用于计算矩阵m的逆矩阵
int Mat_Inverse_4X4(MATRIX4X4_PTR m, MATRIX4X4_PTR mi);
//使用传入的浮点数值以先行后列的顺序初始化矩阵ma
void Mat_Init_4X4(MATRIX4X4_PTR ma,
                  float m00, float m01, float m02, float m03,
                  float m10, float m11, float m12, float m13,
                  float m20, float m21, float m22, float m23,
                  float m30, float m31, float m32, float m33);

//四元数函数
//将两个四元数q1和q2相加，并将结果存储到qsum中
void QUAT_Add(QUAT_PTR q1, QUAT_PTR q2, QUAT_PTR qsum);
//将两个四元数q1和q2相减，并将结果存储到qdiff中
void QUAT_Sub(QUAT_PTR q1, QUAT_PTR q2, QUAT_PTR qdiff);
//计算四元数q的共轭，并将结果存储到qconj中
void QUAT_Conjugate(QUAT_PTR q, QUAT_PTR qconj);
//根据缩放因子scale对四元数q进行缩放，并将结果存储到qs中
void QUAT_Scale(QUAT_PTR q, float scale, QUAT_PTR qs);
//根据缩放因子直接缩放四元数q,即修改q
void QUAT_Scale(QUAT_PTR q, float scale);
//返回四元数q的范数，即长度
float QUAT_Norm(QUAT_PTR q);
//返回四元数q的范数平方，即长度平方
float QUAT_Norm2(QUAT_PTR q);
//将四元数q归一化，并将结果存储到qn中
void QUAT_Normalize(QUAT_PTR q, QUAT_PTR qn);
//直接对四元数q进行归一化，即修改q
void QUAT_Normalize(QUAT_PTR q);
//计算四元数q的逆，并将结果存储到qi中
void QUAT_Unit_Inverse(QUAT_PTR q, QUAT_PTR qi);
//直接计算四元数q的逆
void QUAT_Unit_Inverse(QUAT_PTR q);
//计算非单元四元数q的逆，将结果存储到qi中
void QUAT_Inverse(QUAT_PTR q, QUAT_PTR qi);
//计算非单元四元数q的逆
void QUAT_Inverse(QUAT_PTR q);
//将四元数相乘，并将结果存储到qprod中
void QUAT_Mul(QUAT_PTR q1, QUAT_PTR q2, QUAT_PTR qprod);
//将3ge四元数相乘，并将结果存储到qprod中
void QUAT_Triple_Product(QUAT_PTR q1, QUAT_PTR q2, QUAT_PTR q3, QUAT_PTR qprod);
//根据方向向量v和角度theta创建一个旋转四元数
void VECTOR3D_Theta_To_QUAT(QUAT_PTR q, VECTOR3D_PTR v, float theta);
void VECTOR4D_Theta_To_QUAT(QUAT_PTR q, VECTOR4D_PTR v, float theta);
//根据绕z、y、x旋转的欧拉角创建一个旋转四元数
void EulerZYX_To_QUAT(QUAT_PTR q, float theta_z, float theta_y, float theta_x);
//将一个单位旋转四元数转换为一个单位3d向量和一个绕该向量旋转的角度theta
void QUAT_To_VECTOR3D_Theta(QUAT_PTR q, VECTOR3D_PTR v, float *theta);

//2d参数化直线函数
//根据指定的点计算他们之间的向量，并初始化一条2d参数化直线
void Init_Parm_Line2D(POINT2D_PTR p_init, POINT2D_PTR p_term, PARMLINE2D_PTR p);
//计算2d参数化直线在参数t处的值，并将其返回到制定的点中
void Compute_Parm_Line2D(PARMLINE2D_PTR p, float t, POINT2D_PTR pt);
//计算两条参数化直线p1和p2的交点，并将交点对应的参数值t1和t2分别存储到相应的变量中
int Intersect_Parm_Lines2D(PARMLINE2D_PTR p1, PARMLINE2D_PTR p2, float *t1, float *t2);
//计算参数化直线的交点，但不返回交点对应的参数t1和t2的值，而是返回交点
int Intersect_Parm_Lines2D(PARMLINE2D_PTR p1, PARMLINE2D_PTR p2, POINT2D_PTR pt);

//3d参数化直线函数
//根据制定的点以及它们之间的向量，初始化一个3d参数化直线结构
void Init_Parm_Line3D(POINT3D_PTR p_init, POINT3D_PTR p_term, PARMLINE3D_PTR p);
//计算一条参数化直线在参数t处的值，并将其返回存储到制定的点中
void Compute_Parm_Line3D(PARMLINE3D_PTR p, float t, POINT3D_PTR pt);

//3d平面函数
//使用指定的点和法线来初始化一个3d平面
void PLANE3D_Init(PLANE3D_PTR plane, POINT3D_PTR p0,
                  VECTOR3D_PTR normal, int normalize);
//判断制定点位于哪个半空间
float Compute_Point_In_Plane3D(POINT3D_PTR pt, PLANE3D_PTR plane);
//计算一条3d参数化直线与一个3d平面的交点，将交点处的参数值存储到t中，并将交点存储到pt中
int Intersect_Parm_Line3D_Plane3D(PARMLINE3D_PTR pline, PLANE3D_PTR plane,
                                  float *t, POINT3D_PTR pt);
//定点数函数
FIXP16 FIXP16_MUL(FIXP16 fp1, FIXP16 fp2);
FIXP16 FIXP16_DIV(FIXP16 fp1, FIXP16 fp2);

//六、全局变量
extern float cos_look[361];
extern float sin_look[361];
void Build_Sin_Cos_Tables();
#endif //SOFTWAREDEMO_MYMATH_H
