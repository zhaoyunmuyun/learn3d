//
// Created by chen on 18-5-29.
//
#include <cmath>
#include <iostream>
#include <cstdarg>
#include "myMath.h"

float cos_look[361];
float sin_look[361];

float Fast_Sin(float theta)
{
    //将角度转换为0-359的值
    theta=fmodf(theta,360);
    //将角度转为正值
    if(theta<0) theta+=360.0;
    //提取角度的整数部分和小数部分，以便进行插值计算
    auto theta_int=(int)theta;
    float theta_frac=theta-theta_int;

    return (sin_look[theta_int]+
            theta_frac*(sin_look[theta_int+1]-sin_look[theta_int]));

}
float Fast_Cos(float theta)
{
    theta = fmodf(theta,360);

    if (theta < 0) theta+=360.0;

    auto theta_int    = (int)theta;
    float theta_frac = theta - theta_int;

    return(cos_look[theta_int] +
           theta_frac*(cos_look[theta_int+1] - cos_look[theta_int]));
}

int Fast_Distance_2D(int x, int y)
{
    x = abs(x);
    y = abs(y);

    int mn = MIN(x,y);

    return(x+y-(mn>>1)-(mn>>2)+(mn>>4));

}

float Fast_Distance_3D(float fx, float fy, float fz)
{
    int temp;
    int x,y,z;

    x = static_cast<int>(std::fabs(fx) * 1024);
    y = static_cast<int>(std::fabs(fy) * 1024);
    z = static_cast<int>(std::fabs(fz) * 1024);

    if (y < x) SWAP(x,y,temp)

    if (z < y) SWAP(y,z,temp)

    if (y < x) SWAP(x,y,temp)

    int dist = (z + 11 * (y >> 5) + (x >> 2) );

    return((float)(dist >> 10));
}

////////////////////////////////////////////////////////////

void POLAR2D_To_POINT2D(POLAR2D_PTR polar, POINT2D_PTR rect)
{
    rect->x = polar->r*cosf(polar->theta);
    rect->y = polar->r*sinf(polar->theta);

}
void POLAR2D_To_RectXY(POLAR2D_PTR polar, float *x, float *y)
{
    *x = polar->r*cosf(polar->theta);
    *y = polar->r*sinf(polar->theta);

}
void POINT2D_To_POLAR2D(POINT2D_PTR rect, POLAR2D_PTR polar)
{
    polar->r     = sqrtf((rect->x * rect->x) + (rect->y * rect->y));
    polar->theta = atanf(rect->y/rect->x);

}
void POINT2D_To_PolarRTh(POINT2D_PTR rect, float *r, float *theta)
{
// convert rectangular to polar
    *r=sqrtf((rect->x * rect->x) + (rect->y * rect->y));
    *theta = atanf(rect->y/rect->x);

}
void CYLINDRICAL3D_To_POINT3D(CYLINDRICAL3D_PTR cyl, POINT3D_PTR rect)
{
    rect->x = cyl->r*cosf(cyl->theta);
    rect->y = cyl->r*sinf(cyl->theta);
    rect->z = cyl->z;

}
void CYLINDRICAL3D_To_RectXYZ(CYLINDRICAL3D_PTR cyl,
                              float *x, float *y, float *z)
{
    *x = cyl->r*cosf(cyl->theta);
    *y = cyl->r*sinf(cyl->theta);
    *z = cyl->z;

}
void POINT3D_To_CYLINDRICAL3D(POINT3D_PTR rect,
                              CYLINDRICAL3D_PTR cyl)
{
    cyl->r     = sqrtf((rect->x * rect->x) + (rect->y * rect->y));
    cyl->theta = atanf(rect->y/rect->x);
    cyl->z     = rect->z;

}
void POINT3D_To_CylindricalRThZ(POINT3D_PTR rect,
                                float *r, float *theta, float *z)
{
    *r     = sqrtf((rect->x * rect->x) + (rect->y * rect->y));
    *theta = atanf(rect->y/rect->x);
    *z     = rect->z;

}
void SPHERICAL3D_To_POINT3D(SPHERICAL3D_PTR sph, POINT3D_PTR rect)
{
    float r;

    r       = sph->p*sinf(sph->phi);
    rect->z = sph->p*cosf(sph->phi);

    rect->x = r*cosf(sph->theta);
    rect->y = r*sinf(sph->theta);

}
void SPHERICAL3D_To_RectXYZ(SPHERICAL3D_PTR sph,
                            float *x, float *y, float *z)
{
    float r;

    r  = sph->p*sinf(sph->phi);
    *z = sph->p*cosf(sph->phi);

    *x = r*cosf(sph->theta);
    *y = r*sinf(sph->theta);
}
void POINT3D_To_SPHERICAL3D(POINT3D_PTR rect, SPHERICAL3D_PTR sph)
{
    sph->p = sqrtf((rect->x*rect->x)+(rect->y*rect->y)+(rect->z*rect->z));

    sph->theta = atanf(rect->y/rect->x);

    float r = sph->p*sinf(sph->phi);

    sph->phi   = asinf(r/sph->p);
}
void POINT3D_To_SphericalPThPh(POINT3D_PTR rect,
                               float *p, float *theta, float *phi)
{
    *p     = sqrtf((rect->x*rect->x)+(rect->y*rect->y)+(rect->z*rect->z));
    *theta = atanf(rect->y/rect->x);

    float r = sqrtf((rect->x * rect->x) + (rect->y * rect->y));
    *phi    = asinf(r / (*p));

}
///////////////////////////////////////////////////////

void VECTOR2D_Add(VECTOR2D_PTR va,
                  VECTOR2D_PTR vb,
                  VECTOR2D_PTR vsum)
{
    vsum->x = va->x + vb->x;
    vsum->y = va->y + vb->y;

}
VECTOR2D VECTOR2D_Add(VECTOR2D_PTR va,
                      VECTOR2D_PTR vb)
{
    VECTOR2D vsum;

    vsum.x = va->x + vb->x;
    vsum.y = va->y + vb->y;

    return(vsum);

}
void VECTOR2D_Sub(VECTOR2D_PTR va,
                  VECTOR2D_PTR vb,
                  VECTOR2D_PTR vdiff)
{
    vdiff->x = va->x - vb->x;
    vdiff->y = va->y - vb->y;

}
VECTOR2D VECTOR2D_Sub(VECTOR2D_PTR va,
                      VECTOR2D_PTR vb)
{
    VECTOR2D vdiff;

    vdiff.x = va->x - vb->x;
    vdiff.y = va->y - vb->y;

    return(vdiff);

}
void VECTOR2D_Scale(float k,
                    VECTOR2D_PTR va,
                    VECTOR2D_PTR vscaled)
{
    vscaled->x = k*va->x;
    vscaled->y = k*va->y;
}
void VECTOR2D_Scale(float k, VECTOR2D_PTR va)
{
    va->x*=k;
    va->y*=k;
}
float VECTOR2D_Dot(VECTOR2D_PTR va, VECTOR2D_PTR vb)
{
    return( (va->x * vb->x) + (va->y * vb->y) );
}
float VECTOR2D_Length(VECTOR2D_PTR va)
{
    return(sqrtf(va->x*va->x + va->y*va->y));
}
float VECTOR2D_Length_Fast(VECTOR2D_PTR va)
{
    return( (float)Fast_Distance_2D(va->x, va->y) );

}
void VECTOR2D_Normalize(VECTOR2D_PTR va)
{
    float length = sqrtf(va->x*va->x + va->y*va->y );

    if (length < EPSILON_E5)
        return;

    float length_inv = 1/length;

    va->x = va->x*length_inv;
    va->y = va->y*length_inv;

}
void VECTOR2D_Normalize(VECTOR2D_PTR va, VECTOR2D_PTR vn)
{
    VECTOR2D_ZERO(vn);

    auto length = (float)sqrtf(va->x*va->x + va->y*va->y );

    if (length < EPSILON_E5)
        return;

    float length_inv = 1/length;

    vn->x = va->x*length_inv;
    vn->y = va->y*length_inv;

}
void VECTOR2D_Build(VECTOR2D_PTR init,
                    VECTOR2D_PTR term,
                    VECTOR2D_PTR result)
{

    result->x = term->x - init->x;
    result->y = term->y - init->y;

}
float VECTOR2D_CosTh(VECTOR2D_PTR va, VECTOR2D_PTR vb)
{
    return(VECTOR2D_Dot(va,vb)/(VECTOR2D_Length(va)*VECTOR2D_Length(vb)));
}
void VECTOR3D_Add(VECTOR3D_PTR va,
                  VECTOR3D_PTR vb,
                  VECTOR3D_PTR vsum)
{
    vsum->x = va->x + vb->x;
    vsum->y = va->y + vb->y;
    vsum->z = va->z + vb->z;
}
VECTOR3D VECTOR3D_Add(VECTOR3D_PTR va,
                      VECTOR3D_PTR vb)
{
    VECTOR3D vsum;

    vsum.x = va->x + vb->x;
    vsum.y = va->y + vb->y;
    vsum.z = va->z + vb->z;

    return(vsum);
}
void VECTOR3D_Sub(VECTOR3D_PTR va,
                  VECTOR3D_PTR vb,
                  VECTOR3D_PTR vdiff)
{
    vdiff->x = va->x - vb->x;
    vdiff->y = va->y - vb->y;
    vdiff->z = va->z - vb->z;
}
VECTOR3D VECTOR3D_Sub(VECTOR3D_PTR va, VECTOR3D_PTR vb)
{
    VECTOR3D vdiff;

    vdiff.x = va->x - vb->x;
    vdiff.y = va->y - vb->y;
    vdiff.z = va->z - vb->z;

    return(vdiff);
}
void VECTOR3D_Scale(float k, VECTOR3D_PTR va)
{
    va->x*=k;
    va->y*=k;
    va->z*=k;

}
void VECTOR3D_Scale(float k, VECTOR3D_PTR va, VECTOR3D_PTR vscaled)
{
    vscaled->x = k*va->x;
    vscaled->y = k*va->y;
    vscaled->z = k*va->z;

}
float VECTOR3D_Dot(VECTOR3D_PTR va, VECTOR3D_PTR vb)
{
    return( (va->x * vb->x) + (va->y * vb->y) + (va->z * vb->z) );
}
void VECTOR3D_Cross(VECTOR3D_PTR va, VECTOR3D_PTR vb, VECTOR3D_PTR vn)
{
    vn->x =  ( (va->y * vb->z) - (va->z * vb->y) );
    vn->y = -( (va->x * vb->z) - (va->z * vb->x) );
    vn->z =  ( (va->x * vb->y) - (va->y * vb->x) );
}
VECTOR3D VECTOR3D_Cross(VECTOR3D_PTR va, VECTOR3D_PTR vb)
{
    VECTOR3D vn;

    vn.x =  ( (va->y * vb->z) - (va->z * vb->y) );
    vn.y = -( (va->x * vb->z) - (va->z * vb->x) );
    vn.z =  ( (va->x * vb->y) - (va->y * vb->x) );

    return(vn);
}
float VECTOR3D_Length(VECTOR3D_PTR va)
{
    return( (float)sqrtf(va->x*va->x + va->y*va->y + va->z*va->z) );
}
float VECTOR3D_Length_Fast(VECTOR3D_PTR va)
{
    return( Fast_Distance_3D(va->x, va->y, va->z) );

}
void VECTOR3D_Normalize(VECTOR3D_PTR va)
{
    float length = sqrtf(va->x*va->x + va->y*va->y + va->z*va->z);

    if (length < EPSILON_E5)
        return;

    float length_inv = 1/length;

    va->x*=length_inv;
    va->y*=length_inv;
    va->z*=length_inv;

}
void VECTOR3D_Normalize(VECTOR3D_PTR va, VECTOR3D_PTR vn)
{
    VECTOR3D_ZERO(vn);

    float length = VECTOR3D_Length(va);

    if (length < EPSILON_E5)
        return;

    float length_inv = 1.0/length;

    vn->x = va->x*length_inv;
    vn->y = va->y*length_inv;
    vn->z = va->z*length_inv;
}
void VECTOR3D_Build(VECTOR3D_PTR init,VECTOR3D_PTR term,VECTOR3D_PTR result)
{
    result->x = term->x - init->x;
    result->y = term->y - init->y;
    result->z = term->z - init->z;
}
float VECTOR3D_CosTh(VECTOR3D_PTR va, VECTOR3D_PTR vb)
{
    return(VECTOR3D_Dot(va,vb)/(VECTOR3D_Length(va)*VECTOR3D_Length(vb)));

}
///////////////////////////////////////////////////

void VECTOR4D_Build(VECTOR4D_PTR init, VECTOR4D_PTR term, VECTOR4D_PTR result)
{
    result->x = term->x - init->x;
    result->y = term->y - init->y;
    result->z = term->z - init->z;
    result->w = 1;

}
void VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vsum)
{
    vsum->x = va->x + vb->x;
    vsum->y = va->y + vb->y;
    vsum->z = va->z + vb->z;
    vsum->w = 1;

}
VECTOR4D VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb)
{
    VECTOR4D vsum;

    vsum.x = va->x + vb->x;
    vsum.y = va->y + vb->y;
    vsum.z = va->z + vb->z;
    vsum.w = 1;

    return(vsum);
}
void VECTOR4D_Sub(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vdiff)
{
    vdiff->x = va->x - vb->x;
    vdiff->y = va->y - vb->y;
    vdiff->z = va->z - vb->z;
    vdiff->w = 1;

}
VECTOR4D VECTOR4D_Sub(VECTOR4D_PTR va, VECTOR4D_PTR vb)
{
    VECTOR4D vdiff;

    vdiff.x = va->x - vb->x;
    vdiff.y = va->y - vb->y;
    vdiff.z = va->z - vb->z;
    vdiff.w = 1;

    return(vdiff);
}
void VECTOR4D_Scale(float k, VECTOR4D_PTR va)
{
    va->x*=k;
    va->y*=k;
    va->z*=k;
    va->w = 1;

}
void VECTOR4D_Scale(float k, VECTOR4D_PTR va, VECTOR4D_PTR vscaled)
{
    vscaled->x = k*va->x;
    vscaled->y = k*va->y;
    vscaled->z = k*va->z;
    vscaled->w = 1;

}
float VECTOR4D_Dot(VECTOR4D_PTR va, VECTOR4D_PTR vb)
{
    return( (va->x * vb->x) + (va->y * vb->y) + (va->z * vb->z) );
}
void VECTOR4D_Cross(VECTOR4D_PTR va,
                    VECTOR4D_PTR vb,
                    VECTOR4D_PTR vn)
{
    vn->x =  ( (va->y * vb->z) - (va->z * vb->y) );
    vn->y = -( (va->x * vb->z) - (va->z * vb->x) );
    vn->z =  ( (va->x * vb->y) - (va->y * vb->x) );
    vn->w = 1;

}
VECTOR4D VECTOR4D_Cross(VECTOR4D_PTR va, VECTOR4D_PTR vb)
{
    VECTOR4D vn;

    vn.x =  ( (va->y * vb->z) - (va->z * vb->y) );
    vn.y = -( (va->x * vb->z) - (va->z * vb->x) );
    vn.z =  ( (va->x * vb->y) - (va->y * vb->x) );
    vn.w = 1;

    return(vn);
}
float VECTOR4D_Length(VECTOR4D_PTR va)
{
    return(sqrtf(va->x*va->x + va->y*va->y + va->z*va->z) );
}
float VECTOR4D_Length_Fast(VECTOR4D_PTR va)
{
    return( Fast_Distance_3D(va->x, va->y, va->z) );
}
void VECTOR4D_Normalize(VECTOR4D_PTR va)
{
    float length = sqrtf(va->x*va->x + va->y*va->y + va->z*va->z);

    if (length < EPSILON_E5)
        return;

    auto length_inv = static_cast<float>(1.0 / length);

    va->x*=length_inv;
    va->y*=length_inv;
    va->z*=length_inv;
    va->w = 1;

}
void VECTOR4D_Normalize(VECTOR4D_PTR va, VECTOR4D_PTR vn)
{
    VECTOR4D_ZERO(vn);

    auto length = static_cast<float>(std::sqrt(va->x * va->x + va->y * va->y + va->z * va->z));

    if (length < EPSILON_E5)
        return;

    auto length_inv = static_cast<float>(1.0 / length);

    vn->x = va->x*length_inv;
    vn->y = va->y*length_inv;
    vn->z = va->z*length_inv;
    vn->w = 1;

}
float VECTOR4D_CosTh(VECTOR4D_PTR va, VECTOR4D_PTR vb)
{
    return(VECTOR4D_Dot(va,vb)/(VECTOR4D_Length(va)*VECTOR4D_Length(vb)));

}
////////////////////////////////////////////////////////////

void Mat_Init_2X2(MATRIX2X2_PTR ma,
                  float m00, float m01,
                  float m10, float m11)

{
    ma->M00 = m00; ma->M01 = m01;
    ma->M10 = m10; ma->M11 = m11;

}
float Mat_Det_2X2(MATRIX2X2_PTR m)
{
    return(m->M00*m->M11 - m->M01*m->M10);
}
void Mat_Add_2X2(MATRIX2X2_PTR ma, MATRIX2X2_PTR mb, MATRIX2X2_PTR msum)
{
    msum->M00 = ma->M00+mb->M00;
    msum->M01 = ma->M01+mb->M01;
    msum->M10 = ma->M10+mb->M10;
    msum->M11 = ma->M11+mb->M11;

}
void Mat_Mul_2X2(MATRIX2X2_PTR ma, MATRIX2X2_PTR mb, MATRIX2X2_PTR mprod)
{
    mprod->M00 = ma->M00*mb->M00 + ma->M01*mb->M10;
    mprod->M01 = ma->M00*mb->M01 + ma->M01*mb->M11;

    mprod->M10 = ma->M10*mb->M00 + ma->M11*mb->M10;
    mprod->M11 = ma->M10*mb->M01 + ma->M11*mb->M11;

}
int Mat_Inverse_2X2(MATRIX2X2_PTR m, MATRIX2X2_PTR mi)
{
    float det = (m->M00*m->M11 - m->M01*m->M10);

    if (std::fabs(det) < EPSILON_E5)
        return(0);

    auto det_inv = static_cast<float>(1.0 / det);

    mi->M00 =  m->M11*det_inv;
    mi->M01 = -m->M01*det_inv;
    mi->M10 = -m->M10*det_inv;
    mi->M11 =  m->M00*det_inv;

    return(1);
}
int Solve_2X2_System(MATRIX2X2_PTR A, MATRIX1X2_PTR X, MATRIX1X2_PTR B)
{
    float det_A = Mat_Det_2X2(A);

    if (std::fabs(det_A) < EPSILON_E5)
        return(0);

    MATRIX2X2 work_mat;

    MAT_COPY_2X2(A, &work_mat);

    MAT_COLUMN_SWAP_2X2(&work_mat, 0, B);

    float det_ABx = Mat_Det_2X2(&work_mat);

    X->M00 = det_ABx/det_A;

    MAT_COPY_2X2(A, &work_mat);

    MAT_COLUMN_SWAP_2X2(&work_mat, 1, B);

    float det_ABy = Mat_Det_2X2(&work_mat);

    X->M01 = det_ABy/det_A;

    return(1);

}
/////////////////////////////////////////////
int Mat_Mul_1X2_3X2(MATRIX1X2_PTR ma, MATRIX3X2_PTR mb, MATRIX1X2_PTR mprod)
{
    for (int col=0; col<2; col++)
    {
        float sum = 0;
        int index;
        for (index=0; index<2; index++)
        {
            sum+=(ma->M[index]*mb->M[index][col]);
        }
        sum+= mb->M[index][col];

        mprod->M[col] = sum;
    }
    return(1);
}
int Mat_Mul_1X3_3X3(MATRIX1X3_PTR ma, MATRIX3X3_PTR mb, MATRIX1X3_PTR mprod)
{
    for (int col=0; col<3; col++)
    {
        float sum = 0;

        for (int index=0; index<3; index++)
        {
            sum+=(ma->M[index]*mb->M[index][col]);
        }
        mprod->M[col] = sum;

    }
    return(1);
}
int Mat_Mul_3X3(MATRIX3X3_PTR ma, MATRIX3X3_PTR mb, MATRIX3X3_PTR mprod)
{
    for (int row=0; row<3; row++)
    {
        for (int col=0; col<3; col++)
        {
            float sum = 0;

            for (int index=0; index<3; index++)
            {
                sum+=(ma->M[row][index]*mb->M[index][col]);
            }
            mprod->M[row][col] = sum;

        }
    }
    return(1);
}
int Mat_Init_3X2(MATRIX3X2_PTR ma,
                        float m00, float m01,
                        float m10, float m11,
                        float m20, float m21)
{
    ma->M[0][0] = m00; ma->M[0][1] = m01;
    ma->M[1][0] = m10; ma->M[1][1] = m11;
    ma->M[2][0] = m20; ma->M[2][1] = m21;

    return(1);
}
void Mat_Add_3X3(MATRIX3X3_PTR ma,
                 MATRIX3X3_PTR mb,
                 MATRIX3X3_PTR msum)
{
    for (int row=0; row<3; row++)
    {
        for (int col=0; col<3; col++)
        {
            msum->M[row][col] = ma->M[row][col] + mb->M[row][col];
        }
    }
}
void Mat_Mul_VECTOR3D_3X3(VECTOR3D_PTR  va,
                          MATRIX3X3_PTR mb,
                          VECTOR3D_PTR  vprod)
{
    for (int col=0; col < 3; col++)
    {
        float sum = 0;
        for (int row=0; row<3; row++)
        {
            sum+=(va->M[row]*mb->M[row][col]);
        }
        vprod->M[col] = sum;
    }
}
int Mat_Inverse_3X3(MATRIX3X3_PTR m, MATRIX3X3_PTR mi)
{
    float det = m->M00*(m->M11*m->M22 - m->M21*m->M12) -
                m->M01*(m->M10*m->M22 - m->M20*m->M12) +
                m->M02*(m->M10*m->M21 - m->M20*m->M11);

    if (std::fabs(det) < EPSILON_E5)
        return(0);

    auto det_inv = static_cast<float>(1.0 / det);

    mi->M00 =  det_inv*(m->M11*m->M22 - m->M21*m->M12);
    mi->M10 = -det_inv*(m->M10*m->M22 - m->M20*m->M12);
    mi->M20 =  det_inv*(m->M10*m->M21 - m->M20*m->M11);

    mi->M01 = -det_inv*(m->M01*m->M22 - m->M21*m->M02);
    mi->M11 =  det_inv*(m->M00*m->M22 - m->M20*m->M02);
    mi->M21 = -det_inv*(m->M00*m->M21 - m->M20*m->M01);

    mi->M02 =  det_inv*(m->M01*m->M12 - m->M11*m->M02);
    mi->M12 = -det_inv*(m->M00*m->M12 - m->M10*m->M02);
    mi->M22 =  det_inv*(m->M00*m->M11 - m->M10*m->M01);

    return(1);
}
void Mat_Init_3X3(MATRIX3X3_PTR ma,
                  float m00, float m01, float m02,
                  float m10, float m11, float m12,
                  float m20, float m21, float m22)
{
    ma->M00 = m00; ma->M01 = m01; ma->M02 = m02;
    ma->M10 = m10; ma->M11 = m11; ma->M12 = m12;
    ma->M20 = m20; ma->M21 = m21; ma->M22 = m22;
}
float Mat_Det_3X3(MATRIX3X3_PTR m)
{
    return(m->M00*(m->M11*m->M22 - m->M21*m->M12) -
           m->M01*(m->M10*m->M22 - m->M20*m->M12) +
           m->M02*(m->M10*m->M21 - m->M20*m->M11) );

}
int Solve_3X3_System(MATRIX3X3_PTR A, MATRIX1X3_PTR X, MATRIX1X3_PTR B)
{
    float det_A = Mat_Det_3X3(A);

    if (std::fabs(det_A) < EPSILON_E5)
        return(0);

    MATRIX3X3 work_mat;

    MAT_COPY_3X3(A, &work_mat);

    MAT_COLUMN_SWAP_3X3(&work_mat, 0, B);

    float det_ABx = Mat_Det_3X3(&work_mat);

    X->M00 = det_ABx/det_A;

    MAT_COPY_3X3(A, &work_mat);

    MAT_COLUMN_SWAP_3X3(&work_mat, 1, B);

    float det_ABy = Mat_Det_3X3(&work_mat);

    X->M01 = det_ABy/det_A;

    MAT_COPY_3X3(A, &work_mat);

    MAT_COLUMN_SWAP_3X3(&work_mat, 2, B);

    float det_ABz = Mat_Det_3X3(&work_mat);

    X->M02 = det_ABz/det_A;
    return(1);
}

///////////////////////////////////////////////
void Mat_Add_4X4(MATRIX4X4_PTR ma,
                 MATRIX4X4_PTR mb,
                 MATRIX4X4_PTR msum)
{
    for (int row=0; row<4; row++)
    {
        for (int col=0; col<4; col++)
        {
            msum->M[row][col] = ma->M[row][col] + mb->M[row][col];
        }
    }
}
void Mat_Mul_4X4(MATRIX4X4_PTR ma,
                 MATRIX4X4_PTR mb,
                 MATRIX4X4_PTR mprod)
{
    for (int row=0; row<4; row++)
    {
        for (int col=0; col<4; col++)
        {
            float sum = 0;

            for (int index=0; index<4; index++)
            {
                sum+=(ma->M[row][index]*mb->M[index][col]);
            }
            mprod->M[row][col] = sum;
        }
    }
}
void Mat_Mul_1X4_4X4(MATRIX1X4_PTR ma,
                     MATRIX4X4_PTR mb,
                     MATRIX1X4_PTR mprod)
{
    for (int col=0; col<4; col++)
    {
        float sum = 0;

        for (int row=0; row<4; row++)
        {
            sum+=(ma->M[row] * mb->M[row][col]);
        }
        mprod->M[col] = sum;
    }
}
void Mat_Mul_VECTOR3D_4X4(VECTOR3D_PTR  va,
                          MATRIX4X4_PTR mb,
                          VECTOR3D_PTR  vprod)
{
    for (int col=0; col < 3; col++)
    {
        float sum = 0;
        int row;
        for (row=0; row<3; row++)
        {
            sum+=(va->M[row]*mb->M[row][col]);
        }

        sum+=mb->M[row][col];

        vprod->M[col] = sum;
    }
}
void Mat_Mul_VECTOR3D_4X3(VECTOR3D_PTR  va,
                          MATRIX4X3_PTR mb,
                          VECTOR3D_PTR  vprod)
{
    for (int col=0; col < 3; col++)
    {
        float sum = 0;
        int row;
        for (row=0; row<3; row++)
        {
            sum+=(va->M[row]*mb->M[row][col]);
        }

        sum+=mb->M[row][col];

        vprod->M[col] = sum;

    }
}
void Mat_Mul_VECTOR4D_4X4(VECTOR4D_PTR  va,
                          MATRIX4X4_PTR mb,
                          VECTOR4D_PTR  vprod)
{
    for (int col=0; col < 4; col++)
    {
        float sum = 0;

        for (int row=0; row<4; row++)
        {
            sum+=(va->M[row]*mb->M[row][col]);
        }
        vprod->M[col] = sum;
    }
}
void Mat_Mul_VECTOR4D_4X3(VECTOR4D_PTR  va,
                          MATRIX4X4_PTR mb,
                          VECTOR4D_PTR  vprod)
{
    for (int col=0; col < 3; col++)
    {
        float sum = 0;

        for (int row=0; row<4; row++)
        {
            sum+=(va->M[row]*mb->M[row][col]);
        }
        vprod->M[col] = sum;

    }
    vprod->M[3] = va->M[3];
}
void Mat_Init_4X4(MATRIX4X4_PTR ma,
                  float m00, float m01, float m02, float m03,
                  float m10, float m11, float m12, float m13,
                  float m20, float m21, float m22, float m23,
                  float m30, float m31, float m32, float m33)

{
    ma->M00 = m00; ma->M01 = m01; ma->M02 = m02; ma->M03 = m03;
    ma->M10 = m10; ma->M11 = m11; ma->M12 = m12; ma->M13 = m13;
    ma->M20 = m20; ma->M21 = m21; ma->M22 = m22; ma->M23 = m23;
    ma->M30 = m30; ma->M31 = m31; ma->M32 = m32; ma->M33 = m33;

}
int Mat_Inverse_4X4(MATRIX4X4_PTR m, MATRIX4X4_PTR mi)
{

    float det =  ( m->M00 * ( m->M11 * m->M22 - m->M12 * m->M21 ) -
                   m->M01 * ( m->M10 * m->M22 - m->M12 * m->M20 ) +
                   m->M02 * ( m->M10 * m->M21 - m->M11 * m->M20 ) );

    if (std::fabs(det) < EPSILON_E5)
        return(0);

    float det_inv  = 1.0f / det;

    mi->M00 =  det_inv * ( m->M11 * m->M22 - m->M12 * m->M21 );
    mi->M01 = -det_inv * ( m->M01 * m->M22 - m->M02 * m->M21 );
    mi->M02 =  det_inv * ( m->M01 * m->M12 - m->M02 * m->M11 );
    mi->M03 = 0.0f;

    mi->M10 = -det_inv * ( m->M10 * m->M22 - m->M12 * m->M20 );
    mi->M11 =  det_inv * ( m->M00 * m->M22 - m->M02 * m->M20 );
    mi->M12 = -det_inv * ( m->M00 * m->M12 - m->M02 * m->M10 );
    mi->M13 = 0.0f;

    mi->M20 =  det_inv * ( m->M10 * m->M21 - m->M11 * m->M20 );
    mi->M21 = -det_inv * ( m->M00 * m->M21 - m->M01 * m->M20 );
    mi->M22 =  det_inv * ( m->M00 * m->M11 - m->M01 * m->M10 );
    mi->M23 = 0.0f;

    mi->M30 = -( m->M30 * mi->M00 + m->M31 * mi->M10 + m->M32 * mi->M20 );
    mi->M31 = -( m->M30 * mi->M01 + m->M31 * mi->M11 + m->M32 * mi->M21 );
    mi->M32 = -( m->M30 * mi->M02 + m->M31 * mi->M12 + m->M32 * mi->M22 );
    mi->M33 = 1.0f;

    return(1);
}
////////////////////////////////////////////////
void QUAT_Add(QUAT_PTR q1, QUAT_PTR q2, QUAT_PTR qsum)
{
    qsum->x = q1->x + q2->x;
    qsum->y = q1->y + q2->y;
    qsum->z = q1->z + q2->z;
    qsum->w = q1->w + q2->w;
}
void QUAT_Sub(QUAT_PTR q1, QUAT_PTR q2, QUAT_PTR qdiff)
{
    qdiff->x = q1->x - q2->x;
    qdiff->y = q1->y - q2->y;
    qdiff->z = q1->z - q2->z;
    qdiff->w = q1->w - q2->w;
}
void QUAT_Conjugate(QUAT_PTR q, QUAT_PTR qconj)
{
    qconj->x = -q->x;
    qconj->y = -q->y;
    qconj->z = -q->z;
    qconj->w = q->w;
}
void QUAT_Scale(QUAT_PTR q, float scale, QUAT_PTR qs)
{
    qs->x = scale*q->x;
    qs->y = scale*q->y;
    qs->z = scale*q->z;
    qs->w = scale*q->w;
}
void QUAT_Scale(QUAT_PTR q, float scale)
{
    q->x*=scale;
    q->y*=scale;
    q->z*=scale;
    q->w*=scale;
}
float QUAT_Norm(QUAT_PTR q)
{
    return(sqrtf(q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z));
}
float QUAT_Norm2(QUAT_PTR q)
{
    return(q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z);
}
void QUAT_Normalize(QUAT_PTR q, QUAT_PTR qn)
{
    auto qlength_inv = static_cast<float>(1.0 / (sqrtf(q->w * q->w + q->x * q->x + q->y * q->y + q->z * q->z)));

    qn->w=q->w*qlength_inv;
    qn->x=q->x*qlength_inv;
    qn->y=q->y*qlength_inv;
    qn->z=q->z*qlength_inv;
}
void QUAT_Normalize(QUAT_PTR q)
{
    auto qlength_inv = static_cast<float>(1.0 / (sqrtf(q->w * q->w + q->x * q->x + q->y * q->y + q->z * q->z)));

    q->w*=qlength_inv;
    q->x*=qlength_inv;
    q->y*=qlength_inv;
    q->z*=qlength_inv;
}
void QUAT_Unit_Inverse(QUAT_PTR q, QUAT_PTR qi)
{
    qi->w =  q->w;
    qi->x = -q->x;
    qi->y = -q->y;
    qi->z = -q->z;
}
void QUAT_Unit_Inverse(QUAT_PTR q)
{
    q->x = -q->x;
    q->y = -q->y;
    q->z = -q->z;
}
void QUAT_Inverse(QUAT_PTR q, QUAT_PTR qi)
{
    auto norm2_inv = static_cast<float>(1.0 / (q->w * q->w + q->x * q->x + q->y * q->y + q->z * q->z));

    qi->w =  q->w*norm2_inv;
    qi->x = -q->x*norm2_inv;
    qi->y = -q->y*norm2_inv;
    qi->z = -q->z*norm2_inv;
}
void QUAT_Inverse(QUAT_PTR q)
{
    auto norm2_inv = static_cast<float>(1.0 / (q->w * q->w + q->x * q->x + q->y * q->y + q->z * q->z));

    q->w =  q->w*norm2_inv;
    q->x = -q->x*norm2_inv;
    q->y = -q->y*norm2_inv;
    q->z = -q->z*norm2_inv;
}
void QUAT_Mul(QUAT_PTR q1, QUAT_PTR q2, QUAT_PTR qprod)
{
    float prd_0 = (q1->z - q1->y) * (q2->y - q2->z);
    float prd_1 = (q1->w + q1->x) * (q2->w + q2->x);
    float prd_2 = (q1->w - q1->x) * (q2->y + q2->z);
    float prd_3 = (q1->y + q1->z) * (q2->w - q2->x);
    float prd_4 = (q1->z - q1->x) * (q2->x - q2->y);
    float prd_5 = (q1->z + q1->x) * (q2->x + q2->y);
    float prd_6 = (q1->w + q1->y) * (q2->w - q2->z);
    float prd_7 = (q1->w - q1->y) * (q2->w + q2->z);

    auto prd_8 = prd_5 + prd_6 + prd_7;
    auto prd_9 = static_cast<float>(0.5 * (prd_4 + prd_8));

    qprod->w = prd_0 + prd_9 - prd_5;
    qprod->x = prd_1 + prd_9 - prd_8;
    qprod->y = prd_2 + prd_9 - prd_7;
    qprod->z = prd_3 + prd_9 - prd_6;
}
void QUAT_Triple_Product(QUAT_PTR q1, QUAT_PTR q2, QUAT_PTR q3,
                         QUAT_PTR qprod)
{
    QUAT qtmp;
    QUAT_Mul(q1,q2,&qtmp);
    QUAT_Mul(&qtmp, q3, qprod);
}
void VECTOR3D_Theta_To_QUAT(QUAT_PTR q, VECTOR3D_PTR v, float theta)
{
    auto theta_div_2 = static_cast<float>((0.5) * theta);

    float sinf_theta = sinf(theta_div_2);

    q->x = sinf_theta * v->x;
    q->y = sinf_theta * v->y;
    q->z = sinf_theta * v->z;
    q->w = cosf( theta_div_2 );
}
void VECTOR4D_Theta_To_QUAT(QUAT_PTR q, VECTOR4D_PTR v, float theta)
{
    auto theta_div_2 = static_cast<float>((0.5) * theta);

    float sinf_theta = sinf(theta_div_2);

    q->x = sinf_theta * v->x;
    q->y = sinf_theta * v->y;
    q->z = sinf_theta * v->z;
    q->w = cosf( theta_div_2 );
}
void EulerZYX_To_QUAT(QUAT_PTR q, float theta_z, float theta_y, float theta_x)
{
    auto cos_z_2 = static_cast<float>(0.5 * cosf(theta_z));
    auto cos_y_2 = static_cast<float>(0.5 * cosf(theta_y));
    auto cos_x_2 = static_cast<float>(0.5 * cosf(theta_x));

    auto sin_z_2 = static_cast<float>(0.5 * sinf(theta_z));
    auto sin_y_2 = static_cast<float>(0.5 * sinf(theta_y));
    auto sin_x_2 = static_cast<float>(0.5 * sinf(theta_x));

    q->w = cos_z_2*cos_y_2*cos_x_2 + sin_z_2*sin_y_2*sin_x_2;
    q->x = cos_z_2*cos_y_2*sin_x_2 - sin_z_2*sin_y_2*cos_x_2;
    q->y = cos_z_2*sin_y_2*cos_x_2 + sin_z_2*cos_y_2*sin_x_2;
    q->z = sin_z_2*cos_y_2*cos_x_2 - cos_z_2*sin_y_2*sin_x_2;
}
void QUAT_To_VECTOR3D_Theta(QUAT_PTR q, VECTOR3D_PTR v, float *theta)
{
    *theta = acosf(q->w);

    auto sinf_theta_inv = static_cast<float>(1.0 / sinf(*theta));

    v->x    = q->x*sinf_theta_inv;
    v->y    = q->y*sinf_theta_inv;
    v->z    = q->z*sinf_theta_inv;

    *theta*=2;
}
////////////////////////////////////////////////

void Init_Parm_Line2D(POINT2D_PTR p_init,
                      POINT2D_PTR p_term, PARMLINE2D_PTR p)
{
    VECTOR2D_INIT(&(p->p0), p_init);

    VECTOR2D_INIT(&(p->p1), p_term);

    VECTOR2D_Build(p_init, p_term, &(p->v));
}void Compute_Parm_Line2D(PARMLINE2D_PTR p, float t, POINT2D_PTR pt)
{
    pt->x = p->p0.x + p->v.x*t;
    pt->y = p->p0.y + p->v.y*t;
}
int Intersect_Parm_Lines2D(PARMLINE2D_PTR p1, PARMLINE2D_PTR p2,
                           float *t1, float *t2)
{
    float det_p1p2 = (p1->v.x*p2->v.y - p1->v.y*p2->v.x);
    if (std::fabs(det_p1p2) <= EPSILON_E5)
    {
        return(PARM_LINE_NO_INTERSECT);

    }
    *t1 = (p2->v.x*(p1->p0.y - p2->p0.y) - p2->v.y*(p1->p0.x - p2->p0.x))
          /det_p1p2;

    *t2 = (p1->v.x*(p1->p0.y - p2->p0.y) - p1->v.y*(p1->p0.x - p2->p0.x))
          /det_p1p2;

    if ((*t1>=0) && (*t1<=1) && (*t2>=0) && (*t2<=1))
        return(PARM_LINE_INTERSECT_IN_SEGMENT);
    else
        return(PARM_LINE_INTERSECT_OUT_SEGMENT);
}
int Intersect_Parm_Lines2D(PARMLINE2D_PTR p1, PARMLINE2D_PTR p2, POINT2D_PTR pt)
{
    float t1, t2, det_p1p2 = (p1->v.x*p2->v.y - p1->v.y*p2->v.x);

    if (std::fabs(det_p1p2) <= EPSILON_E5)
    {
        return(PARM_LINE_NO_INTERSECT);

    }
    t1 = (p2->v.x*(p1->p0.y - p2->p0.y) - p2->v.y*(p1->p0.x - p2->p0.x))
         /det_p1p2;

    t2 = (p1->v.x*(p1->p0.y - p2->p0.y) - p1->v.y*(p1->p0.x - p2->p0.x))
         /det_p1p2;

    pt->x = p1->p0.x + p1->v.x*t1;
    pt->y = p1->p0.y + p1->v.y*t1;

    if ((t1>=0) && (t1<=1) && (t2>=0) && (t2<=1))
        return(PARM_LINE_INTERSECT_IN_SEGMENT);
    else
        return(PARM_LINE_INTERSECT_OUT_SEGMENT);
}
void Init_Parm_Line3D(POINT3D_PTR p_init,
                      POINT3D_PTR p_term, PARMLINE3D_PTR p)
{
    VECTOR3D_INIT(&(p->p0), p_init);

    VECTOR3D_INIT(&(p->p1),p_term);

    VECTOR3D_Build(p_init, p_term, &(p->v));
}
void Compute_Parm_Line3D(PARMLINE3D_PTR p, float t, POINT3D_PTR pt)
{
    pt->x = p->p0.x + p->v.x*t;
    pt->y = p->p0.y + p->v.y*t;
    pt->z = p->p0.z + p->v.z*t;
}
void PLANE3D_Init(PLANE3D_PTR plane, POINT3D_PTR p0,
                  VECTOR3D_PTR normal, int normalize=0)
{
    POINT3D_COPY(&plane->p0, p0);

    if (!normalize)
        VECTOR3D_COPY(&plane->n, normal);
    else
    {
        VECTOR3D_Normalize(normal,&plane->n);
    }
}
float Compute_Point_In_Plane3D(POINT3D_PTR pt, PLANE3D_PTR plane)
{
    float hs = plane->n.x*(pt->x - plane->p0.x) +
               plane->n.y*(pt->y - plane->p0.y) +
               plane->n.z*(pt->z - plane->p0.z);

    return(hs);
}
int Intersect_Parm_Line3D_Plane3D(PARMLINE3D_PTR pline,
                                  PLANE3D_PTR plane,
                                  float *t, POINT3D_PTR pt)
{
    float plane_dot_line = VECTOR3D_Dot(&pline->v, &plane->n);

    if (std::fabs(plane_dot_line) <= EPSILON_E5)
    {
        if (std::fabs(Compute_Point_In_Plane3D(&pline->p0, plane)) <= EPSILON_E5)
            return(PARM_LINE_INTERSECT_EVERYWHERE);
        else
            return(PARM_LINE_NO_INTERSECT);
    }
    *t = -(plane->n.x*pline->p0.x +
           plane->n.y*pline->p0.y +
           plane->n.z*pline->p0.z -
           plane->n.x*plane->p0.x -
           plane->n.y*plane->p0.y -
           plane->n.z*plane->p0.z) / (plane_dot_line);

    pt->x = pline->p0.x + pline->v.x*(*t);
    pt->y = pline->p0.y + pline->v.y*(*t);
    pt->z = pline->p0.z + pline->v.z*(*t);

    if (*t>=0.0 && *t<=1.0)
        return(PARM_LINE_INTERSECT_IN_SEGMENT );
    else
        return(PARM_LINE_INTERSECT_OUT_SEGMENT);
}

/////////////////////////////////////////////////

FIXP16 FIXP16_MUL(FIXP16 fp1, FIXP16 fp2)
{
    FIXP16 fp_prod;
    __asm__ {
            mov eax, fp1
            imul fp2
            shrd eax, edx, 16
    }
}
FIXP16 FIXP16_DIV(FIXP16 fp1, FIXP16 fp2)
{
    __asm__ {
            mov eax, fp1
            cdq
            shld edx, eax, 16
            sal eax, 16
            idiv fp2
    }
}

//////////////////////////////////////////////
void Build_Sin_Cos_Tables()
{
    for (int ang = 0; ang <= 360; ang++)
    {
        float theta = (float)ang*PI/(float)180;
        cos_look[ang] = static_cast<float>(std::cos(theta));
        sin_look[ang] = static_cast<float>(std::sin(theta));
    }
}
