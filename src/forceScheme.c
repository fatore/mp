#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#if defined(i386) || defined(i486) || \
    defined(intel) || defined(x86) || defined(i86pc) || \
    defined(__alpha) || defined(__osf__)
#define __LITTLE_ENDIAN
#endif

#ifdef __LITTLE_ENDIAN
#define __HI(x) *(1+(int*)&x)
#define __LO(x) *(int*)&x
#define __HIp(x) *(1+(int*)x)
#define __LOp(x) *(int*)x
#else
#define __HI(x) *(int*)&x
#define __LO(x) *(1+(int*)&x)
#define __HIp(x) *(int*)x
#define __LOp(x) *(1+(int*)x)
#endif

double __ieee754_hypot(double x, double y)
{
    double a=x,b=y,t1,t2,y1,y2,w;
    int j,k,ha,hb;

    ha = __HI(x)&0x7fffffff;    /* high word of  x */
    hb = __HI(y)&0x7fffffff;    /* high word of  y */
    if(hb > ha) {a=y;b=x;j=ha; ha=hb;hb=j;} else {a=x;b=y;}
    __HI(a) = ha;    /* a <- |a| */
    __HI(b) = hb;    /* b <- |b| */
    if((ha-hb)>0x3c00000) {return a+b;} /* x/y > 2**60 */
    k=0;
    if(ha > 0x5f300000) {    /* a>2**500 */
       if(ha >= 0x7ff00000) {    /* Inf or NaN */
           w = a+b;            /* for sNaN */
           if(((ha&0xfffff)|__LO(a))==0) w = a;
           if(((hb^0x7ff00000)|__LO(b))==0) w = b;
           return w;
       }
       /* scale a and b by 2**-600 */
       ha -= 0x25800000; hb -= 0x25800000;    k += 600;
       __HI(a) = ha;
       __HI(b) = hb;
    }
    if(hb < 0x20b00000) {    /* b < 2**-500 */
        if(hb <= 0x000fffff) {    /* subnormal b or 0 */
        if((hb|(__LO(b)))==0) return a;
        t1=0;
        __HI(t1) = 0x7fd00000;    /* t1=2^1022 */
        b *= t1;
        a *= t1;
        k -= 1022;
        } else {        /* scale a and b by 2^600 */
            ha += 0x25800000;     /* a *= 2^600 */
        hb += 0x25800000;    /* b *= 2^600 */
        k -= 600;
           __HI(a) = ha;
           __HI(b) = hb;
        }
    }
    /* medium size a and b */
    w = a-b;
    if (w>b) {
        t1 = 0;
        __HI(t1) = ha;
        t2 = a-t1;
        w  = sqrt(t1*t1-(b*(-b)-t2*(a+t1)));
    } else {
        a  = a+a;
        y1 = 0;
        __HI(y1) = hb;
        y2 = b - y1;
        t1 = 0;
        __HI(t1) = ha+0x00100000;
        t2 = a - t1;
        w  = sqrt(t1*y1-(w*(-w)-(t1*y2+t2*b)));
    }
    if(k!=0) {
        t1 = 1.0;
        __HI(t1) += (k<<20);
        return t1*w;
    } else return w;
}

/*
 * Arrange the N elements of ARRAY in random order.
 * Only effective if N is much smaller than RAND_MAX;
 * if this may not be the case, use a better random
 * number generator.
 */
static void shuffle(size_t *array, size_t n)
{
    if (n > 1) {
        size_t i;
        for (i = 0; i < n - 1; i++) {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
            int t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}

void force_scheme(double p[],
        const double dmat[],
        const int *n,
        const int *max_it,
        const double *tol,
        const double *EPSILON,
        const double *fraction)
{
    double prev_delta_sum = 1. / 0., delta_sum, d2, dn, delta, diff_x, diff_y;
    size_t i, j, k, p1_index, p2_index;
    size_t size = (size_t) *n;

    size_t s_j[size], s_k[size];
    for (i = 0; i < size; i++)
        s_j[i] = s_k[i] = i;

    for (i = 0; i < *max_it; i++) {
        shuffle(s_j, size);
        for (j = 0; j < size; j++) {
            p1_index = s_j[j];

            delta_sum = 0;
            shuffle(s_k, size);
            for (k = 0; k < size; k++) {
                p2_index = s_k[k];

                if (p1_index == p2_index)
                    continue;

                diff_x = p[p2_index] - p[p1_index];
                diff_y = p[p2_index + size] - p[p1_index + size];
                d2 = __ieee754_hypot(diff_x, diff_y);
                if (d2 < *EPSILON)
                    d2 = *EPSILON;

                dn = dmat[p1_index * size + p2_index];
                delta = (dn - d2) / *fraction;
                delta_sum += fabs(delta);
                p[p2_index]        += delta * (diff_x / d2);
                p[p2_index + size] += delta * (diff_y / d2);
            }

        }

        if (fabs(prev_delta_sum - delta_sum) < *tol)
            break;
        prev_delta_sum = delta_sum;
    }
}
