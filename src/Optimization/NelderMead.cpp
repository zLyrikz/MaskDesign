#include "NelderMead.h"
#include <ctime>

#include <iostream>
using std::cout;
using std::endl;

NelderMead::NelderMead()
{
}

NelderMead::~NelderMead()
{
}

void NelderMead::Optimize(std::function<double(double[])> fn, int n, double start[], double xmin[],
    double* ynewlo, double reqmin, double step[], int konvge, int kcount,
    int* icount, int* numres)
{
    clock_t time_start, time_end;
    time_start = clock();

    double ccoeff = 0.5;// contraction coefficient
    double del;
    double dn;
    double dnn;
    double ecoeff = 2.0;// extension coefficient
    double eps = 0.001;
    int i;
    int ihi;// index of max y
    int ilo;// index of min y
    int j;
    int jcount;
    int l;// rank of ystar in all y
    int nn;
    double* p;// the simplex
    double* p2star;// extension or contraction point
    double* pbar;// the centroid of the simplex vertices excepting the vertex with Y value YNEWLO.
    double* pstar;// reflected point
    double rcoeff = 1.0;// reflection coefficient
    double rq;
    double x;
    double* y;
    double y2star;
    double ylo;
    double ystar;
    double z;
    //
    //  Check the input parameters.
    //
    if (reqmin <= 0.0)
    {
        ifault_ = 1;
        return;
    }

    if (n < 1)
    {
        ifault_ = 1;
        return;
    }

    if (konvge < 1)
    {
        ifault_ = 1;
        return;
    }

    p = new double[n * (n + 1)];// the simplex
    pstar = new double[n];// reflected point
    p2star = new double[n];// extension or contraction point
    pbar = new double[n];// the centroid of the simplex vertices excepting the vertex with Y value YNEWLO.
    y = new double[n + 1];// the function value

    *icount = 0;
    *numres = 0;

    jcount = konvge;
    dn = (double)(n);
    nn = n + 1;
    dnn = (double)(nn);
    del = 1.0;
    rq = reqmin * dnn;
    //
    //  Initial or restarted loop.
    //
    for (; ; )
    {
        // put initial point at the end of the simplex
        for (i = 0; i < n; i++)
        {
            p[i + n * n] = start[i];
        }
        y[n] = fn(start);
        *icount = *icount + 1;

        for (j = 0; j < n; j++) 
        {
            // construct j-th point in simplex

            x = start[j];
            start[j] = start[j] + step[j] * del;
            for (i = 0; i < n; i++)
            {
                p[i + j * n] = start[i];
            }
            y[j] = fn(start);
            *icount = *icount + 1;
            start[j] = x;// keep start unchanged
        }
        //                    
        //  The simplex construction is complete.
        //                    
        //  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
        //  the vertex of the simplex to be replaced.
        //                
        ylo = y[0];// min y
        ilo = 0;// index of min y

        for (i = 1; i < nn; i++)
        {
            if (y[i] < ylo)
            {
                ylo = y[i];
                ilo = i;
            }
        }
        //
        //  Inner loop.
        //
        for (; ; )
        {
            if (kcount <= *icount)
            {
                break;
            }
            *ynewlo = y[0];// max y
            ihi = 0;// index of max y

            for (i = 1; i < nn; i++)
            {
                if (*ynewlo < y[i])
                {
                    *ynewlo = y[i];
                    ihi = i;
                }
            }
            //
            //  Calculate PBAR, the centroid of the simplex vertices
            //  excepting the vertex with Y value YNEWLO.
            //
            for (i = 0; i < n; i++)
            {
                // i-th element in the x dimension

                z = 0.0;
                for (j = 0; j < nn; j++)
                {
                    // j-th point in simplex

                    z = z + p[i + j * n];
                }
                z = z - p[i + ihi * n];
                pbar[i] = z / dn;
            }
            //
            //  Reflection through the centroid.
            //
            for (i = 0; i < n; i++)
            {
                pstar[i] = pbar[i] + rcoeff * (pbar[i] - p[i + ihi * n]);
            }
            ystar = fn(pstar);
            *icount = *icount + 1;
            //
            //  Successful reflection, so extension.
            //
            if (ystar < ylo)
            {
                for (i = 0; i < n; i++)
                {
                    p2star[i] = pbar[i] + ecoeff * (pstar[i] - pbar[i]);
                }
                y2star = fn(p2star);
                *icount = *icount + 1;
                //
                //  Check extension.
                //
                if (ystar < y2star) 
                {
                    for (i = 0; i < n; i++)
                    {
                        p[i + ihi * n] = pstar[i];
                    }
                    y[ihi] = ystar;
                }
                //
                //  Retain extension or contraction.
                //
                else
                {
                    for (i = 0; i < n; i++)
                    {
                        p[i + ihi * n] = p2star[i];
                    }
                    y[ihi] = y2star;
                }
            }
            //
            //  No extension.
            //
            else
            {
                l = 0;// decresing order rank of ystar in all y; (0 represent first)
                for (i = 0; i < nn; i++)
                {
                    if (ystar < y[i])
                    {
                        l = l + 1;
                    }
                }

                if (1 < l) // ystar better(smaller) than somebody
                {
                    for (i = 0; i < n; i++)
                    {
                        p[i + ihi * n] = pstar[i];
                    }
                    y[ihi] = ystar;
                }
                //
                //  Contraction on the Y(IHI) side of the centroid.
                //
                else if (l == 0)// ystar still the worst
                {
                    for (i = 0; i < n; i++)
                    {
                        p2star[i] = pbar[i] + ccoeff * (p[i + ihi * n] - pbar[i]);
                    }
                    y2star = fn(p2star);
                    *icount = *icount + 1;
                    //
                    //  Contract the whole simplex.
                    //
                    if (y[ihi] < y2star)
                    {
                        for (j = 0; j < nn; j++)
                        {
                            for (i = 0; i < n; i++)
                            {
                                p[i + j * n] = (p[i + j * n] + p[i + ilo * n]) * 0.5;
                                xmin[i] = p[i + j * n];
                            }
                            y[j] = fn(xmin);
                            *icount = *icount + 1;
                        }
                        ylo = y[0];
                        ilo = 0;

                        for (i = 1; i < nn; i++)
                        {
                            if (y[i] < ylo)
                            {
                                ylo = y[i];
                                ilo = i;
                            }
                        }
                        continue;
                    }
                    //
                    //  Retain contraction.
                    //
                    else
                    {
                        for (i = 0; i < n; i++)
                        {
                            p[i + ihi * n] = p2star[i];
                        }
                        y[ihi] = y2star;
                    }
                }
                //
                //  Contraction on the reflection side of the centroid.
                //
                else if (l == 1)
                {
                    for (i = 0; i < n; i++)
                    {
                        p2star[i] = pbar[i] + ccoeff * (pstar[i] - pbar[i]);
                    }
                    y2star = fn(p2star);
                    *icount = *icount + 1;
                    //
                    //  Retain reflection?
                    //
                    if (y2star <= ystar)
                    {
                        for (i = 0; i < n; i++)
                        {
                            p[i + ihi * n] = p2star[i];
                        }
                        y[ihi] = y2star;
                    }
                    else
                    {
                        for (i = 0; i < n; i++)
                        {
                            p[i + ihi * n] = pstar[i];
                        }
                        y[ihi] = ystar;
                    }
                }
            }
            //
            //  Check if YLO improved.
            //
            if (y[ihi] < ylo)
            {
                ylo = y[ihi];
                ilo = ihi;
            }
            jcount = jcount - 1;

            if (0 < jcount)
            {
                continue;
            }
            //
            //  Check to see if minimum reached.
            //
            if (*icount <= kcount)
            {
                jcount = konvge;

                // average of all y
                z = 0.0;
                for (i = 0; i < nn; i++)
                {
                    z = z + y[i];
                }
                x = z / dnn;

                // variance of all y
                z = 0.0;
                for (i = 0; i < nn; i++)
                {
                    z = z + pow(y[i] - x, 2);
                }

                if (z <= rq)
                {
                    //cout << "variance =" << z/dnn << ", required = " << reqmin << endl;
                    break;// break the inner loop
                }
            }
        }
        //
        //  Factorial tests to check that YNEWLO is a local minimum.
        //
        for (i = 0; i < n; i++)
        {
            xmin[i] = p[i + ilo * n];
        }
        *ynewlo = y[ilo];

        if (kcount < *icount)
        {
            ifault_ = 2;
            break;// break the whole process
        }

        ifault_ = 0;

        //for (i = 0; i < n; i++)
        //{
        //    del = step[i] * eps;
        //    xmin[i] = xmin[i] + del;
        //    z = fn(xmin);
        //    *icount = *icount + 1;
        //    if (z < *ynewlo)
        //    {
        //        ifault_ = 2;// ready for a restart
        //        break;
        //    }
        //    xmin[i] = xmin[i] - del - del;
        //    z = fn(xmin);
        //    *icount = *icount + 1;
        //    if (z < *ynewlo)
        //    {
        //        ifault_ = 2;// ready for a restart
        //        break;
        //    }
        //    xmin[i] = xmin[i] + del;
        //}

        if (ifault_ == 0)
        {
            break;// break the whole process
        }
        ////
        ////  Restart the procedure.
        ////
        //for (i = 0; i < n; i++)
        //{
        //    start[i] = xmin[i];
        //}
        //del = eps;
        //*numres = *numres + 1;
    }
    delete[] p;
    delete[] pstar;
    delete[] p2star;
    delete[] pbar;
    delete[] y;

    time_end = clock();
    //cout << " NelderMead Optimize time = " << double(time_end - time_start) / CLOCKS_PER_SEC << "s" << endl;

    //cout << "\n";
    //cout << "  optimal F(X*) = " << *ynewlo << "\n";

    //cout << "\n";
    //cout << "  Number of iterations = " << *icount << "\n";
    //cout << "  Number of restarts =   " << *numres << "\n";

    return;
}
void NelderMead::PrintFaultInfo() const
{
//    Output, int *_ifault, error indicator.
//    0, no errors detected.
//    1, _min_f_variance, N, or _convergence_check has an illegal value.
//    2, iteration terminated because _max_f_evaluation was exceeded without convergence.
    cout << "trajectory curve adjustment finished: ";
    if (ifault_ == 0)
    {
        cout << "no errors detected." << endl;
    }
    else if (ifault_ == 1)
    {
        cout << "[WARNING NelderMead::PrintFaultInfo] _min_f_variance, N, or _convergence_check has an illegal value." << endl;
    }
    else if (ifault_ == 2)
    {
        cout << "iteration terminated because _max_f_evaluation was exceeded without convergence." << endl;
    }
}
//****************************************************************************80

