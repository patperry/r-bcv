
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "bcv-vector-private.h"
#include "bcv-matrix-private.h"
#include "bcv-svd-impute.h"


/*
 * Return the length of the work array needed for the SVD computation
 * for matrix x of the given dimensions.
 */
static bcv_index_t
bcv_svd_impute_svd_lwork (bcv_index_t m, bcv_index_t n);

/*
 * Impute values for @x at the given indices.  Store:
 *   - the completed matrix in impute->xhat,
 *   - the column means in the first column of impute->vt,
 *   - a vector of ones in the first column of impute->ud.
 */
static void
bcv_svd_col_mean_impute (bcv_svd_impute_t *impute,
                         bcv_matrix_t *xhat, const bcv_matrix_t *x, 
                         const bcv_index_t *indices, bcv_index_t num_indices);

/*
 * Count the number of missing indices from each row and column.  The
 * @row_counts and @col_counts arguments are allowed to be NULL.
 */
static void
bcv_matrix_miss_counts (bcv_index_t m, bcv_index_t n,
                        const bcv_index_t *indices, bcv_index_t num_indices,
                        bcv_index_t *row_counts, bcv_index_t *col_counts);

/*
 * Perform one step of the EM algorithm with initial values in impute->xhat.
 * Return error on failure to compute the SVD.
 */
static bcv_error_t
bcv_svd_impute_step (bcv_svd_impute_t *impute, 
                     bcv_matrix_t *xhat, const bcv_matrix_t *x,
                     const bcv_index_t *indices, bcv_index_t num_indices);

/*
 * Decompose xhat = u d vt.  Set ud, d, and vt.  Destroy the memory in xhat.
 */
static bcv_error_t
bcv_svd_impute_decompose_xhat (bcv_svd_impute_t *impute, bcv_matrix_t *xhat);

/*
 * Replace non-missing values in @xhat with the values from @x.  Return the
 * RSS between @xhat and @x at the non-missing values.  This assumes
 * that @indices is sorted in ascending order.
 */
static double
bcv_impute_replace_nonmissing (bcv_matrix_t *xhat,
                               const bcv_matrix_t *x,
                               const bcv_index_t *indices, 
                               bcv_index_t num_indices);

struct _bcv_svd_impute
{
    double rss;
    bcv_index_t iter;

    bcv_index_t k;
    bcv_matrix_t *ud;
    bcv_matrix_t *vt;
    double *d;
    
    void *work;
    bcv_index_t lwork;
};


bcv_svd_impute_t *
bcv_svd_impute_alloc (bcv_index_t m, bcv_index_t n)
{
    bcv_svd_impute_t *result = NULL;
    bcv_svd_impute_t *impute = NULL;
    bcv_index_t mn = BCV_MIN (m,n);
    bcv_index_t svd_work_size, colmean_work_size, work_size;
    
    assert (m >= 0);
    assert (n >= 0);
    
    colmean_work_size = BCV_MAX (n, 1) * sizeof (bcv_index_t);
    svd_work_size     = bcv_svd_impute_svd_lwork (m, n) * sizeof (double);
    work_size         = BCV_MAX (colmean_work_size, svd_work_size);
    
    if (colmean_work_size > 0 && svd_work_size > 0
        && (impute       = calloc (1, sizeof (bcv_svd_impute_t)))
        && (impute->ud   = calloc (1, sizeof (bcv_matrix_t)))
        && (impute->vt   = calloc (1, sizeof (bcv_matrix_t)))
        && (impute->work = malloc (work_size)))
    {    
        if (mn > 0)
        {
            if ((impute->ud->data = malloc (m * mn * sizeof (double)))
                && (impute->vt->data = malloc (mn * n * sizeof (double)))
                && (impute->d        = malloc (mn * sizeof (double))))
            {
                result = impute;
            }
        }
        else
        {
            result = impute;
        }
    }
    
    if (result == NULL)
        bcv_svd_impute_free (impute);
    
    return result;
}


void
bcv_svd_impute_free (bcv_svd_impute_t *impute)
{
    if (impute)
    {
        if (impute->ud)
        {
            if (impute->ud->data) free (impute->ud->data);
            free (impute->ud);
        }
        if (impute->vt)
        {
            if (impute->vt->data) free (impute->vt->data);
            free (impute->vt);
        }
        
        if (impute->d)    free (impute->d);
        if (impute->work) free (impute->work);
        
        free (impute);
    }
}


bcv_error_t
bcv_svd_impute (bcv_svd_impute_t *impute,
                bcv_matrix_t *xhat, const bcv_matrix_t *x, 
                const bcv_index_t *indices, bcv_index_t num_indices,
                bcv_index_t k, double tol, bcv_index_t max_iter)
{
    bcv_error_t err = 0;
    bcv_index_t m, n, mn, lwork;
    bcv_index_t iter = 0;
    double rss0, rss1 = BCV_DBL_POSINF, delta;
    
    assert (impute);
    _bcv_assert_valid_matrix (x);
    
    m     = x->m;
    n     = x->n;
    mn    = BCV_MIN (m,n);
    lwork = bcv_svd_impute_svd_lwork (m, n);

    impute->k       = k;
    impute->ud->m   = m;
    impute->ud->n   = mn;
    impute->ud->lda = m;
    impute->vt->m   = mn;
    impute->vt->n   = n;
    impute->vt->lda = mn;
    impute->lwork   = lwork;
    
    bcv_svd_col_mean_impute (impute, xhat, x, indices, num_indices);
    do
    {
        rss0 = rss1;
        iter++;
        
        err   = bcv_svd_impute_step (impute, xhat, x, indices, num_indices);
        rss1  = impute->rss;
        delta = fabs (rss1 - rss0) / (BCV_DBL_EPSILON + rss1);
        
        printf ("iter: %d rss: %g  delta: %g\n", iter, rss1, delta);
    }
    while (!err && delta > tol && iter < max_iter);
    
    impute->iter = iter;
    
    return err;
}


bcv_index_t
bcv_svd_impute_get_iter (const bcv_svd_impute_t *impute)
{
    assert (impute);
    return impute->iter;
}


double
bcv_svd_impute_get_rss (const bcv_svd_impute_t *impute)
{
    assert (impute);
    return impute->rss;
}


bcv_index_t
bcv_svd_impute_svd_lwork (bcv_index_t m, bcv_index_t n)
{
    bcv_index_t lwork = 0;
    bcv_matrix_svdjob_t jobu, jobvt;
    
    assert (m >= 0);
    assert (n >= 0);
    
    jobu  = BCV_MATRIX_SVDJOB_SOME;
    jobvt = BCV_MATRIX_SVDJOB_SOME;
    lwork = _bcv_lapack_dgesvd_work_len (jobu, jobvt, m, n);
    
    return lwork;
}


void
bcv_svd_col_mean_impute (bcv_svd_impute_t *impute, 
                         bcv_matrix_t *xhat,  const bcv_matrix_t *x, 
                         const bcv_index_t *indices, bcv_index_t num_indices)
{
    bcv_index_t m, n, i, j, idx, incmu;
    const bcv_index_t *indices_start = indices;
    const bcv_index_t *indices_end   = indices_start + num_indices;
    const bcv_index_t *ptr;
    
    assert (impute);
    assert (num_indices >= 0);
    
    _bcv_matrix_copy (xhat, x);
    _bcv_matrix_set_indices (xhat, 0.0, indices, num_indices);
    
    m     = x->m;
    n     = x->n;
    incmu = impute->vt->lda;
    
    if (m > 0 && n > 0)
    {
        bcv_vector_t one     = { m, impute->ud->data, 1 };
        bcv_vector_t mu      = { n, impute->vt->data, incmu };
        bcv_index_t *missing = impute->work;
        bcv_index_t count;

        _bcv_vector_set_constant (&one, 1.0);
        _bcv_blas_dgemv (BCV_MATRIX_TRANS, 1.0, xhat, &one, 0.0, &mu);
        bcv_matrix_miss_counts (m, n, indices, num_indices, NULL, missing);
        
        for (j = 0; j < n; j++)
        {
            count = m - missing[j];
            
            if (count > 0)
            {
                mu.data[ j*incmu ] = mu.data[ j*incmu ] / count;
            }
            else
            {
                mu.data[ j*incmu ] = 0.0;
            }
        }
        
        for (ptr = indices_start; ptr < indices_end; ptr++)
        {
            idx = *ptr;
            i   = idx % m;
            j   = idx / m;
            
            assert (0 <= idx && idx < m*n);
            
            xhat->data[idx] = mu.data[ j*incmu ];
        }
    }
}

void
bcv_matrix_miss_counts (bcv_index_t m, bcv_index_t n,
                        const bcv_index_t *indices, bcv_index_t num_indices,
                        bcv_index_t *row_counts, bcv_index_t *col_counts)
{
    bcv_index_t i, j, idx;
    const bcv_index_t *ptr, *start, *end;
    
    assert (num_indices >= 0);
    assert (indices || num_indices == 0);
    
    if (row_counts)
        memset (row_counts, 0, m * sizeof (bcv_index_t));
    if (col_counts)
        memset (col_counts, 0, n * sizeof (bcv_index_t));
    
    start = indices;
    end   = start + num_indices;
    
    for (ptr = start; ptr < end; ptr++)
    {
        idx = *ptr;
        i   = idx % m;
        j   = idx / m;
        
        assert (0 <= idx && idx < m*n);
        
        if (row_counts)
            row_counts[i] = row_counts[i] + 1;
            
        if (col_counts)
            col_counts[j] = col_counts[j] + 1;
    }
}

bcv_error_t
bcv_svd_impute_step (bcv_svd_impute_t *impute, 
                     bcv_matrix_t *xhat, const bcv_matrix_t *x,
                     const bcv_index_t *indices, bcv_index_t num_indices)
{
    bcv_error_t err;
    double rss;
    
    err = bcv_svd_impute_decompose_xhat (impute, xhat);
    
    if (!err)
    {
        bcv_svd_impute_get_svd (impute, xhat);
        rss = bcv_impute_replace_nonmissing (xhat, x, indices, num_indices);
        impute->rss = rss;
    }
    
    return err;
}


bcv_error_t
bcv_svd_impute_decompose_xhat (bcv_svd_impute_t *impute, bcv_matrix_t *xhat)
{
    bcv_error_t info = 0;
    bcv_index_t i, m, n, mn, ldu;
    bcv_vector_t u_i;
    
    assert (impute);
    _bcv_assert_valid_matrix (xhat);

    m  = xhat->m;
    n  = xhat->n;
    mn = BCV_MIN (m,n);

    if (mn > 0)
    {
        _bcv_assert_valid_matrix (impute->ud);

        info = _bcv_lapack_dgesvd (BCV_MATRIX_SVDJOB_SOME, BCV_MATRIX_SVDJOB_SOME,
                                   xhat, impute->d, impute->ud,
                                   impute->vt, impute->work, impute->lwork);

        u_i.n    = m;
        u_i.data = impute->ud->data;
        u_i.inc  = 1;
    
        ldu = impute->ud->lda;
    
        for (i = 0; i < mn; i++, u_i.data += ldu)
        {
            _bcv_blas_dscal (impute->d[i], &u_i);
        }
    }
    
    return info;
}


void
bcv_svd_impute_get_svd (const bcv_svd_impute_t *impute, bcv_matrix_t *udvt)
{
    bcv_index_t k;
    
    assert (impute);
    assert (udvt);
    
    k = impute->k;
    
    assert (0 <= k && k <= BCV_MIN (udvt->m, udvt->n));
    
    if (k > 0)
    {
        bcv_matrix_t ud_k = *(impute->ud);
        bcv_matrix_t vt_k = *(impute->vt);
    
        ud_k.n = k;
        vt_k.m = k;
    
        _bcv_blas_dgemm (BCV_MATRIX_NOTRANS, BCV_MATRIX_NOTRANS,  
                         1.0, &ud_k, &vt_k, 0.0, udvt);
    }
    else
    {
        _bcv_matrix_set_constant (udvt, 0.0);
    }
}


double
bcv_impute_replace_nonmissing (bcv_matrix_t *xhat,
                               const bcv_matrix_t *x,
                               const bcv_index_t *indices, 
                               bcv_index_t num_indices)
{
    double rss = 1.0, scale = 0.0;
    bcv_index_t m, n, ldx, ldxhat, i, i_start, i_end, col_start, col_end;
    bcv_index_t idx, idx_end, next_miss = 0;
    double *x_j, *xhat_j;
    double x_ij, xhat_ij, d_ij, abs_d_ij;
    
    _bcv_assert_valid_matrix (xhat);
    _bcv_assert_valid_matrix (x);
    assert (indices || num_indices == 0);
    assert (num_indices >= 0);

    m      = x->m;
    n      = x->n;
    ldx    = x->lda;
    ldxhat = xhat->lda;

    assert (xhat->m == m);
    assert (xhat->n == n);

    /* initialize idx to the first index and idx_end to the last */
    idx     = 0;
    idx_end = m*n;

    /* initialize the column variables to refer to the first column */
    x_j       = x->data;
    xhat_j    = xhat->data;
    col_start = 0;
    col_end   = m;
    
    while (idx < idx_end)
    {
        /* get the index of the next missing value */
        if (num_indices == 0)
        {
            next_miss = idx_end;
        }
        else
        {
            assert (0 <= *indices && *indices < idx_end);
            assert (*indices >= next_miss);
            
            next_miss = *indices;

            /* step over duplicates */
            do 
            {
                indices++;
                num_indices--;
            }
            while (num_indices > 0 && *indices == next_miss);
        }
        
        /* copy all of the indices up to the next missing value */
        while (idx < next_miss)
        {
            /* advance to the column of the current index */
            while (idx >= col_end)
            {
                x_j      += ldx;
                xhat_j   += ldxhat;
                col_start = col_end;
                col_end  += m;
            }
            
            i_start = idx - col_start;
            
            /* case1: part of the column is missing */
            if (next_miss < col_end)
            {
                i_end = next_miss - col_start;
                idx   = next_miss;
            }
            /* case2: the entire column is non-missing */
            else
            {
                i_end = m;
                idx   = col_end;
            }
            
            for (i = i_start; i < i_end; i++)
            {
                /* replace xhat[i,j] with x[i,j] */
                x_ij      = x_j[i];
                xhat_ij   = xhat_j[i];
                xhat_j[i] = x_ij;

                /* numerically-stable RSS update based on LAPACK's dnrm2 */
                d_ij = x_ij - xhat_ij;
                if (d_ij != 0)
                {
                    abs_d_ij = fabs (d_ij);
                    if (scale < abs_d_ij)
                    {
                        double ratio = scale / abs_d_ij;
                        rss   = 1.0 + rss * (ratio * ratio);
                        scale = abs_d_ij;
                    }
                    else
                    {
                        double ratio_inv = abs_d_ij / scale;
                        rss   = rss + ratio_inv * ratio_inv;
                    }
                }
            }
        }
        
        /* step over the missing element */
        if (idx != idx_end)
            idx++; 
    }
    
    assert (idx == idx_end);
    assert (num_indices == 0);
    
    return rss * (scale * scale);
}
