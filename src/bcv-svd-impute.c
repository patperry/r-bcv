
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "bcv-align-private.h"
#include "bcv-vector-private.h"
#include "bcv-matrix-private.h"
#include "bcv-svd-impute.h"


/*
 * Return the length of the work arrays needed for the SVD computation
 * for matrix x of the given dimensions.
 */
static bcv_index_t
bcv_svd_impute_svd_lwork (bcv_index_t m, bcv_index_t n);

static bcv_index_t
bcv_svd_impute_svd_liwork (bcv_index_t m, bcv_index_t n);

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
    bcv_index_t k;

    bcv_matrix_t *ud;
    bcv_matrix_t *vt;
    double *d;
    
    double *work;
    bcv_index_t lwork;
    
    bcv_index_t *iwork;
    bcv_index_t liwork;
};


size_t
bcv_svd_impute_size (bcv_index_t m, bcv_index_t n)
{
    size_t result = 0, total;
    bcv_index_t mn = BCV_MIN (m,n);
    bcv_index_t lwork, liwork;
    
    /* space for the bcv_svd_impute_t and ud, vt */
    total = (sizeof (bcv_svd_impute_t) 
             + (_bcv_alignof (bcv_matrix_t) - 1) 
             + 2 * sizeof (bcv_matrix_t)
             + (_bcv_alignof (double) - 1)
             + (_bcv_alignof (bcv_index_t) - 1));

    if (mn > 0)
    {
        /* space for ud */
        if (m <= SIZE_MAX / mn
            && m * mn <= (SIZE_MAX - total) / sizeof (double))
        {
            total +=  m * mn * sizeof (double);
        
            /* space for vt */
            if (mn <= SIZE_MAX / n
                && mn * n <= (SIZE_MAX - total) / sizeof (double))
            {
                total += mn * n * sizeof (double);
        
                /* space for d */
                if (mn <= (SIZE_MAX - total) / sizeof (double))
                {
                    total += mn * sizeof (double);

                    /* space for work */
                    lwork = bcv_svd_impute_svd_lwork (m, n);
                    if (lwork <= (SIZE_MAX - total) / sizeof (double))
                    {
                        total += lwork * sizeof (double);
                        
                        /* space for iwork */
                        liwork = bcv_svd_impute_svd_liwork (m, n);
                        if (liwork <= (SIZE_MAX - total) / sizeof (bcv_index_t))
                        {
                            total += liwork * sizeof (bcv_index_t);
                            result = total;
                        }
                    }
                }
            }
        }
    }
    /* mn == 0 */
    else
    {
        result = total;
    }
    
    return result;
}


size_t
bcv_svd_impute_align ()
{
    return _bcv_alignof (bcv_svd_impute_t);
}


bcv_svd_impute_t *
bcv_svd_impute_alloc (bcv_index_t m, bcv_index_t n)
{
    bcv_svd_impute_t *result = NULL;
    size_t size = bcv_svd_impute_size (m, n);

    if (size > 0)
        result = calloc (1, size);
    
    return result;
}


void
bcv_svd_impute_free (bcv_svd_impute_t *impute)
{
    if (impute) 
        free (impute);
}


void
bcv_svd_impute_init (bcv_svd_impute_t *impute, 
                     bcv_matrix_t *xhat, const bcv_matrix_t *x, 
                     const bcv_index_t *indices, bcv_index_t num_indices)
{
    bcv_index_t m, n, mn, lwork, liwork;
    size_t bcv_matrix_align = _bcv_alignof (bcv_matrix_t);
    size_t double_align     = _bcv_alignof (double);
    size_t index_align      = _bcv_alignof (bcv_index_t);
    char *mem;
    
    assert (impute);
    _bcv_assert_valid_matrix (x);
    
    m      = x->m;
    n      = x->n;
    mn     = BCV_MIN (m,n);
    lwork  = bcv_svd_impute_svd_lwork (m, n);
    liwork = bcv_svd_impute_svd_liwork (m, n);

    mem = (char *)impute; mem += sizeof (bcv_svd_impute_t);

    mem = BCV_ALIGN_PTR (mem, bcv_matrix_align);
    impute->ud = (void *)mem; mem += sizeof (bcv_matrix_t);
    impute->vt = (void *)mem; mem += sizeof (bcv_matrix_t);
    
    mem = BCV_ALIGN_PTR (mem, double_align);
    impute->ud->data = (void *)mem; mem += m  * mn * sizeof (double);
    impute->vt->data = (void *)mem; mem += mn * n  * sizeof (double);
    impute->d        = (void *)mem; mem += mn * sizeof (double);
    impute->work     = (void *)mem; mem += lwork * sizeof (double);

    mem = BCV_ALIGN_PTR (mem, index_align);
    impute->iwork = (void *)mem; mem += liwork * sizeof (bcv_index_t);
    
    impute->ud->m   = m;
    impute->ud->n   = mn;
    impute->ud->lda = m;
    
    impute->vt->m   = mn;
    impute->vt->n   = n;
    impute->vt->lda = mn;
    
    impute->k      = 0;
    impute->rss    = 0.0;
    impute->lwork  = lwork;
    impute->liwork = liwork;
    
    bcv_svd_col_mean_impute (impute, xhat, x, indices, num_indices);    
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
    bcv_matrix_svdjob_t jobz;
    
    assert (m >= 0);
    assert (n >= 0);
    
    jobz  = BCV_MATRIX_SVDJOB_SOME;
    lwork = _bcv_lapack_dgesdd_work_len (jobz, m, n);
    
    return lwork;
}


bcv_index_t
bcv_svd_impute_svd_liwork (bcv_index_t m, bcv_index_t n)
{
    bcv_index_t result = 0;
    bcv_index_t col_mean_work_len = BCV_MAX (n, 1);
    bcv_index_t svd_iwork_len     = _bcv_lapack_dgesdd_iwork_len (m, n);
    
    if (svd_iwork_len > 0 && col_mean_work_len >= 0)
    {
        result = BCV_MAX (col_mean_work_len, svd_iwork_len);
    }
    
    return result;
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
        bcv_index_t *missing = impute->iwork;
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
                     const bcv_index_t *indices, bcv_index_t num_indices,
                     bcv_index_t k)
{
    bcv_error_t err;
    double rss;
    
    _bcv_assert_valid_matrix (x);
    assert (k >= 0 && k <= BCV_MIN (x->m, x->n));
    
    impute->k = k;
    err       = bcv_svd_impute_decompose_xhat (impute, xhat);
    
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
    bcv_index_t i, m, n, mn, ldu, k;
    bcv_vector_t u_i;
    
    assert (impute);
    _bcv_assert_valid_matrix (xhat);

    m  = xhat->m;
    n  = xhat->n;
    mn = BCV_MIN (m,n);
    k  = impute->k;

    if (mn > 0)
    {
        _bcv_assert_valid_matrix (impute->ud);

        info = _bcv_lapack_dgesdd (BCV_MATRIX_SVDJOB_SOME, 
                                   xhat, impute->d, impute->ud, impute->vt, 
                                   impute->work, impute->lwork,
                                   impute->iwork);

        u_i.n    = m;
        u_i.data = impute->ud->data;
        u_i.inc  = 1;
    
        ldu = impute->ud->lda;
    
        for (i = 0; i < k; i++, u_i.data += ldu)
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
