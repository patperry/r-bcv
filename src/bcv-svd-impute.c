
#include <assert.h>
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
bcv_svd_col_mean_impute (bcv_svd_impute_t *impute, const bcv_matrix_t *x, 
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
bcv_svd_impute_step (bcv_svd_impute_t *impute, const bcv_matrix_t *x,
                     const bcv_index_t *indices, bcv_index_t num_indices,
                     bcv_index_t k);

/*
 * Decompose xhat = u d vt.  Set ud, d, and vt.  Destroy the memory in xhat.
 */
static bcv_error_t
bcv_svd_impute_decompose_xhat (bcv_svd_impute_t *impute);

/*
 * Reconstruct xhat := u[k] d[k] vt[k] (the first k terms of the SVD).
 * This assumes that ud and vt are stored in @impute.
 */
static void
bcv_svd_impute_reconstruct_xhat (bcv_svd_impute_t *impute, bcv_index_t k);

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
    bcv_matrix_t *xhat;
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
    
    colmean_work_size = n * sizeof (bcv_index_t);
    svd_work_size     = bcv_svd_impute_svd_lwork (m, n) * sizeof (double);
    work_size         = BCV_MAX (colmean_work_size, svd_work_size);
    
    if (colmean_work_size > 0 && svd_work_size > 0
        && (impute = calloc (1, sizeof (bcv_svd_impute_t)))
        && (impute->xhat = calloc (1, sizeof (bcv_matrix_t)))
        && (impute->ud   = calloc (1, sizeof (bcv_matrix_t)))
        && (impute->vt   = calloc (1, sizeof (bcv_matrix_t)))
        && (impute->work = malloc (work_size)))
    {    
        if (mn > 0)
        {
            if ((impute->xhat->data = malloc (m * n * sizeof (double)))
                && (impute->ud->data = malloc (m * mn * sizeof (double)))
                && (impute->vt->data = malloc (n * mn * sizeof (double)))
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
        if (impute->xhat)
        {
            if (impute->xhat->data) free (impute->xhat->data);
            free (impute->xhat);
        }
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
bcv_svd_impute_init (bcv_svd_impute_t *impute, const bcv_matrix_t *x, 
                     const bcv_index_t *indices, bcv_index_t num_indices,
                     bcv_index_t k, double tol, bcv_index_t max_iter)
{
    bcv_error_t err = 0;
    bcv_index_t m, n, lwork;
    bcv_index_t iter = 0;
    double rss0, rss1 = 1.0/0.0, delta;
    
    assert (impute);
    _bcv_assert_valid_matrix (x);
    
    m              = x->m;
    n              = x->n;
    lwork          = bcv_svd_impute_svd_lwork (m, n);
    impute->lwork  = lwork;
    
    bcv_svd_col_mean_impute (impute, x, indices, num_indices);
    do
    {
        rss0 = rss1;
        iter++;
        
        err   = bcv_svd_impute_step (impute, x, indices, num_indices, k);
        rss1  = impute->rss;
        delta = abs (rss1 - rss0) / (2.220446e-16 + rss1);
    }
    while (!err && delta > tol && iter < max_iter);
    
    impute->iter = iter;
    
    return err;
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
bcv_svd_col_mean_impute (bcv_svd_impute_t *impute, const bcv_matrix_t *x, 
                         const bcv_index_t *indices, bcv_index_t num_indices)
{
    bcv_index_t m, n, i, j, idx;
    const bcv_index_t *indices_start = indices;
    const bcv_index_t *indices_end   = indices_start + num_indices;
    const bcv_index_t *ptr;
    
    assert (impute);
    assert (num_indices >= 0);
    
    _bcv_matrix_copy (impute->xhat, x);
    _bcv_matrix_set_indices (impute->xhat, 0.0, indices, num_indices);
    
    m = x->m;
    n = x->n;
    
    if (m > 0 && n > 0)
    {
        bcv_vector_t one     = { m, impute->ud->data, 1 };
        bcv_vector_t mu      = { n, impute->vt->data, 1 };
        bcv_index_t *missing = impute->work;
        bcv_index_t count;

        _bcv_vector_set_constant (&one, 1.0);
        _bcv_blas_dgemv (BCV_MATRIX_TRANS, 1.0, impute->xhat, &one, 0.0, &mu);
        bcv_matrix_miss_counts (m, n, indices, num_indices, NULL, missing);
        
        for (j = 0; j < n; j++)
        {
            count = m - missing[j];
            
            if (count > 0)
            {
                mu.data[j] = mu.data[j] / count;
            }
            else
            {
                mu.data[j] = 0.0;
            }
        }
        
        for (ptr = indices_start; ptr < indices_end; ptr++)
        {
            idx = *ptr;
            i   = idx % m;
            j   = idx / m;
            
            assert (0 <= idx && idx < m*n);
            
            impute->xhat->data[idx] = mu.data[j];
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
bcv_svd_impute_step (bcv_svd_impute_t *impute, const bcv_matrix_t *x,
                     const bcv_index_t *indices, bcv_index_t num_indices,
                     bcv_index_t k)
{
    bcv_error_t err;
    double rss;
    
    err = bcv_svd_impute_decompose_xhat (impute);
    
    if (!err)
    {
        bcv_svd_impute_reconstruct_xhat (impute, k);
        rss = bcv_impute_replace_nonmissing (impute->xhat, x, indices,
                                             num_indices);
        impute->rss = rss;
    }
    
    return err;
}


bcv_error_t
bcv_svd_impute_decompose_xhat (bcv_svd_impute_t *impute)
{
    bcv_index_t i, m, n, mn, ldu;
    bcv_error_t info;
    bcv_vector_t u_i;
    
    assert (impute);

    info = _bcv_lapack_dgesvd (BCV_MATRIX_SVDJOB_SOME, BCV_MATRIX_SVDJOB_SOME,
                               impute->xhat, impute->d, impute->ud,
                               impute->vt, impute->work, impute->lwork);

    m   = impute->xhat->m;
    n   = impute->xhat->n;
    mn  = BCV_MIN (m,n);
    ldu = impute->ud->lda;
    
    u_i.n    = m;
    u_i.data = impute->ud->data;
    u_i.inc  = 1;
    
    for (i = 0; i < mn; i++, u_i.data += ldu)
    {
        _bcv_blas_dscal (impute->d[i], &u_i);
    }
    
    return info;
}


void
bcv_svd_impute_reconstruct_xhat (bcv_svd_impute_t *impute, bcv_index_t k)
{
    assert (impute);
    assert (impute->xhat);
    assert (0 <= k && k < BCV_MIN (impute->xhat->m, impute->xhat->n));
    
    if (k > 0)
    {
        bcv_matrix_t ud_k = *(impute->ud);
        bcv_matrix_t vt_k = *(impute->vt);
    
        ud_k.n = k;
        vt_k.n = k;
    
        _bcv_blas_dgemm (BCV_MATRIX_NOTRANS, BCV_MATRIX_TRANS,  
                         1.0, &ud_k, &vt_k, 0.0, impute->xhat);
    }
    else
    {
        _bcv_matrix_set_constant (impute->xhat, 0.0);
    }
}

double
bcv_impute_replace_nonmissing (bcv_matrix_t *xhat,
                               const bcv_matrix_t *x,
                               const bcv_index_t *indices, 
                               bcv_index_t num_indices)
{
    double rss = 0.0;
    bcv_index_t m, n, ldx, ldxhat, i, i_start, i_end, col_start, col_end;
    bcv_index_t idx, idx_end, next_miss = 0;
    double *x_j, *xhat_j;
    double x_ij, xhat_ij, d_ij, s_ij;
    
    _bcv_assert_valid_matrix (xhat);
    _bcv_assert_valid_matrix (x);
    assert (indices);
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
            assert (0 <= *indices && *indices < m * n);
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

                /* update RSS */
                d_ij = x_ij - xhat_ij;
                s_ij = d_ij * d_ij;
                rss += s_ij;
            }
        }
        
        /* step over the missing element */
        idx++; 
    }
    
    if (m > 0 && n > 0)
    {
        assert (idx == idx_end + 1);
        assert (num_indices == 0);
        assert (col_start == idx_end - m);
        assert (col_end == idx_end);
    }
    
    return rss;
}
