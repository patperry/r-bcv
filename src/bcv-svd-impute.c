
#include <assert.h>
#include <string.h>
#include "bcv-vector-private.h"
#include "bcv-matrix-private.h"
#include "bcv-svd-impute.h"

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
};


bcv_svd_impute_t *
bcv_svd_impute_alloc (bcv_index_t m, bcv_index_t n)
{
    bcv_svd_impute_t *result = NULL;
    bcv_svd_impute_t *impute = NULL;
    bcv_index_t mn = BCV_MIN (m,n);
    bcv_matrix_svdjob_t jobu, jobvt;
    bcv_index_t svd_work_size, colmean_work_size, work_size;
    
    assert (m >= 0);
    assert (n >= 0);
    
    jobu  = BCV_MATRIX_SVDJOB_SOME;
    jobvt = BCV_MATRIX_SVDJOB_SOME;
    
    colmean_work_size = n * sizeof (bcv_index_t);
    svd_work_size = (_bcv_lapack_dgesvd_work_len (jobu, jobvt, m, n) 
                     * sizeof (double));
    work_size = BCV_MAX (colmean_work_size, svd_work_size);
    
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
                     const bcv_index_t *indices, bcv_index_t num_indices)
{
    bcv_svd_col_mean_impute (impute, x, indices, num_indices);
    return 0;
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
