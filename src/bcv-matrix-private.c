
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "bcv-vector-private.h"
#include "bcv-matrix-private.h"

#include <R_ext/Lapack.h>

#define _bcv_dnrm2_f77  dnrm2_
#define _bcv_dscal_f77  dscal_
#define _bcv_dcopy_f77  dcopy_

#define _bcv_dgemv_f77  dgemv_
#define _bcv_dger_f77   dger_

#define _bcv_dgemm_f77  dgemm_

#define _bcv_dlacpy_f77 dlacpy_
#define _bcv_dlange_f77 dlange_
#define _bcv_dgesvd_f77 dgesvd_
#define _bcv_dgesdd_f77 dgesdd_
#define _bcv_dgebrd_f77 dgebrd_
#define _bcv_dormbr_f77 dormbr_
#define _bcv_dbdsqr_f77 dbdsqr_


static char *BCV_BLAS_TRANS_CODES[] = { "N", "T", "C" };
#define _BCV_F77_TRANS(x) \
    BCV_BLAS_TRANS_CODES[(int) (x) - (int) BCV_MATRIX_NOTRANS]

static char *BCV_BLAS_UPLO_CODES[] = { "U", "L", "NA" };
#define _BCV_F77_UPLO(x) \
    BCV_BLAS_UPLO_CODES[(int) (x) - (int) BCV_MATRIX_UPPER]

static char *BCV_BLAS_SIDE_CODES[] = { "L", "R" };
#define _BCV_F77_SIDE(x) \
    BCV_BLAS_SIDE_CODES[(int) (x) - (int) BCV_MATRIX_LEFT]

static char *BCV_LAPACK_NORM_CODES[] = { "M", "O", "I", "F" };
#define _BCV_F77_NORM(x) \
    BCV_LAPACK_NORM_CODES[(int) (x) - (int) BCV_MATRIX_NORM_MAXABS]

static char *BCV_LAPACK_SVDJOB_CODES[] = { "A", "S", "O", "N" };
#define _BCV_F77_SVDJOB(x) \
    BCV_LAPACK_SVDJOB_CODES[(int) (x) - (int) BCV_MATRIX_SVDJOB_ALL ]

static char *BCV_LAPACK_VECT_CODES[] = { "P", "Q" };
#define _BCV_F77_VECT(x) \
    BCV_LAPACK_VECT_CODES[(int) (x) - (int) BCV_MATRIX_VECT_P]

/*
static char *BCV_BLAS_DIAG_CODES[] = { "N", "U" };
#define _BCV_DIAG(x) BCV_BLAS_DIAG_CODES[(int) (x) - (int) BCV_NONUNIT]

static char *BCV_LAPACK_JOB_CODES[] = { "A", "S", "O", "N" };
#define _BCV_JOB(x) BCV_LAPACK_JOB_CODES[(int) (x) - (int) BCV_JOB_ALL]
*/

void
_bcv_matrix_set_constant (bcv_matrix_t *a, double value)
{
    _bcv_assert_valid_matrix (a);

    bcv_index_t m, n, lda;
    
    m   = a->m;
    n   = a->n;
    lda = a ->lda;
    
    if (lda == m)
    {
        bcv_vector_t vec_a = { m*n, a->data, 1 };
        _bcv_vector_set_constant (&vec_a, value);
    }
    else
    {
        bcv_index_t j;
        bcv_vector_t a_j = { m, a->data, 1 };
        
        for (j = 0; j < n; j++, a_j.data += lda)
        {
            _bcv_vector_set_constant (&a_j, value);
        }
    }
}

void
_bcv_matrix_set_identity (bcv_matrix_t *a)
{
    bcv_index_t m, n, mn, lda, j;
    double *data;
    
    _bcv_assert_valid_matrix (a);
    m    = a->m;
    n    = a->n;
    data = a->data;
    lda  = a->lda;
    
    if (m > 0 && n > 0)
    {
        if (lda == m)
        {
            memset (data, 0, m * n * sizeof (double));

            mn = BCV_MIN (m, n);
            for (j = 0; j < mn; j++, data += lda + 1)
                *data = 1.0;
        }
        else
        {
            for (j = 0; j < n; j++, data += lda)
            {
                memset (data, 0, m * sizeof (double));
                if (j < m) data[j] = 1.0;
            }
        }
    }
}


void
_bcv_matrix_set_indices (bcv_matrix_t *a, double value, 
                         const bcv_index_t *indices, bcv_index_t num_indices)
{
    bcv_index_t idx, i, j, m, n, lda;
    double *data;
    const bcv_index_t *indices_end = indices + num_indices;
    const bcv_index_t *pidx;

    _bcv_assert_valid_matrix (a);
    assert (indices || num_indices == 0);
    
    m    = a->m;
    n    = a->n;
    data = a->data;
    lda  = a->lda;
    
    if (lda == m)
    {
        for (pidx = indices; pidx < indices_end; pidx++)
        {
            data[ *pidx ] = value;
        }
    }
    else
    {
        for (pidx = indices; pidx < indices_end; pidx++)
        {
            idx = *pidx;
            i   = idx % m;
            j   = idx / m;
            data[ i + j*lda ] = value;
        }
    }
}


void
_bcv_matrix_copy (bcv_matrix_t *dst, const bcv_matrix_t *src)
{
    _bcv_assert_valid_matrix (dst);
    _bcv_assert_valid_matrix (src);

    _bcv_lapack_dlacpy (BCV_MATRIX_UPLO_NA, src, dst);
}


void
_bcv_matrix_permute_copy (bcv_matrix_t *dst, const bcv_matrix_t *src,
                          bcv_index_t *p, bcv_index_t *q)
{
    bcv_index_t m, n, i, j;
    
    _bcv_assert_valid_matrix (dst);
    _bcv_assert_valid_matrix (src);
    
    m = dst->m;
    n = dst->n;
    
    assert (m == src->m);
    assert (n == src->n);
    
    if (n > 0 && n > 0) 
    {
        if (p == NULL && q == NULL) 
        {
            _bcv_matrix_copy (dst, src);
        } 
        else if (p == NULL) /* permute the columns only */
        { 
            for (j = 0; j < n; j++) 
            {
                bcv_vector_t src_c = { m, src->data + j * src->lda,    1 };
                bcv_vector_t dst_c = { m, dst->data + q[j] * dst->lda, 1 };

                assert (0 <= q[j] && q[j] < n);
                _bcv_blas_dcopy (&src_c, &dst_c);
            }
        }
        else if (q == NULL) /* permute the rows only */
        {
            for (i = 0; i < m; i++)
            {
                bcv_vector_t src_r = { n, src->data + i,    src->lda };
                bcv_vector_t dst_r = { n, dst->data + p[i], dst->lda };
                
                assert (0 <= p[i] && p[i] < m);
                _bcv_blas_dcopy (&src_r, &dst_r);
            }
        }
        else /* permute both rows and columns */
        {
            for (i = 0; i < m; i++) 
                assert (p[i] >= 0 && p[i] < m);

            for (j = 0; j < n; j++)
            {
                double *src_c_data = src->data + j * src->lda;
                double *dst_c_data = dst->data + q[j] * dst->lda;
            
                assert (0 <= q[j] && q[j] < n);
                for (i = 0; i < m; i++)
                    dst_c_data[p[i]] = src_c_data[i];
            }
        }
    }
}


double
_bcv_matrix_norm_frob (const bcv_matrix_t *a)
{
    return _bcv_lapack_dlange (BCV_MATRIX_NORM_FROB, a, NULL);
}


double
_bcv_blas_dnrm2 (const bcv_vector_t *x)
{
    double result;
    
    _bcv_assert_valid_vector (x);
    result = _bcv_dnrm2_f77 (&(x->n), x->data, &(x->inc));
    
    return result;
}


void
_bcv_blas_dscal (double alpha, bcv_vector_t *x)
{
    _bcv_assert_valid_vector (x);
    _bcv_dscal_f77 (&(x->n), &alpha, x->data, &(x->inc));
}


void
_bcv_blas_dcopy (const bcv_vector_t *x, bcv_vector_t *y)
{
    _bcv_assert_valid_vector (x);    
    _bcv_assert_valid_vector (y);    
    assert (x->n == y->n);
    
    _bcv_dcopy_f77 (&(y->n), x->data, &(x->inc), y->data, &(y->inc));
}


void
_bcv_blas_dgemv (bcv_matrix_transpose_t transA,
                 double alpha, const bcv_matrix_t *a, const bcv_vector_t *x,
                 double beta, bcv_vector_t *y)
{
    _bcv_assert_valid_matrix (a);
    _bcv_assert_valid_vector (x);
    _bcv_assert_valid_vector (y);
    
    if (transA == BCV_MATRIX_NOTRANS) {
        assert (a->n == x->n);
        assert (a->m == y->n);
    } else {
        assert (a->m == x->n);
        assert (a->n == y->n);
    }

    _bcv_dgemv_f77 (_BCV_F77_TRANS (transA), &(a->m), &(a->n),
                    &alpha, a->data, &(a->lda), x->data, &(x->inc), 
                    &beta, y->data, &(y->inc));
}


void
_bcv_blas_dger (double alpha, const bcv_vector_t *x, const bcv_vector_t *y, 
                bcv_matrix_t *a)
{
    _bcv_assert_valid_vector (x);
    _bcv_assert_valid_vector (y);    
    _bcv_assert_valid_matrix (a);
    
    assert (a->m == x->n);
    assert (a->n == y->n);
    
    _bcv_dger_f77 (&(a->m), &(a->n), &alpha, x->data, &(x->inc),
                   y->data, &(y->inc), a->data, &(a->lda));
}


void
_bcv_blas_dgemm (bcv_matrix_transpose_t transA,
                 bcv_matrix_transpose_t transB,
                 double alpha, const bcv_matrix_t *a,
                 const bcv_matrix_t *b, double beta,
                 bcv_matrix_t *c)
{
    bcv_index_t m, n, k, k2;
    
    _bcv_assert_valid_matrix (a);
    _bcv_assert_valid_matrix (b);
    _bcv_assert_valid_matrix (c);

    m  = (transA == BCV_MATRIX_NOTRANS) ? a->m : a->n;
    k  = (transA == BCV_MATRIX_NOTRANS) ? a->n : a->m;
    k2 = (transB == BCV_MATRIX_NOTRANS) ? b->m : b->n;
    n  = (transB == BCV_MATRIX_NOTRANS) ? b->n : b->m;
        
    assert (c->m == m);
    assert (c->n == n);
    assert (k == k2);
    
    _bcv_dgemm_f77 (_BCV_F77_TRANS (transA), _BCV_F77_TRANS (transB),
                    &m, &n, &k, &alpha, 
                    a->data, &(a->lda), 
                    b->data, &(b->lda),
                    &beta,
                    c->data, &(c->lda));
}

void
_bcv_lapack_dlacpy (bcv_matrix_uplo_t uplo, const bcv_matrix_t *a, 
                    bcv_matrix_t *b)
{
    _bcv_assert_valid_matrix (a);
    _bcv_assert_valid_matrix (b);
    _bcv_dlacpy_f77 (_BCV_F77_UPLO (uplo), &(b->m), &(b->n), 
                     a->data, &(a->lda), b->data, &(b->lda));
}


double
_bcv_lapack_dlange (bcv_matrix_norm_t norm, const bcv_matrix_t *a, 
                    double *work)
{
    double result = 0.0;
    
    _bcv_assert_valid_matrix (a);
    
    if (a->m > 0 && a->n > 0)
    {
        result = _bcv_dlange_f77 (_BCV_F77_NORM (norm), &(a->m), &(a->n),
                                  a->data, &(a->lda), work);
    }
    
    return result;
}


bcv_index_t
_bcv_lapack_dlange_work_len (bcv_matrix_norm_t norm, 
                             bcv_index_t m, bcv_index_t n)
{
    bcv_index_t result = 0;
    
    assert (m >= 0);
    assert (n >= 0);
        
    if (norm == BCV_MATRIX_NORM_INF)
    {
        result = BCV_MAX (1, m);
    }
    
    return result;
}

bcv_error_t
_bcv_lapack_dgesvd (bcv_matrix_svdjob_t jobu, bcv_matrix_svdjob_t jobvt,
                    bcv_matrix_t *a, double *s, bcv_matrix_t *u,
                    bcv_matrix_t *vt, double *work, bcv_index_t lwork)
{
    bcv_index_t m, n, mn, lda, ldu = 1, ldvt = 1;
    bcv_error_t info = 0;
    double *u_data = NULL, *vt_data = NULL;
    
    _bcv_assert_valid_matrix (a);
    
    m   = a->m;
    n   = a->n;
    mn  = BCV_MIN (m, n);
    lda = a->lda;
    
    assert (!(jobu == BCV_MATRIX_SVDJOB_OVERWRITE
              && jobvt == BCV_MATRIX_SVDJOB_OVERWRITE));
    
    if (jobu == BCV_MATRIX_SVDJOB_ALL || jobu == BCV_MATRIX_SVDJOB_SOME) {
        _bcv_assert_valid_matrix (u);
        
        assert (u->m == m);
        
        if (jobu == BCV_MATRIX_SVDJOB_ALL)
            assert (u->n == m);
        else
            assert (u->n == mn);
        
        if (mn == 0)
            _bcv_matrix_set_identity (u);
            
        u_data = u->data;
        ldu    = u->lda;
    }

    if (jobvt == BCV_MATRIX_SVDJOB_ALL || jobvt == BCV_MATRIX_SVDJOB_SOME) {
        _bcv_assert_valid_matrix (vt);
     
        assert (vt->n == n);
        
        if (jobu == BCV_MATRIX_SVDJOB_ALL)
            assert (vt->m == n);
        else
            assert (vt->m == mn);
        
        if (mn == 0)
            _bcv_matrix_set_identity (vt);
            
        vt_data = vt->data;
        ldvt    = vt->lda;
    }

    if (mn > 0)
    {
        assert (s);
        assert (work);
    
        _bcv_dgesvd_f77 (_BCV_F77_SVDJOB (jobu), _BCV_F77_SVDJOB (jobvt),
                         &m, &n, a->data, &lda, s, u_data, &ldu, vt_data, 
                         &ldvt, work, &lwork, &info);

        assert (info >= 0);
    }
    
    return info;
}


bcv_index_t
_bcv_lapack_dgesvd_work_len (bcv_matrix_svdjob_t jobu, 
                             bcv_matrix_svdjob_t jobvt, bcv_index_t m,
                             bcv_index_t n)
{
    bcv_index_t result = 0;
    bcv_index_t lda  = BCV_MAX (1, m);
    bcv_index_t ldu  = BCV_MAX (1, m);
    bcv_index_t ldvt = BCV_MAX (1, n);
    double work;
    bcv_index_t lwork = -1;
    bcv_error_t info;

    assert (m >= 0);
    assert (n >= 0);
    
    if (m > 0 && n > 0)
    {
        _bcv_dgesvd_f77 (_BCV_F77_SVDJOB (jobu), _BCV_F77_SVDJOB (jobvt),
                         &m, &n, NULL, &lda, NULL, NULL, &ldu, NULL, &ldvt,
                         &work, &lwork, &info);
        assert (info == 0);
    
        if (work <= (double) BCV_MAX_INDEX)
        {
            result = (bcv_index_t) work;
        }
    } 
    else 
    {
        result = 1;
    }
    
    return result;
}


bcv_error_t
_bcv_lapack_dgesdd (bcv_matrix_svdjob_t jobz,
                    bcv_matrix_t *a, double *s, bcv_matrix_t *u,
                    bcv_matrix_t *vt, double *work, bcv_index_t lwork,
                    bcv_index_t *iwork)
{
    bcv_index_t m, n, mn, lda, ldu = 1, ldvt = 1;
    bcv_error_t info = 0;
    double *u_data = NULL, *vt_data = NULL;
    
    _bcv_assert_valid_matrix (a);
    
    m   = a->m;
    n   = a->n;
    mn  = BCV_MIN (m, n);
    lda = a->lda;
    
    if (jobz == BCV_MATRIX_SVDJOB_ALL || jobz == BCV_MATRIX_SVDJOB_SOME)
    {
        _bcv_assert_valid_matrix (u);
        _bcv_assert_valid_matrix (vt);
        
        assert (u->m  == m);
        assert (vt->n == n);
        
        if (jobz == BCV_MATRIX_SVDJOB_ALL)
        {
            assert (u->n  == m);
            assert (vt->m == n);
        }
        else
        {
            assert (u->n  == mn);
            assert (vt->m == mn);
        }
        
        if (mn == 0)
        {
            _bcv_matrix_set_identity (u);
            _bcv_matrix_set_identity (vt);
        }
        else
        {
            u_data  = u->data;
            ldu     = u->lda;
            vt_data = vt->data;
            ldvt    = vt->lda;
        }
    }
    else if (jobz == BCV_MATRIX_SVDJOB_OVERWRITE)
    {
        if (m >= n)
        {
            _bcv_assert_valid_matrix (vt);
            assert (vt->m == n);
            assert (vt->n == n);
            
            vt_data = vt->data;
            ldvt    = vt->lda;
        }
        else
        {
            _bcv_assert_valid_matrix (u);
            assert (u->m == m);
            assert (u->n == m);
            
            u_data = u->data;
            ldu    = u->lda;
        }
    }

    if (mn > 0)
    {
        assert (s);
        assert (work);
    
        _bcv_dgesdd_f77 (_BCV_F77_SVDJOB (jobz),
                         &m, &n, a->data, &lda, s, u_data, &ldu, vt_data, 
                         &ldvt, work, &lwork, iwork, &info);

        assert (info >= 0);
    }
    
    return info;
}


bcv_index_t
_bcv_lapack_dgesdd_work_len (bcv_matrix_svdjob_t jobz, 
                             bcv_index_t m, bcv_index_t n)
{
    bcv_index_t result = 0;
    bcv_index_t lda  = BCV_MAX (1,m);
    bcv_index_t ldu  = BCV_MAX (1,m);
    bcv_index_t ldvt = BCV_MAX (1,n);
    double work;
    bcv_index_t lwork = -1;
    bcv_error_t info;
    
    assert (m >= 0);
    assert (n >= 0);
    
    if (m > 0 && n > 0)
    {
        _bcv_dgesdd_f77 (_BCV_F77_SVDJOB (jobz),
                         &m, &n, NULL, &lda, NULL, NULL, &ldu, NULL, &ldvt,
                         &work, &lwork, NULL, &info);
        assert (info == 0);
    
        if (work <= (double) BCV_MAX_INDEX)
        {
            result = (bcv_index_t) work;
        }
    } 
    else 
    {
        result = 1;
    }
    
    return result;
}

bcv_index_t
_bcv_lapack_dgesdd_iwork_len (bcv_index_t m, bcv_index_t n)
{
    bcv_index_t result = 0;
    
    assert (m >= 0);
    assert (n >= 0);
    
    result = BCV_MAX (1, 8 * BCV_MIN (m, n));
    
    return result;
}


void
_bcv_lapack_dgebrd (bcv_matrix_t *a, double *d, double *e, double *tauq,
                    double *taup, double *work, bcv_index_t lwork)
{
    bcv_error_t info  = 0;
    
    _bcv_assert_valid_matrix (a);
    
    if (a->m > 0 && a->n > 0)
    {
        assert (d);
        assert (e);
        assert (tauq);
        assert (taup);
        assert (work);
        
        _bcv_dgebrd_f77 (&(a->m), &(a->n), a->data, &(a->lda), 
                         d, e, tauq, taup, work, &lwork, &info);
        assert (info == 0);
    }
}


bcv_index_t
_bcv_lapack_dgebrd_work_len (bcv_index_t m, bcv_index_t n)
{
    bcv_error_t info  = 0;
    bcv_index_t lwork = -1;
    double size;
    bcv_index_t result = 0;
    
    assert (m >= 0);
    assert (n >= 0);    
    
    if (m > 0 && n > 0)
    {
        _bcv_dgebrd_f77 (&m, &n, NULL, &m,
                         NULL, NULL, NULL, NULL, &size, &lwork, &info);
        assert (info == 0);
        
        if (size <= (double) BCV_MAX_INDEX)
        {
            result = (bcv_index_t) size;
        }
    }
    else
    {
        result = 1;
    }

    return result;
}


void
_bcv_lapack_dormbr (bcv_matrix_vect_t vect, bcv_matrix_side_t side,
                    bcv_matrix_transpose_t trans, 
                    const bcv_matrix_t *a, double *tau,
                    bcv_matrix_t *c, double *work, bcv_index_t lwork)
{
    bcv_index_t k;
    bcv_error_t info = 0;
    
    _bcv_assert_valid_matrix (a);
    _bcv_assert_valid_matrix (c);
    
    assert (((vect == BCV_MATRIX_VECT_Q) ? a->m : a->n)
            == ((side == BCV_MATRIX_LEFT) ? c->m : c->n));

    k = (vect == BCV_MATRIX_VECT_Q) ? a->n : a->m;

    if (k > 0 && c->m > 0 && c->n > 0)
    {
        _bcv_dormbr_f77 (_BCV_F77_VECT (vect), _BCV_F77_SIDE (side),
                         _BCV_F77_TRANS (trans), &(c->m), &(c->n), &k,
                         a->data, &(a->lda), tau, c->data, &(c->lda),
                         work, &lwork, &info);
        assert (info == 0);
    }
}


bcv_index_t
_bcv_lapack_dormbr_work_len (bcv_matrix_vect_t vect, bcv_matrix_side_t side,
                             bcv_index_t ma, bcv_index_t na,
                             bcv_index_t mc, bcv_index_t nc)
{
    bcv_index_t result = 0;
    bcv_matrix_transpose_t trans = BCV_MATRIX_NOTRANS;
    bcv_index_t r   = (side == BCV_MATRIX_LEFT) ? mc : nc;
    bcv_index_t k   = (vect == BCV_MATRIX_VECT_Q) ? na : ma;
    bcv_index_t lda = (vect == BCV_MATRIX_VECT_Q) ? BCV_MAX (1,r) 
                                                  : BCV_MAX (1, BCV_MIN (r,k));
    bcv_index_t ldc = BCV_MAX (1,mc);
    double work;
    bcv_index_t lwork = -1;
    bcv_error_t info;

    assert (ma >= 0);
    assert (na >= 0);
    assert (mc >= 0);
    assert (nc >= 0);
    
    if (mc > 0 && nc > 0 && k > 0 && r > 0)
    {
        _bcv_dormbr_f77 (_BCV_F77_VECT (vect), _BCV_F77_SIDE (side),
                         _BCV_F77_TRANS (trans), &mc, &nc, &k,
                         NULL, &lda, NULL, NULL, &ldc,
                         &work, &lwork, &info);
        assert (info == 0);
    
        if (work <= (double) BCV_MAX_INDEX)
        {
            result = (bcv_index_t) work;
        }
    }
    else
    {
        result = 1;
    }
    
    return result;
}


bcv_error_t
_bcv_lapack_dbdsqr (bcv_matrix_uplo_t uplo, bcv_index_t n,
                    double *d, double *e, bcv_matrix_t *vt, bcv_matrix_t *u,
                    bcv_matrix_t *c, double *work)
{
    double *vt_data = NULL, *u_data = NULL, *c_data = NULL;
    bcv_index_t ncvt = 0, nru = 0, ncc = 0, ldvt = 1, ldu = 1, ldc = 1;
    bcv_error_t info = 0;
    
    if (n > 0)
    {
        assert (d);
        assert (e);
        
        if (vt)
        {
            _bcv_assert_valid_matrix (vt);
            assert (n == vt->m);
            
            ncvt    = vt->n;
            vt_data = vt->data;
            ldvt    = vt->lda;
        }
        if (u)
        {
            _bcv_assert_valid_matrix (u);
            assert (n == u->n);
            
            nru    = u->m;
            u_data = u->data;
            ldu    = u->lda;
        }
        if (c)
        {
            _bcv_assert_valid_matrix (c);
            assert (n == c->m);
            
            ncc    = c->n;
            c_data = c->data;
            ldc    = c->lda;
        }
        
        _bcv_dbdsqr_f77 (_BCV_F77_UPLO (uplo), &n, &ncvt, &nru, &ncc,
                         d, e, vt_data, &ldvt, u_data, &ldu, c_data, &ldc,
                         work, &info);
                         
        assert (info >= 0);
    }
    
    return info;
}

bcv_index_t
_bcv_lapack_dbdsqr_work_len (bcv_index_t n, bcv_bool_t only_values)
{
    bcv_index_t result = 0;
    
    assert (n >= 0);
    
    if (only_values) 
    {
        if (n <= BCV_MAX_INDEX / 2)
        {
            result = BCV_MAX (2 * n, 1);
        }
    }
    else
    {
        if (n <= BCV_MAX_INDEX / 4) 
        {
            result = BCV_MAX (4 * n, 1);
        }
    }
    
    return result;
}
