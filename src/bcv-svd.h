/* bcv-svd.h
 * 
 */
 
#ifndef _BCV_SVD_H
#define _BCV_SVD_H


typedef struct _bcv_svd bcv_svd_t;


bcv_svd_t *bcv_svd_alloc (int M_max, int N_max);
void bcv_svd_init (bcv_svd_t *bcv, int M, int N, int m, int n, double *x, int ldx);
void bcv_svd_initp (bcv_svd_t *bcv, int M, int N, int m, int n, double *x, int ldx, int *p, int *q);
void bcv_svd_free (bcv_svd_t *bcv);

int bcv_svd_decompose (bcv_svd_t *bcv);

void bcv_svd_get_resid (const bcv_svd_t *bcv, int *m2, int *n2, double **resid, int *ldr);
int bcv_svd_get_max_rank (bcv_svd_t *bcv);

void bcv_svd_update_resid (bcv_svd_t *bcv, double scale, int k);
double bcv_svd_get_resid_mse (const bcv_svd_t *bcv);

#endif /* _BCV_SVD_H */
