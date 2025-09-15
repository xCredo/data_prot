#include <stdio.h>
#include <gmp.h>

// Быстрое возведение в степень по модулю
void fast_pow(mpz_t res, const mpz_t a, const mpz_t exp, const mpz_t mod) {
    mpz_t base, exponent;
    mpz_init_set(base, a);
    mpz_init_set(exponent, exp);
    mpz_set_ui(res, 1);
    
    mpz_mod(base, base, mod);
    
    while (mpz_cmp_ui(exponent, 0) > 0) {
        if (mpz_odd_p(exponent)) {
            mpz_mul(res, res, base);
            mpz_mod(res, res, mod);
        }
        mpz_mul(base, base, base);
        mpz_mod(base, base, mod);
        mpz_fdiv_q_ui(exponent, exponent, 2);
    }
    
    mpz_clear(base);
    mpz_clear(exponent);
}

// Тест Ферма на простоту
int fermat_test(const mpz_t n, int k) {
    if (mpz_cmp_ui(n, 2) == 0) return 1;
    if (mpz_cmp_ui(n, 1) <= 0 || mpz_even_p(n)) return 0;
    
    mpz_t a, n_minus_1, res;
    mpz_init(a);
    mpz_init(n_minus_1);
    mpz_init(res);
    mpz_sub_ui(n_minus_1, n, 1);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    
    for (int i = 0; i < k; i++) {
        mpz_urandomm(a, state, n_minus_1);
        if (mpz_cmp_ui(a, 2) < 0) mpz_set_ui(a, 2);
        
        fast_pow(res, a, n_minus_1, n);
        if (mpz_cmp_ui(res, 1) != 0) {
            mpz_clear(a);
            mpz_clear(n_minus_1);
            mpz_clear(res);
            gmp_randclear(state);
            return 0;
        }
    }
    
    mpz_clear(a);
    mpz_clear(n_minus_1);
    mpz_clear(res);
    gmp_randclear(state);
    return 1;
}

// Расширенный алгоритм Евклида
void extended_gcd(mpz_t gcd, mpz_t x, mpz_t y, const mpz_t a, const mpz_t b) {
    if (mpz_cmp_ui(b, 0) == 0) {
        mpz_set(gcd, a);
        mpz_set_ui(x, 1);
        mpz_set_ui(y, 0);
        return;
    }
    
    mpz_t x1, y1;
    mpz_init(x1);
    mpz_init(y1);
    mpz_t mod;
    mpz_init(mod);
    
    mpz_mod(mod, a, b);
    extended_gcd(gcd, x1, y1, b, mod);
    
    mpz_set(x, y1);
    mpz_t temp;
    mpz_init(temp);
    mpz_fdiv_q(temp, a, b);
    mpz_mul(temp, temp, y1);
    mpz_sub(y, x1, temp);
    
    mpz_clear(x1);
    mpz_clear(y1);
    mpz_clear(mod);
    mpz_clear(temp);
}

int main() {
    mpz_t a, b, res, x, y;
    mpz_init(a);
    mpz_init(b);
    mpz_init(res);
    mpz_init(x);
    mpz_init(y);
    
    int choice;
    printf("1. Fast exponentiation\n");
    printf("2. Fermat test\n");
    printf("3. Extended GCD\n");
    printf("Choose: ");
    scanf("%d", &choice);
    
    if (choice == 1) {
        gmp_scanf("%Zd %Zd %Zd", a, b, res);
        fast_pow(x, a, b, res);
        gmp_printf("Result: %Zd\n", x);
    }
    else if (choice == 2) {
        gmp_scanf("%Zd", a);
        int is_prime = fermat_test(a, 5);
        printf("Is prime: %d\n", is_prime);
    }
    else if (choice == 3) {
        gmp_scanf("%Zd %Zd", a, b);
        extended_gcd(res, x, y, a, b);
        gmp_printf("GCD: %Zd\n", res);
        gmp_printf("x: %Zd, y: %Zd\n", x, y);
    }
    
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(res);
    mpz_clear(x);
    mpz_clear(y);
    
    return 0;
}
