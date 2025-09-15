#include <stdio.h>
#include <gmp.h>
#include <math.h>

// Быстрое возведения в степень
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

// Алгоритм Шаг младенца - шаг великана
void baby_step_giant_step(mpz_t result, const mpz_t a, const mpz_t y, const mpz_t p) {
    mpz_t m, j, current, a_inv, a_inv_m, temp, i_val;
    mpz_init(m);
    mpz_init(j);
    mpz_init(current);
    mpz_init(a_inv);
    mpz_init(a_inv_m);
    mpz_init(temp);
    mpz_init(i_val);
    
    // Вычисляем m = ceil(sqrt(p-1))
    mpz_sub_ui(temp, p, 1);
    double sqrt_val = sqrt(mpz_get_d(temp));
    mpz_set_ui(m, (long)ceil(sqrt_val));
    
    // Baby steps - сохраняем a^j mod p для j от 0 до m-1
    mpz_set_ui(current, 1);
    mpz_t *table = malloc(mpz_get_ui(m) * sizeof(mpz_t));
    
    for (mpz_set_ui(j, 0); mpz_cmp(j, m) < 0; mpz_add_ui(j, j, 1)) {
        mpz_init(table[mpz_get_ui(j)]);
        mpz_set(table[mpz_get_ui(j)], current);
        mpz_mul(current, current, a);
        mpz_mod(current, current, p);
    }
    
    // Находим a^(-m) mod p
    mpz_t x, y_temp;
    mpz_init(x);
    mpz_init(y_temp);
    extended_gcd(temp, x, y_temp, a, p);
    fast_pow(a_inv_m, x, m, p);
    
    // Giant steps
    mpz_set(current, y);
    for (mpz_set_ui(i_val, 0); mpz_cmp(i_val, m) < 0; mpz_add_ui(i_val, i_val, 1)) {
        // Ищем current в таблице
        for (mpz_set_ui(j, 0); mpz_cmp(j, m) < 0; mpz_add_ui(j, j, 1)) {
            if (mpz_cmp(current, table[mpz_get_ui(j)]) == 0) {
                // Нашли решение: x = i*m + j
                mpz_mul(temp, i_val, m);
                mpz_add(result, temp, j);
                
                // Очистка
                for (unsigned long k = 0; k < mpz_get_ui(m); k++) {
                    mpz_clear(table[k]);
                }
                free(table);
                mpz_clear(m);
                mpz_clear(j);
                mpz_clear(current);
                mpz_clear(a_inv);
                mpz_clear(a_inv_m);
                mpz_clear(temp);
                mpz_clear(i_val);
                mpz_clear(x);
                mpz_clear(y_temp);
                return;
            }
        }
        mpz_mul(current, current, a_inv_m);
        mpz_mod(current, current, p);
    }
    
    // Решение не найдено
    mpz_set_si(result, -1);
    
    // Очистка
    for (unsigned long k = 0; k < mpz_get_ui(m); k++) {
        mpz_clear(table[k]);
    }
    free(table);
    mpz_clear(m);
    mpz_clear(j);
    mpz_clear(current);
    mpz_clear(a_inv);
    mpz_clear(a_inv_m);
    mpz_clear(temp);
    mpz_clear(i_val);
    mpz_clear(x);
    mpz_clear(y_temp);
}

// Расширенный алгоритм Евклида (нужен для поиска обратного элемента)
void extended_gcd(mpz_t gcd, mpz_t x, mpz_t y, const mpz_t a, const mpz_t b) {
    if (mpz_cmp_ui(b, 0) == 0) {
        mpz_set(gcd, a);
        mpz_set_ui(x, 1);
        mpz_set_ui(y, 0);
        return;
    }
    
    mpz_t x1, y1, mod;
    mpz_init(x1);
    mpz_init(y1);
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
    mpz_t a, y, p, result;
    mpz_init(a);
    mpz_init(y);
    mpz_init(p);
    mpz_init(result);
    
    printf("Baby-Step Giant-Step Algorithm\n");
    printf("Solve a^x ≡ y mod p\n");
    
    gmp_printf("Enter a: ");
    gmp_scanf("%Zd", a);
    gmp_printf("Enter y: ");
    gmp_scanf("%Zd", y);
    gmp_printf("Enter p: ");
    gmp_scanf("%Zd", p);
    
    baby_step_giant_step(result, a, y, p);
    
    if (mpz_cmp_si(result, -1) == 0) {
        printf("No solution found\n");
    } else {
        gmp_printf("Solution found: x = %Zd\n", result);
        
        // Проверка
        mpz_t check;
        mpz_init(check);
        fast_pow(check, a, result, p);
        gmp_printf("Check: a^x mod p = %Zd\n", check);
        gmp_printf("Expected: y = %Zd\n", y);
        mpz_clear(check);
    }
    
    mpz_clear(a);
    mpz_clear(y);
    mpz_clear(p);
    mpz_clear(result);
    
    return 0;
}
