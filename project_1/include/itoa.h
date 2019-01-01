#ifndef MY_ITOA
#define MY_ITOA

#include <stdint.h>

/* implementation of itoa to convert an Integer to C-String
 * taken from: http://www.strudel.org.uk/itoa/
 */
char *my_itoa(int64_t val, char *res, int base)
{
    /* check if base is valid */
    if (base < 2 || base > 36) {
        *res = '\0';
        return res;
    }

    char *ptr = res, *ptr_orig = res, tmp_char;
    int64_t tmp_val;

    /* write val to ptr in reverse order */
    do {
        tmp_val = val;
        val /= base;
        *ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz" [35 + (tmp_val - val * base)];
    } while (val);

    /* negative sign */
    if (tmp_val < 0)
        *ptr++ = '-';
    /* add NULL terminator */
    *ptr-- = '\0';

    /* revert string */
    while (ptr_orig < ptr) {
        tmp_char = *ptr;
        *ptr-- = *ptr_orig;
        *ptr_orig++ = tmp_char;
    }

    return res;
}

#endif
