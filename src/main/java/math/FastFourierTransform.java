package math;

import static math.Complex.*;
import static java.lang.Math.*;

public class FastFourierTransform {
    public static Complex[] fft(Complex[] f) {
        int n = f.length;
        if (n == 1) {
            return new Complex[]{f[0]};
        }

        if (n % 2 != 0) {
            throw new IllegalArgumentException("Count of elements is not a power of 2");
        }

        Complex[] even = new Complex[n / 2];
        for (int i = 0; i < n / 2; i++) {
            even[i] = f[2 * i];
        }
        Complex[] leftImage = fft(even);

        Complex[] odd = new Complex[n / 2];
        for (int i = 0; i < n / 2; i++) {
            odd[i] = f[2 * i + 1];
        }
        Complex[] rightImage = fft(odd);

        Complex[] image = new Complex[n];
        for (int i = 0; i < n / 2; i++) {
            double exponent = -2 * i * PI / n;
            Complex term = exp(new Complex(0, exponent));
            image[i] = add(leftImage[i], multiply(term, rightImage[i]));
            image[i + n / 2] = subtract(leftImage[i], multiply(term, rightImage[i]));
        }
        return image;
    }

    public static Complex[] ifft(Complex[] f) {
        int n = f.length;
        Complex[] original = new Complex[n];
        for (int i = 0; i < n; i++) {
            original[i] = f[i].conjugate();
        }
        original = fft(original);
        for (int i = 0; i < n; i++) {
            original[i] = original[i].conjugate();
        }
        for (int i = 0; i < n; i++) {
            original[i] = divide(original[i], n);
        }
        return original;
    }

    public static Complex[] fftShift(Complex[] f) {
        int n = f.length;
        Complex[] shift = new Complex[n];
        for (int i = 0; i < n / 2; i++) {
            shift[i] = f[n / 2 + i];
            shift[n / 2 + i] = f[i];
        }
        return shift;
    }
}
