package math;

public class Complex {
    public static final Complex EIN = new Complex(0, 1);
    public static final Complex ZERO = new Complex(0, 0);

    private double real;
    private double image;

    public Complex() {
        this.real = 0;
        this.image = 0;
    }

    public Complex(double real, double image) {
        this.real = real;
        this.image = image;
    }

    public double getReal() {
        return real;
    }

    public void setReal(long real) {
        this.real = real;
    }

    public double getImage() {
        return image;
    }

    public void setImage(long image) {
        this.image = image;
    }

    public double getModulus() {
        return Math.sqrt(getSquaredModulus());
    }

    public double getSquaredModulus() {
        return Math.pow(real, 2) + Math.pow(image, 2);
    }

    public static Complex add(Complex c1, Complex c2) {
        return new Complex(c1.getReal() + c2.getReal(), c1.getImage() + c2.getImage());
    }

    public static Complex add(Complex... c) {
        Complex sum = c[0];
        for (int i = 1; i < c.length; i++) {
            sum = add(sum, c[i]);
        }
        return sum;
    }

    public static Complex add(double d, Complex c) {
        return new Complex(d + c.getReal(), c.getImage());
    }

    public Complex conjugate() {
        return new Complex(real, -image);
    }

    public static Complex conjugate(Complex c) {
        return new Complex(c.real, -c.image);
    }

    public static Complex subtract(Complex c1, Complex c2) {
        return add(c1, negative(c2));
    }

    public static Complex subtract(double d, Complex c) {
        return subtract(new Complex(d, 0), c);
    }

    public static Complex negative(Complex c) {
        return new Complex(-c.getReal(), -c.getImage());
    }

    public static Complex multiply(Complex c1, Complex c2) {
        return new Complex(c1.getReal() * c2.getReal() - c1.getImage() * c2.getImage(),
                c1.getReal() * c2.getImage() + c1.getImage() * c2.getReal());
    }

    public static Complex multiply(Complex... c) {
        Complex product = c[0];
        for (int i = 1; i < c.length; i++) {
            product = multiply(product, c[i]);
        }
        return product;
    }

    public static Complex multiply(double d, Complex c) {
        return new Complex(d * c.getReal(), d * c.getImage());
    }

    public static Complex sqr(Complex c) {
        return multiply(c, c);
    }

    public static Complex divide(Complex c1, Complex c2) {
        return divide(multiply(c1, new Complex(c2.getReal(), -c2.getImage())), c2.getSquaredModulus());
    }

    public static Complex divide(Complex c, double d) {
        return new Complex(c.getReal() / d, c.getImage() / d);
    }

    public static Complex exp(Complex c) {
        return new Complex(Math.exp(c.getReal()) * Math.cos(c.getImage()),
                Math.exp(c.getReal()) * Math.sin(c.getImage()));
    }

    @Override
    public String toString() {
        return "Real: " + real + ", Image: " + image;
    }

    @Override
    public boolean equals(Object obj) {
        Complex c = (Complex) obj;
        return (real == c.getReal() && image == c.getImage()) ? true : false;
    }
}