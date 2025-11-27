/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
    /**
     * Epsilon value for numerical computation, it serves as a "close enough" threshold.
     */
    public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
    /**
     * The zero polynomial function is represented as an array with a single (0) entry.
     */
    public static final double[] ZERO = {0};

    /**
     * Computes the f(x) value of the polynomial function at x.
     *
     * @param poly - polynomial function
     * @param x
     * @return f(x) - the polynomial function value at x.
     */
    public static double f(double[] poly, double x) {
        double ans = 0;
        for (int i = 0; i < poly.length; i++) {
            double c = Math.pow(x, i);
            ans += c * poly[i];
        }
        return ans;
    }

    /**
     * Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
     * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps,
     * assuming p(x1)*p(x2) <= 0.
     * This function should be implemented recursively.
     *
     * @param p   - the polynomial function
     * @param x1  - minimal value of the range
     * @param x2  - maximal value of the range
     * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
     * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
     */
    public static double root_rec(double[] p, double x1, double x2, double eps) {
        double f1 = f(p, x1);
        double x12 = (x1 + x2) / 2;
        double f12 = f(p, x12);
        if (Math.abs(f12) < eps) {
            return x12;
        }
        if (f12 * f1 <= 0) {
            return root_rec(p, x1, x12, eps);
        } else {
            return root_rec(p, x12, x2, eps);
        }
    }

    /**
     * This function computes a polynomial representation from a set of 2D points on the polynom.
     * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
     * Note: this function only works for a set of points containing up to 3 points, else returns null.
     *
     * @param xx
     * @param yy
     * @return an array of doubles representing the coefficients of the polynom.
     */
    public static double[] PolynomFromPoints(double[] xx, double[] yy) {
        double[] ans = null;
        int lx = xx.length;
        int ly = yy.length;
        double x1, x2, x3, y1, y2, y3, denom, a, b, c, m;
        if (xx != null && yy != null && lx == ly && lx > 1 && lx < 4) {

            if (lx == 3) {
                x1 = xx[0];
                x2 = xx[1];
                y1 = yy[0];
                y2 = yy[1];
                y3 = yy[2];
                x3 = xx[2];
                denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
                a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
                b = (Math.pow(x3, 2) * (y1 - y2) + Math.pow(x2, 2) * (y3 - y1) + Math.pow(x1, 2) * (y2 - y3)) / denom;
                c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;
                ans = new double[]{c, b, a};

            }
            if (lx == 2) {
                x1 = xx[0];
                x2 = xx[1];
                y1 = yy[0];
                y2 = yy[1];
                m = (y1 - y2) / (x1 - x2);
                a = m;
                b = y1 - (x1 * m);
                ans = new double[]{b, a};
            }
        }
        return ans;
    }

    /**
     * Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
     * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
     *
     * @param p1 first polynomial function
     * @param p2 second polynomial function
     * @return true iff p1 represents the same polynomial function as p2.
     */
    public static boolean equals(double[] p1, double[] p2) {
        boolean ans = true;
        double [] fix1= fix0(p1);
        double [] fix2= fix0(p2);
        int powP1 = fix1.length, powP2 = fix2.length;
        if (powP1 != powP2)
            return false;
        for (int i = 0; i < powP1 + 1; i++) {
            if (f(fix1, i) != f(fix2, i)) {
                return false;
            }
        }
        return ans;
    }

    /**
     * Computes a String representing the polynomial function.
     * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
     *
     * @param poly the polynomial function represented as an array of doubles
     * @return String representing the polynomial function:
     */
    public static String poly(double[] poly) {
        String ans = "";
        if (poly.length == 0) {
            ans = "0";
        } else {
            ans += poly[poly.length - 1] + "x" + "^" + (poly.length - 1);
            for (int i = 1; i < poly.length - 1; i++) {
                if (poly.length - 2 > i) {
                    if (poly[poly.length - i - 1] < 0) {
                        ans += " " + poly[i] + "x" + "^" + (poly.length - i - 1);
                    } else {
                        ans += " +" + poly[poly.length - i - 1] + "x" + "^" + (poly.length - i - 1);
                    }
                } else {
                    if (poly[poly.length - i - 1] < 0) {
                        ans += " " + poly[1] + "x";
                    } else {
                        ans += " +" + poly[1] + "x";
                    }
                    if (poly[poly.length - i - 2] < 0) {
                        ans += " " + poly[0];
                    } else {
                        ans += " +" + poly[0];
                    }

                }
            }

        }
        return ans;
    }

    /**
     * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
     * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
     *
     * @param p1  - first polynomial function
     * @param p2  - second polynomial function
     * @param x1  - minimal value of the range
     * @param x2  - maximal value of the range
     * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
     * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
     */
    public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        double ans = x1, chek;

        for (double i = x1; i <= x2; i += EPS) {
            chek = Math.abs(f(p1, i) - f(p2, i));
            if (Math.abs(f(p1, i) - f(p2, i)) <= 0.01) {
                return i;
            }
        }
        return ans;


    }

    /**
     * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
     * This function computes an approximation of the length of the function between f(x1) and f(x2)
     * using n inner sample points and computing the segment-path between them.
     * assuming x1 < x2.
     * This function should be implemented iteratively (none recursive).
     *
     * @param p                - the polynomial function
     * @param x1               - minimal value of the range
     * @param x2               - maximal value of the range
     * @param numberOfSegments - (A positive integer value (1,2,...).
     * @return the length approximation of the function between f(x1) and f(x2).
     */
    public static double length(double[] p, double x1, double x2, int numberOfSegments) {
        double ans = 0,l;
        l=Math.abs(x1 - x2)/numberOfSegments;
        for(double i=x1;i<x2;i+=l)
        {
            ans+=Math.sqrt(Math.pow((l),2)+Math.pow((f(p,i)-f(p,i+l)),2));
        }
        return ans;
    }

    /**
     * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
     * This function computes an approximation of the area between the polynomial functions within the x-range.
     * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
     *
     * @param p1                - first polynomial function
     * @param p2                - second polynomial function
     * @param x1                - minimal value of the range
     * @param x2                - maximal value of the range
     * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
     * @return the approximated area between the two polynomial functions within the [x1,x2] range.
     */
    public static double area(double[] p1, double[] p2, double x1, double x2, int numberOfTrapezoid) {
        double ans = 0, hT = (Math.abs(x1) + Math.abs(x2)) / numberOfTrapezoid;
        double d1, d2, sameX;

        for (double i = x1; i < x2; i += hT) {
            d1 = Math.round(Math.abs(f(p1, i) - f(p2, i)) * 1000) * EPS;
            d2 = Math.round(Math.abs(f(p1, i + hT) - f(p2, i + hT)) * 1000) * EPS;
            if (sameValue(p1, p2, i, (i + hT), EPS) != i) {
                sameX = sameValue(p1, p2, i, (i + hT), EPS);
                ans += (d1 * (Math.abs(sameX - i))) / 2;
                ans += (d2 * (Math.abs(sameX - (i + hT)))) / 2;
            } else {
                ans += ((d1 + d2) * hT) / 2;
            }
        }
        return ans;
    }

    /**
     * This function computes the array representation of a polynomial function from a String
     * representation. Note:given a polynomial function represented as a double array,
     * getPolynomFromString(poly(p)) should return an array equals to p.
     *
     * @param p - a String representing polynomial function.
     * @return
     */
    public static double[] getPolynomFromString(String p) {
        double[] ans = ZERO;//  -1.0x^2 +3.0x +2.0
        int count = 0, deg = 0;
        char[] poly = p.toCharArray();
        int d = 0;
        boolean negative = false;
        for (int j = 0; j < poly.length; j++) {
            if (poly[j] == 'x' && poly[j + 1] == '^') {
                deg = Character.getNumericValue(poly[j + 2]);
                break;
            }
            if (poly[j] == 'x' && poly[j + 1] != '^') {
                deg = 1;
                break;
            }
        }
        ans = new double[deg + 1];
        if (deg >= 2) {
            for (int k = deg; k > 1; k--) {
                for (int i = count; i < poly.length; i++) {

                    if (poly[i] == '-') {
                        negative = true;
                        i++;
                    }
                    if (poly[i] == '.') {
                        i++;
                    }
                    if (poly[i] == 'x') {
                        count += i;
                        if (poly[i + 1] == '^') {
                            count = count + 4;
                        } else {
                            count = count + 2;
                        }
                        break;
                    }
                    ans[k] += (Character.getNumericValue(poly[i])) / (Math.pow(10, d));
                    d++;
                }
                if (negative) {
                    ans[k] *= -1;
                }
                negative = false;
                d = 0;
            }
        }
        if (deg >= 1) {
            for (int i = count; i < poly.length; i++) {
                if (poly[i] == '-') {
                    negative = true;
                    i++;
                }
                if (poly[i] == '+') {
                    i++;
                }
                if (poly[i] == '.') {
                    i++;
                }
                if (poly[i] == 'x') {
                    count = i;
                    if (poly[i + 1] == '^') {
                        count = count + 4;
                    } else {
                        count = count + 2;
                    }
                    break;
                }
                ans[1] += (Character.getNumericValue(poly[i]) / Math.pow(10, d));
                d++;
            }
            if (negative) {
                ans[1] *= -1;
            }
            negative = false;
            d = 0;
        }
        for (int i = count; i < poly.length; i++) {
            if (poly[i] == '-') {
                negative = true;
                i++;
            }
            if (poly[i] == '+') {
                i++;
            }
            if (poly[i] == '.') {
                i++;
            }
            ans[0] += (Character.getNumericValue(poly[i]) / Math.pow(10, d));
            d++;
        }
        if (negative) {
            ans[0] *= -1;
        }
        return ans;
    }

    /**
     * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
     *
     * @param p1
     * @param p2
     * @return
     */
    public static double[] add(double[] p1, double[] p2) {
        double[] ans = ZERO;//
        int l1 = p1.length, l2 = p2.length, maxl;
        maxl = Math.max(l1, l2);
        ans = new double[maxl];
        for (int i = 0; i < maxl; i++) {
            if (i >= l1) {
                ans[i] = p2[i];
            }
            if (i >= l2) {
                ans[i] = p1[i];
            }
            if (i < l1 && i < l2) {
                ans[i] = p1[i] + p2[i];
            }
        }
        return fix0(ans);
    }

    /**
     * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
     *
     * @param p1
     * @param p2
     * @return
     */
    public static double[] mul(double[] p1, double[] p2) {
        double[] ans = new double[(p1.length - 1) + (p2.length - 1) + 1];
        if (!(p1 != ZERO && p2 != ZERO)) {
            return ZERO;
        }
        int c = 0;
        for (int i = 0; i < p1.length; i++) {
            for (int j = 0; j < p2.length; j++) {
                ans[i + j] += p1[i] * p2[j];
                c++;
            }
        }
        return fix0(ans);
    }

    /**
     * This function computes the derivative of the p0 polynomial function.
     *
     * @param po
     * @return
     */
    public static double[] derivative(double[] po) {
        double[] ans = ZERO;//
        if (po.length == 1) {
            return ZERO;
        }
        double[] ptag = new double[po.length - 1];
        for (int i = 0; i < ptag.length; i++) {
            ptag[i] = po[i + 1] * (i + 1);
        }

        return ptag;
    }

    public static double[] fix0(double[] p) {
        int i = 0, waste = 0;
        for (int i1 = 0; i1 < p.length; i1++) {
            if (p[i1] == 0 ) {
                waste++;
            } else {
                waste = 0;
            }
        }
        if (waste == p.length) {
            return ZERO;
        }
        double[] fix = new double[p.length - waste];
        for (int j = 0; j < fix.length; j++) {
            fix[j] = fixr(p[j]);
        }
        return fix;
    }
    public static double fixr(double p)
    {
        double ans = p;
        ans -=(int)p;
        if (Math.abs(ans) < EPS) {return (int)p;}
        else {return p;}
    }



    }


