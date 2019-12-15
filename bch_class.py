import numpy as np

from pyfinite import ffield, genericmatrix


def degree(p):
    poly = np.poly1d(np.flipud(p))
    return poly.order


class GF(ffield.FField):
    def PolyToString(self, poly):
        if type(poly) is int:
            poly = np.array([poly])
        elif type(poly) is str:
            return poly
        s = ''

        index_list = range(0, poly.size)

        for i in index_list:
            coeff_element = poly[i]
            if coeff_element != 0:
                if s != '':
                    s += ' + '
                if coeff_element > 1:
                    s += self.ElementToString(coeff_element)
                    if i > 0:
                        s += '*'
                if i > 0:
                    s += 'X'
                    if i > 1:
                        s += '^' + str(i)
                elif coeff_element == 1:
                    s += '1'

        if s == '':
            s = '0'

        return s

    def ElementToString(self, a):
        return '(' + self.ShowPolynomial(a) + ')'

    def SolveMatrix(self, A, b):
        return b

    def Pow(self, a, pow):
        res = 1
        for i in range(pow):
            res = self.Multiply(res, a)

        return res

    def GetPrimitive(self):
        for a in range(2, 2 ** self.n):
            is_primitive = True
            for power in range(2, 2 ** self.n - 1):  # no need for 2^n (a^2^n == 1)
                if self.Pow(a, power) == 0:
                    is_primitive = False
                    break

            if is_primitive:
                return a

    def FindRoots(self, poly):
        roots = []
        for a in range(2 ** self.n):
            result = self.SubstituteElementIntoPoly(poly, a)
            if result == 0:
                roots.append(a)
        return roots

    def SubstituteElementIntoPoly(self, poly, a):
        result = 0
        for i in range(len(poly)):
            tmp = self.Pow(a, i)
            tmp = self.Multiply(tmp, poly[i])
            result = self.Add(result, tmp)

        return result

    def ConjugateRoots(self, root):
        roots = []
        s = 1
        new_root = 1
        while new_root != root:
            new_root = self.Pow(root, 2**s)
            roots.append(new_root)
            s += 1

        return sorted(roots)

    def MinimalPoly(self, roots):
        result = 1
        for root in roots:
            root_poly = np.array([root, 1]) # (root + X)
            result = self.MultPoly(result, root_poly)
        return result.astype(int)

    def AddPoly(self, poly1, poly2):
        if type(poly1) is int:
            poly1 = np.array([poly1])
        if type(poly2) is int:
            poly2 = np.array([poly2])

        summ_degree = max(degree(poly1), degree(poly2))
        summ = np.zeros(summ_degree + 1)

        for i in range(0, summ.size):
            if i > poly1.size - 1:
                summ[i] = poly2[i]
            elif i > poly2.size - 1:
                summ[i] = poly1[i]
            else:
                summ[i] = self.Add(poly1[i], poly2[i])

        summ = summ[: degree(summ) + 1]  # remove dangling zeros

        return summ.astype(int)

    def MultPoly(self, poly1, poly2):
        if type(poly1) is int:
            poly1 = np.array([poly1])
        if type(poly2) is int:
            poly2 = np.array([poly1])

        product_max_degree = degree(poly1) + degree(poly2)
        product = np.zeros(product_max_degree + 1)

        for i in range(degree(poly1) + 1):
            for j in range(degree(poly2) + 1):
                product_tmp = self.Multiply(poly1[i], poly2[j])
                product[i + j] = self.Add(int(product[i + j]), product_tmp)

        product = product[: degree(product) + 1]  # remove dangling zeros
        return product.astype(int)

    def DivModPoly(self, poly1, poly2):
        if type(poly1) is int:
            poly1 = np.array([poly1])
        if type(poly2) is int:
            poly2 = np.array([poly2])

        assert np.count_nonzero(poly2) > 0, "Division by zero!"

        q = np.zeros(0)

        while degree(poly1) >= degree(poly2):
            exp_a = degree(poly1)
            exp_b = degree(poly2)
            coeff_element_a = poly1[exp_a]
            coeff_element_b = poly2[exp_b]
            multiplier = np.zeros(degree(poly1) - degree(poly2) + 1)
            multiplier[exp_a - exp_b] = self.Divide(coeff_element_a, coeff_element_b)
            subtrahend = self.MultPoly(poly2, multiplier.astype(int))
            poly1 = self.AddPoly(poly1, subtrahend.astype(int))  # addition = subtraction
            q = self.AddPoly(q.astype(int), multiplier.astype(int))

        remainder = poly1
        return [q, remainder]

    def DivPoly(self, a, b):
        return self.DivModPoly(a, b)[0]

    def ModPoly(self, a, b):
        return self.DivModPoly(a, b)[1]


GF2 = GF(1)


class BCH:
    def __init__(self, GF, t):
        self.t = t
        self.GF = GF
        self.n = 2 ** self.GF.n - 1

        self.primitive = self.GF.GetPrimitive()

        self.generator = self.GetGenerator()

        self.k = self.n - degree(self.generator)

    def GetGenerator(self):
        power_list = list(range(1, 2 * self.t + 1, 2))

        print()
        print('Generator polynomial g(X) of BCH code with t = ' + str(self.t) + ':')
        print()

        total_roots = []
        for power in power_list:
            current_cong_roots = self.GF.ConjugateRoots(self.GF.Pow(self.primitive, power))
            total_roots += current_cong_roots

        g = self.GF.MinimalPoly(set(total_roots))

        print('g(X) = ' + self.GF.PolyToString(g))
        print()

        return g

    def Encode(self, word):
        r = degree(self.generator)   # r = n - k

        # Xr = X^r
        Xr = np.zeros(r + 1)
        Xr[r] = 1
        Xr = Xr.astype(int)

        XrmX = GF2.MultPoly(Xr, word)  # X^(n-k) * word(X)
        p = GF2.ModPoly(XrmX, self.generator)  # p(X) = (X^(n-k) * word(X)) mod g(X)
        c = GF2.AddPoly(p, XrmX).tolist()  # c(X) = p(X) + (X^(n-k) * m(X))
        #
        # c = GF2.MultPoly(word, self.generator).tolist()

        c = np.array(c + [0] * (self.n - len(c)))

        return c.astype(int)

    def GetSyndrome(self, r):
        syndrome = np.zeros(2 * self.t)

        for power in range(1, 2 * self.t + 1):  # 1 <= i <= 2t
            a = self.GF.Pow(self.primitive, power)
            s = self.GF.SubstituteElementIntoPoly(r, a)
            syndrome[power - 1] = s

        return syndrome.astype(int)

    def Decode(self, r):
        # print('\tcode to decode: ', r)

        s = self.GetSyndrome(r)
        # print('\tsyndrome:', s)

        if s.tolist() == np.zeros(2 * self.t).tolist():
            return r[-self.k:]

        errors = None
        for n_errors in reversed(range(1, self.t + 1)):  # 0 already checked (should me minimum of 1 error)
            syndrome_matrix = genericmatrix.GenericMatrix(
                size=(n_errors, n_errors),
                zeroElement=0,
                identityElement=1,
                add=self.GF.Add,
                mul=self.GF.Multiply,
                sub=self.GF.Subtract,
                div=self.GF.Divide,
            )

            for row in range(n_errors):
                syndrome_matrix.SetRow(row, s[row: row + n_errors])

            if syndrome_matrix.Determinant() != 0:
                # print('\tn errors: ', n_errors)
                b = s[n_errors: 2 * n_errors]
                # print('\tmatrix: ', syndrome_matrix)
                # print('\tb: ', b)

                locators = np.array(syndrome_matrix.Solve(b)).astype(int).tolist()[:: -1]
                locators_poly = np.array([1] + locators).astype(int)

                # print('\tlocators: ', locators)
                # print('\tlocators poly: ', locators_poly)

                roots = self.GF.FindRoots(locators_poly)
                # print('\troots: ', roots)
                X_i = [self.GF.Inverse(root) for root in roots]
                # print('\tX_i: ', X_i)

                errors = np.zeros(self.n)

                for x in X_i:
                    for power in range(self.n):
                        if self.GF.Pow(self.primitive, power) == x:
                            errors[power] = 1
                            break

                # print('\terrors: ', errors)
                break

        if errors is None:
            return r[-self.k:]

        v = GF2.AddPoly(r, errors.astype(int)).tolist()
        v = np.array(v + [0] * (self.n - len(v)))

        return v[-self.k:]

