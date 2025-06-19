from functools import reduce

q = 3329
zetas = [1, 1729, 2580, 3289, 2642, 630, 1897, 848, 1062, 1919, 193, 797, 2786, 3260, 569, 1746, 296, 2447, 1339, 1476, 3046, 56, 2240, 1333, 1426, 2094, 535, 2882, 2393, 2879, 1974, 821, 289, 331, 3253, 1756, 1197, 2304, 2277, 2055, 650, 1977, 2513, 632, 2865, 33, 1320, 1915, 2319, 1435, 807, 452, 1438, 2868, 1534, 2402, 2647, 2617, 1481, 648, 2474, 3110, 1227, 910, 17, 2761, 583, 2649, 1637, 723, 2288, 1100, 1409, 2662, 3281, 233, 756, 2156, 3015, 3050, 1703, 1651, 2789, 1789, 1847, 952, 1461, 2687, 939, 2308, 2437, 2388, 733, 2337, 268, 641, 1584, 2298, 2037, 3220, 375, 2549, 2090, 1645, 1063, 319, 2773, 757, 2099, 561, 2466, 2594, 2804, 1092, 403, 1026, 1143, 2150, 2775, 886, 1722, 1212, 1874, 1029, 2110, 2935, 885, 2154]
zetas_2 = [17, 3312, 2761, 568, 583, 2746, 2649, 680, 1637, 1692, 723, 2606, 2288, 1041, 1100, 2229, 1409, 1920, 2662, 667, 3281, 48, 233, 3096, 756, 2573, 2156, 1173, 3015, 314, 3050, 279, 1703, 1626, 1651, 1678, 2789, 540, 1789, 1540, 1847, 1482, 952, 2377, 1461, 1868, 2687, 642, 939, 2390, 2308, 1021, 2437, 892, 2388, 941, 733, 2596, 2337, 992, 268, 3061, 641, 2688, 1584, 1745, 2298, 1031, 2037, 1292, 3220, 109, 375, 2954, 2549, 780, 2090, 1239, 1645, 1684, 1063, 2266, 319, 3010, 2773, 556, 757, 2572, 2099, 1230, 561, 2768, 2466, 863, 2594, 735, 2804, 525, 1092, 2237, 403, 2926, 1026, 2303, 1143, 2186, 2150, 1179, 2775, 554, 886, 2443, 1722, 1607, 1212, 2117, 1874, 1455, 1029, 2300, 2110, 1219, 2935, 394, 885, 2444, 2154, 1175]

def NTT(f):
    """
    Aplica la transformada NTT al polinomio f en R_q.

    Esta función implementa la NTT sobre un polinomio de 256 coeficientes,
    según el Algoritmo 9 del estándar FIPS 203. El resultado es el polinomio evaluado en T_q, usando
    constantes precomputadas zeta en orden BitRev7(i).

    Entrada:
    - f: lista de 256 enteros módulo q (coeficientes del polinomio en R_q).

    Salida:
    - Lista de 256 enteros módulo q (representación del polinomio en T_q).
    """
    assert(len(f) == 256)

    f_gorro = f.copy()
    i = 1         # Índice para recorrer el array de constantes zetas
    l = 128       # Longitud inicial de cada mitad de la mariposa

    while l >= 2:
        # Procesa grupos de 256/2l elementos
        for start in range(0, 256, 2 * l):
            zeta = zetas[i]
            i = i + 1
            for j in range(start, start + l):
                # Mariposa de Cooley–Tukey
                t = (zeta * f_gorro[j + l]) % q
                f_gorro[j + l] = (f_gorro[j] - t) % q
                f_gorro[j] = (f_gorro[j] + t) % q

        l = l // 2

    return f_gorro


def INTT(f_gorro):
    """
    Aplica la transformada inversa NTT^{-1} al polinomio f_gorro en T_q.

    Esta función revierte la transformación NTT, según el Algoritmo 10 del estándar FIPS 203.
    Utiliza las constantes zeta en orden inverso (BitRev7(i) con i de 127 a 1) y una multiplicación
    final por 3303 ≡ 256^{-1} mod q.

    Entrada:
    - f_gorro: lista de 256 enteros módulo q (polinomio en T_q).

    Salida:
    - Lista de 256 enteros módulo q (polinomio original en R_q).
    """
    assert(len(f_gorro) == 256)

    f = f_gorro.copy()
    i = 127       # Índice descendente para acceder a zetas en orden inverso
    l = 2         # Longitud inicial de cada mitad de la mariposa

    while l <= 128:
        # Procesa grupos de 256/2l elementos
        for start in range(0, 256, 2 * l):
            zeta = zetas[i]
            i = i - 1
            for j in range(start, start + l):
                # Mariposa inversa Gentleman–Sande
                t = f[j]
                f[j] = (t + f[j + l]) % q
                f[j + l] = (zeta * (f[j + l] - t)) % q

        l = 2 * l

    # Normalización final: multiplicar por n^{-1} mod q = 3303
    return [(x * 3303) % q for x in f]


def BaseCaseMultiply(a0, a1, b0, b1, gamma):
    """
    Multiplica dos polinomios lineales módulo X² - gamma, como en Algoritmo 12 (BaseCaseMultiply).

    Entrada:
    - a0, a1: coeficientes del primer polinomio (a0 + a1·X)
    - b0, b1: coeficientes del segundo polinomio (b0 + b1·X)
    - gamma: valor de zeta^(2·BitRev7(i)) + 1, precomputado

    Salida:
    - c0, c1: coeficientes del producto módulo X² - gamma
    """
    # Producto de polinomios módulo X² - gamma
    c0 = (a0 * b0 + a1 * b1 * gamma) % q
    c1 = (a0 * b1 + a1 * b0) % q

    return c0, c1


def MultiplyNTTs(f_gorro, g_gorro):
    """
    Multiplica dos elementos en el dominio NTT (en T_q), como en el Algoritmo 11 (MultiplyNTTs).

    Cada elemento de T_q se representa mediante 128 pares de coeficientes (grado ≤1).
    Se realiza una multiplicación por coordenadas, módulo X² - gamma_i, con gamma_i = zeta^{2·BitRev7(i) + 1}.

    Entrada:
    - f_gorro: lista de 256 enteros módulo q (representación NTT de f ∈ T_q)
    - g_gorro: lista de 256 enteros módulo q (representación NTT de g ∈ T_q)

    Salida:
    - h_gorro: lista de 256 enteros módulo q (representación NTT del producto h = f ×_T_q g)
    """
    assert(len(f_gorro) == 256 and len(g_gorro) == 256)

    h_gorro = [0 for _ in range(256)]
    for i in range(128):
        # Multiplica los polinomios (f[2i] + f[2i+1]·X) y (g[2i] + g[2i+1]·X) módulo X² - gamma
        (h_gorro[2 * i], h_gorro[2 * i + 1]) = BaseCaseMultiply(
            f_gorro[2 * i], f_gorro[2 * i + 1],
            g_gorro[2 * i], g_gorro[2 * i + 1],
            zetas_2[i]  # gamma = zeta^{2·BitRev7(i)} + 1
        )

    return h_gorro


def SumNTTs(f_gorro, g_gorro):
    """
    Suma dos elementos del dominio NTT (en T_q), componente a componente.

    Entrada:
    - f_gorro, g_gorro: listas de 256 coeficientes módulo q (representación NTT)

    Salida:
    - Lista de 256 coeficientes módulo q, resultado de f_gorro + g_gorro en T_q
    """
    assert(len(f_gorro) == 256 and len(g_gorro) == 256)

    # Suma componente a componente módulo q
    return list(map(lambda x, y: (x + y) % q, f_gorro, g_gorro))


def SubtractNTTs(f_gorro, g_gorro):
    """
    Resta dos elementos del dominio NTT (en T_q), componente a componente.

    Entrada:
    - f_gorro, g_gorro: listas de 256 coeficientes módulo q (representación NTT)

    Salida:
    - Lista de 256 coeficientes módulo q, resultado de f_gorro - g_gorro en T_q
    """
    assert(len(f_gorro) == 256 and len(g_gorro) == 256)

    # Resta componente a componente módulo q
    return list(map(lambda x, y: (x - y) % q, f_gorro, g_gorro))


def NTT_vector_vector_multiply(f_gorro, g_gorro):
    """
    Realiza el producto escalar de dos vectores cuyos elementos están en el dominio NTT (T_q).

    Entrada:
    - f_gorro, g_gorro: listas de la misma longitud, cada una conteniendo elementos de T_q 
      (cada uno representado por 256 coeficientes)

    Salida:
    - Elemento de T_q (256 coeficientes) resultado del producto escalar ∑ f[i] × g[i]
    """
    assert(len(f_gorro) == len(g_gorro))

    # Multiplica componente a componente y acumula la suma de los productos
    return reduce(SumNTTs, [MultiplyNTTs(f_gorro[i], g_gorro[i]) for i in range(len(f_gorro))])


def NTT_matrix_vector_multiply(A_gorro, s_gorro):
    """
    Realiza la multiplicación de una matriz por un vector en el dominio NTT (T_q).

    Entrada:
    - A_gorro: matriz de elementos de T_q (cada uno con 256 coeficientes)
    - s_gorro: vector de elementos de T_q (cada uno con 256 coeficientes)

    Salida:
    - Vector de elementos de T_q resultado de A_gorro × s_gorro, cada uno con 256 coeficientes
    """
    # Aplica producto escalar fila por fila
    return [NTT_vector_vector_multiply(A_gorro[i], s_gorro) for i in range(len(A_gorro))]
