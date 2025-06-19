from keccak import XOF
from conversions import BytesToBits

q = 3329

def SampleNTT(B):
    """
    Realiza el muestreo uniforme de una representación en el dominio NTT.

    Esta función implementa el algoritmo SampleNTT (Algoritmo 7 del estándar FIPS 203),
    que convierte una semilla de 32 bytes junto con dos bytes de índice (total: 34 bytes)
    en un polinomio de grado 256 con coeficientes en ℤ_q, directamente en el dominio NTT.

    Si la semilla de entrada es uniformemente aleatoria, entonces el polinomio resultante
    es indistinguible computacionalmente de uno tomado de la distribución uniforme sobre T_q
    (conjunto de representaciones NTT válidas). Esto se logra usando una función XOF para 
    generar valores pseudoaleatorios que son mapeados a enteros < q.

    Entrada:
    - B: lista de 34 bytes (32 de semilla y 2 de índice).

    Salida:
    - a: lista de 256 enteros en ℤ_q que representan un polinomio muestreado uniformemente en T_q.
    """
    assert(len(B) == 34)

    xof = XOF()
    xof.absorb(B)

    a = [0 for _ in range(256)]
    j = 0
    while j < 256:
        # Extraemos 3 bytes para producir hasta dos coeficientes por iteración
        C = xof.squeeze(3)

        # Primer candidato: usa los primeros 8 bits y 4 bits del segundo byte
        d1 = C[0] + 256 * (C[1] % 16)

        # Segundo candidato: usa 4 bits del segundo byte y 8 bits del tercero
        d2 = C[1] // 16 + 16 * C[2]

        if d1 < q:
            a[j] = d1
            j += 1

        if d2 < q and j < 256:
            a[j] = d2
            j += 1

    return a


def SamplePolyCBD(eta, B):
    """
    Muestra un polinomio con coeficientes pequeños según la distribución binomial centrada D_η(ℤ_q).

    Esta función implementa el algoritmo SamplePolyCBD (Algoritmo 8 del estándar FIPS 203),
    que genera un polinomio en R_q con coeficientes distribuidos según una distribución binomial
    centrada con parámetro η ∈ {2, 3}. Esta distribución se usa para generar los polinomios de error ("noise").

    Para cada uno de los 256 coeficientes, se generan 2η bits: η de ellos se suman para obtener x, y los otros η
    se suman para obtener y. El coeficiente es entonces (x - y) mod q, dando un valor en ℤ_q centrado alrededor de 0.

    Entrada:
    - eta: parámetro de la distribución binomial centrada, debe ser 2 o 3.
    - B: lista de 64 * eta bytes aleatorios (uniformemente distribuidos).

    Salida:
    - f: lista de 256 enteros en ℤ_q que forman el polinomio muestreado.
    """
    assert(eta == 2 or eta == 3)
    assert(len(B) == 64 * eta)

    # Convertimos los bytes a una lista de bits (little-endian por byte)
    b = BytesToBits(B)

    f = [0 for _ in range(256)]
    for i in range(256):
        # Sumamos η bits para x y otros η bits para y
        x = sum(b[2 * i * eta + j] for j in range(eta))
        y = sum(b[2 * i * eta + eta + j] for j in range(eta))

        # La diferencia centrada (x - y) se reduce módulo q
        f[i] = (x - y) % q

    return f
