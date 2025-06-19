import math

q = 3329
def h2b(H, n=None):
    """
    Convierte un array hexadecimal H en un array de bits de longitud n,
    extrayendo los bits en orden b0 (LSB) a b7 (MSB) por cada byte.
    
    Entrada: 
    - H: Array hexadecimal (2m caracteres hexadecimales)
    - n: Número de bits a devolver (n ≤ 8m)
    
    Salida:
    - Array de bits de longitud n
    """
    # Paso 1: Convertir cada carácter hexadecimal a su valor numérico
    hex_values = [int(ch, 16) for ch in H]

    m = len(H) // 2  # Porque H tiene 2m dígitos hexadecimales

    # Paso 2: Agrupar cada par de dígitos hexadecimales en un byte (valor entre 0 y 255)
    bytes_list = []
    for i in range(m):
        hi = 16 * hex_values[2 * i] + hex_values[2 * i + 1]
        bytes_list.append(hi)

    # Paso 3: Convertir cada byte a su representación binaria de 8 bits
    # PERO en orden b0 (LSB) a b7 (MSB)
    bit_list = []
    for hi in bytes_list:
        bits = [((hi >> bit) & 1) for bit in range(8)]  # De b0 (LSB) a b7 (MSB)
        bit_list.extend(bits)

    if n != None:
        # Paso 4: Truncar la lista de bits a n bits y convertir a string
        S = bit_list[:n]
    else:
        S = bit_list

    return S

def b2h(S):
    """
    Convierte un array de bits S en un array hexadecimal H,
    considerando que S está en orden b0 (LSB) a bN (MSB).
    
    Entrada:
    - S: Array de bits ('0' y '1')
    
    Salida:
    - Cadena hexadecimal H
    """
    n = len(S)

    # Paso 1: Padding con ceros al final si es necesario
    padding = (8 - (n % 8)) % 8
    T = S + [0] * padding
    m = math.ceil(n / 8)

    # Paso 2 y 3: Procesar cada bloque de 8 bits
    H = []
    for i in range(m):
        byte_bits = T[8 * i: 8 * (i + 1)]  # Extraer los 8 bits (b0 a b7)

        # Interpretar cada bloque en el orden b0 * 2^0 + b1 * 2^1 + ... + b7 * 2^7
        hi = sum(byte_bits[j] << j for j in range(8))

        # Convertir hi en dos dígitos hexadecimales
        H2i = hi // 16
        H2i1 = hi % 16

        H += [format(H2i, 'X'), format(H2i1, 'X')]

    return H

def BytesToBits(B):
    """
    Convierte una lista de bytes (enteros entre 0 y 255) en una lista de bits.
    
    Cada byte se descompone en sus 8 bits en orden little-endian (bit menos significativo primero).
    
    Entrada: 
    - B: lista de enteros (bytes).
    
    Salida: 
    - Lista de bits (enteros 0 o 1) de longitud 8 * len(B).
    """
    # Creamos una copia de la lista de bytes para no modificar el argumento original
    C = B.copy()

    # Inicializamos una lista de bits de tamaño 8 veces el número de bytes
    b = [0 for _ in range(8 * len(B))]

    # Recorremos cada byte
    for i in range(len(B)):
        # Para cada byte, extraemos sus 8 bits (del menos significativo al más significativo)
        for j in range(8):
            # Guardamos el bit menos significativo del byte actual
            b[8 * i + j] = C[i] % 2
            # Eliminamos el bit ya procesado dividiendo entre 2
            C[i] = C[i] // 2

    # Devolvemos la lista de bits
    return b


def BitsToBytes(b):
    """
    Convierte una lista de bits (enteros 0 o 1) en una lista de bytes (enteros entre 0 y 255).

    La entrada debe tener una longitud múltiplo de 8. Los bits se agrupan de 8 en 8 en orden little-endian
    (el bit menos significativo primero), formando un byte por cada grupo.
    
    Entrada: 
    - b: lista de bits (enteros 0 o 1).
    
    Salida: 
    - Lista de enteros (bytes) de longitud len(b) // 8.
    """
    # Comprobamos que la longitud de la lista sea múltiplo de 8
    assert((len(b) % 8) == 0)

    # Inicializamos la lista de bytes con ceros
    B = [0 for _ in range(len(b) // 8)]

    # Recorremos la lista de bits
    for i in range(len(b)):
        # Calculamos el valor del byte correspondiente acumulando potencias de 2
        B[i // 8] += b[i] * (2 ** (i % 8))

    # Devolvemos la lista de bytes
    return B


def ByteEncode(d, F):
    """
    Codifica una lista de 256 enteros F, cada uno en el rango [0, m) con m = 2^d si d < 12 o m = q si d = 12, en una secuencia compacta de bytes.

    Cada entero se representa con exactamente d bits en orden little-endian (bit menos significativo primero),
    y los 256*d bits resultantes se empaquetan en una lista de bytes.

    Entrada:
    - d: número de bits por entero (1 ≤ d ≤ 12).
    - F: lista de 256 enteros, cada uno en [0, 2^d).

    Salida:
    - Lista de bytes (enteros entre 0 y 255) que codifican los bits de F.
    """
    assert(len(F) == 256)
    assert(1 <= d <= 12)

    # Creamos una lista de bits de tamaño 256*d
    b = [0 for _ in range(256 * d)]

    # Para cada entero en F, lo convertimos en su representación binaria de d bits
    for i in range(256):
        a = F[i]
        for j in range(d):
            # Almacenamos el bit menos significativo
            b[i * d + j] = a % 2
            # Eliminamos ese bit para la siguiente iteración
            a = (a - b[i * d + j]) // 2

    # Convertimos la lista de bits en una lista de bytes
    return BitsToBytes(b)


def ByteDecode(d, B):
    """
    Decodifica una secuencia de bytes que representa 256 enteros codificados con d bits cada uno,
    y reconstruye la lista original de enteros módulo m.

    Cada grupo de d bits consecutivos en orden little-endian representa un entero.
    Si d < 12, m = 2^d; si d == 12, m = q (un valor predefinido en el contexto de Kyber).

    Entrada:
    - d: número de bits por entero (1 ≤ d ≤ 12).
    - B: lista de bytes (enteros entre 0 y 255) de longitud 32*d, que codifican 256 enteros.

    Salida:
    - F: lista de 256 enteros, cada uno en [0, m), reconstruidos a partir de los bits.
    """
    assert(len(B) == 32 * d)
    assert(1 <= d <= 12)

    # Inicializamos la lista de enteros reconstruidos
    F = [0 for _ in range(256)]

    # Definimos el módulo m según el valor de d
    m = 2 ** d if d < 12 else q

    # Convertimos los bytes en una lista de bits
    b = BytesToBits(B)

    # Para cada grupo de d bits, reconstruimos el entero correspondiente
    for i in range(256):
        F[i] = sum(b[i * d + j] << j for j in range(d)) % m

    return F


def transpose(A):
    """
    Transpone una matriz representada como lista de listas.

    Entrada:
    - A: matriz de tamaño n × m representada como lista de n listas (filas), cada una de longitud m.

    Salida:
    - Matriz transpuesta de tamaño m × n (columnas se convierten en filas).
    """
    return list(map(list, zip(*A)))


def Compress(d, x):
    """
    Comprime un entero x ∈ [0, q) a un valor representable con d bits.

    Esta función se utiliza en Kyber para reducir la precisión de los coeficientes de los polinomios
    antes de ser empaquetados. El valor x se escala al rango [0, 2^d), se redondea y se reduce módulo 2^d.

    Entrada:
    - d: número de bits de precisión objetivo (0 < d < 12).
    - x: entero en el rango [0, q).

    Salida:
    - Entero comprimido en el rango [0, 2^d).
    """
    assert(d < 12)

    # Escalamos x del rango [0, q) al rango [0, 2^d), redondeamos y reducimos módulo 2^d
    return round(((2 ** d) / q) * x) % (2 ** d)


def Decompress(d, y):
    """
    Descomprime un entero comprimido y de d bits al rango original [0, q).

    Esta función es la inversa de `Compress`. Recibe un valor representado con d bits y
    lo escala al rango completo [0, q), aproximando el valor original antes de la compresión.

    Entrada:
    - d: número de bits de precisión (0 < d < 12).
    - y: entero en el rango [0, 2^d).

    Salida:
    - Entero aproximado en el rango [0, q).
    """
    assert(d < 12)

    # Escalamos y del rango [0, 2^d) al rango [0, q), redondeando al entero más cercano
    return round((q / (2 ** d)) * y)
