import math
from conversions import BytesToBits, BitsToBytes, b2h
    
class Keccak_p:
    
    def __init__(self, b, nr):
        """
        Inicializa una instancia de la permutación Keccak-p[b, nr].

        Entrada:
        - b: tamaño del estado en bits (uno de {25, 50, 100, 200, 400, 800, 1600})
        - nr: número de rondas que se aplicarán

        Se calculan los parámetros derivados w (ancho del plano), l (log2(w)) y se guarda el número de rondas.
        """
        assert(b in [25, 50, 100, 200, 400, 800, 1600])
        self.__b = b
        self.__w = b // 25
        self.__l = int(math.log(self.__w, 2))
        self.__nr = nr
        
    def __theta(self, A):
        """
        Aplica la transformación θ al estado A.

        Entrada:
        - A: estado tridimensional de tamaño 5×5×w

        Salida:
        - A': estado tras aplicar θ
        """
        A_prime = [[[0 for _ in range (self.__w)] for _ in range(5)] for _ in range(5)]
        C = [[0 for _ in range (self.__w)] for _ in range(5)]
        D = [[0 for _ in range (self.__w)] for _ in range(5)]

        # Calcula las paridades por columna
        for x in range(5):
            for z in range(self.__w):
                C[x][z] = A[x][0][z] ^ A[x][1][z] ^ A[x][2][z] ^ A[x][3][z] ^ A[x][4][z]
        
        # Aplica rotaciones y XORs para cada columna
        for x in range(5):
            for z in range(self.__w):
                D[x][z] = C[(x - 1) % 5][z] ^ C[(x + 1) % 5][(z - 1) % self.__w]
        
        # XOR con el valor D calculado
        for x in range(5):
            for y in range(5):
                for z in range(self.__w):
                    A_prime[x][y][z] = A[x][y][z] ^ D[x][z]
                    
        return A_prime

    def __rho(self, A):
        """
        Aplica la transformación ρ (rotación de bits en el eje z).

        Entrada:
        - A: estado tridimensional 5×5×w

        Salida:
        - A': estado tras aplicar ρ
        """
        A_prime = [[[0 for _ in range (self.__w)] for _ in range(5)] for _ in range(5)]

        # La celda (0,0) no se rota
        for z in range(self.__w):
            A_prime[0][0][z] = A[0][0][z]
        
        # Aplica rotaciones sucesivas
        (x, y) = (1, 0)
        for t in range(24):
            for z in range(self.__w):
                A_prime[x][y][z] = A[x][y][(z - (t + 1)*(t + 2) // 2) % self.__w]
            (x, y) = (y, (2*x + 3*y) % 5)
                
        return A_prime

    def __pi(self, A):
        """
        Aplica la transformación π (permuta las coordenadas x, y del estado).

        Entrada:
        - A: estado tridimensional 5×5×w

        Salida:
        - A': estado tras aplicar π
        """
        A_prime = [[[0 for _ in range (self.__w)] for _ in range(5)] for _ in range(5)]

        # Aplica la permutación de coordenadas
        for x in range(5):
            for y in range(5):
                for z in range(self.__w):
                    A_prime[x][y][z] = A[(x + 3*y) % 5][x][z]
        
        return A_prime

    def __chi(self, A):
        """
        Aplica la transformación χ (no linealidad).

        Entrada:
        - A: estado tridimensional 5×5×w

        Salida:
        - A': estado tras aplicar χ
        """
        A_prime = [[[0 for _ in range (self.__w)] for _ in range(5)] for _ in range(5)]

        # XOR con AND de valores de la fila (función no lineal)
        for x in range(5):
            for y in range(5):
                for z in range(self.__w):
                    A_prime[x][y][z] = A[x][y][z] ^ ((A[(x + 1) % 5][y][z] ^ 1) & A[(x + 2) % 5][y][z])
                    
        return A_prime

    def __rc(self, t):
        """
        Devuelve el bit t-ésimo del polinomio LFSR usado en iota.

        Entrada:
        - t: índice del bit deseado

        Salida:
        - Bit (0 o 1) correspondiente al paso t del generador LFSR
        """
        if t % 255 == 0:
            return 1

        R = [1, 0, 0, 0, 0, 0, 0, 0]
        for i in range(1, (t % 255) + 1):
            R = [0] + R
            R[0] = R[0] ^ R[8]
            R[4] = R[4] ^ R[8]
            R[5] = R[5] ^ R[8]
            R[6] = R[6] ^ R[8]
            R = R[:8]
            
        return R[0]

    def __iota(self, A, ir):
        """
        Aplica la transformación ι (mezcla constante de ronda).

        Entrada:
        - A: estado tridimensional 5×5×w
        - ir: índice de ronda actual

        Salida:
        - A': estado tras mezclar con constante de ronda RC
        """
        A_prime = [[[A[x][y][z] for z in range (self.__w)] for y in range(5)] for x in range(5)]

        # Calcula la constante de ronda RC
        RC = [0 for _ in range(self.__w)]
        for j in range(self.__l + 1):
            RC[2 ** j - 1] = self.__rc(j + 7*ir)

        # Aplica la constante al bit (0,0)
        for z in range(self.__w):
            A_prime[0][0][z] = A_prime[0][0][z] ^ RC[z]
            
        return A_prime

    def __Rnd(self, A, ir):
        """
        Aplica una ronda completa del Keccak-p: θ → ρ → π → χ → ι

        Entrada:
        - A: estado actual
        - ir: número de ronda

        Salida:
        - A': estado tras una ronda de transformación
        """
        return self.__iota(self.__chi(self.__pi(self.__rho(self.__theta(A)))), ir)

    def __string_to_state(self, S):
        """
        Convierte una lista plana de bits en un estado tridimensional 5×5×w.

        Entrada:
        - S: lista de bits de tamaño b = 25 × w

        Salida:
        - A: estado tridimensional
        """
        A = [[[S[self.__w * (5 * y + x) + z] for z in range (self.__w)] for y in range(5)] for x in range(5)]
        return A

    def __state_to_string(self, A):
        """
        Convierte un estado tridimensional en una lista plana de bits.

        Entrada:
        - A: estado tridimensional 5×5×w

        Salida:
        - S: lista de bits linealizada
        """
        S = []
        for j in range(5):
            for i in range(5):
                S = S + A[i][j]
                
        return S

    def keccak(self, S):
        """
        Ejecuta la permutación Keccak-p sobre la entrada S.

        Entrada:
        - S: lista de bits de tamaño b

        Salida:
        - S': lista de bits tras aplicar nr rondas de Keccak-p
        """
        A = self.__string_to_state(S)
        
        # Aplica las rondas especificadas
        for ir in range(12 + 2*self.__l - self.__nr, 12 + 2*self.__l):
            A = self.__Rnd(A, ir)
        
        S_prime = self.__state_to_string(A)
        
        return S_prime
    
class Keccak_f:
    """
    Clase que representa la permutación Keccak-f[b], utilizada como núcleo en SHA-3.
    
    Esta versión fija el número de rondas como 12 + 2l, donde l = log2(w), y w = b / 25,
    tal como se especifica en la familia Keccak-f.
    """
    
    def __init__(self, b):
        """
        Inicializa la permutación Keccak-f[b] a partir del valor de b.

        Entrada:
        - b: tamaño del estado en bits (debe ser uno de los valores válidos en Keccak: 25, 50, ..., 1600)

        Internamente, se crea una instancia de Keccak-p con el número de rondas 12 + 2 * log2(b / 25).
        """
        self.__keccak = Keccak_p(b, 12 + 2*int(math.log(b // 25, 2)))
        
    def keccak(self, S):
        """
        Aplica la permutación Keccak-f[b] sobre la cadena de bits de entrada.

        Entrada:
        - S: lista de bits (tamaño b)

        Salida:
        - Lista de bits resultante tras aplicar la permutación Keccak-f[b]
        """
        return self.__keccak.keccak(S)
    
class Sponge:
    """
    Implementa la construcción de esponja (sponge construction), utilizada en funciones hash y KDFs.

    Entrada:
    - f: función de permutación sobre bloques de tamaño b
    - pad: función de padding dependiente de r y del tamaño del mensaje
    - r: tasa de absorción
    - b: tamaño total del estado interno (b = r + c)
    """
    def __init__(self, f, pad, r, b):
        self.__f = f            # Se guarda la función de permutación
        self.__pad = pad        # Se guarda la función de padding
        self.__r = r            # Tasa de absorción (bytes procesados por iteración)
        self.__b = b            # Tamaño total del estado interno
        self.__c = b - r        # Capacidad (parte oculta del estado)
        self.__S = [0] * b      # Estado interno inicializado a ceros
        self.__pos = 0          # Posición actual dentro de la fase de extracción

    def absorb(self, N):
        """
        Absorbe los datos de entrada N en el estado interno del esponjado.
        
        Entrada:
        - N: lista de enteros (bytes) a absorber
        """
        # Se aplica el padding a la entrada para que su longitud sea múltiplo de r
        P = N + self.__pad(self.__r, len(N))
        
        # Se divide el mensaje en bloques de r bytes
        n = len(P) // self.__r

        for i in range(n):
            # Se construye un bloque de tamaño b = r + c rellenando con ceros la parte oculta
            P_padded = P[self.__r * i : self.__r * (i + 1)] + [0] * self.__c
            
            # Se mezcla el bloque con el estado actual mediante XOR
            S_xor = [self.__S[j] ^ P_padded[j] for j in range(self.__b)]
            
            # Se aplica la función de permutación a todo el estado
            self.__S = self.__f(S_xor)

    def squeeze(self, d):
        """
        Extrae d bytes del estado interno en modo de goteo (squeeze).
        
        Entrada:
        - d: número de bytes de salida a generar
        
        Salida:
        - Lista de d bytes extraídos del estado
        """
        assert d >= 0  # La cantidad de bytes a extraer debe ser no negativa
        
        # Se extraen los bytes disponibles en la tasa desde la posición actual
        Z = self.__S[self.__pos : self.__r]

        # Mientras no se tengan suficientes bytes, se aplica f para generar más
        while d > len(Z):
            self.__S = self.__f(self.__S)     # Nueva permutación del estado
            Z += self.__S[:self.__r]          # Se añaden los primeros r bytes del nuevo estado

        # Se actualiza la posición para la siguiente extracción parcial (si se repite squeeze)
        self.__pos = self.__r - (len(Z) - d)
        
        # Se devuelve exactamente d bytes de salida
        return Z[:d]
    
class SHA_3_Keccak:
    """
    Implementación general de SHA-3 usando la construcción Keccak-f[1600] con padding pad10*1.
    """

    def __init__(self, c):
        """
        Inicializa el objeto SHA_3_Keccak con una capacidad dada.

        Entrada:
        - c: capacidad del algoritmo SHA-3 en bits (por ejemplo: 448, 512, 768, 1024)
        """
        # Se define el objeto esponja con Keccak-f, padding pad10*1, tasa r = 1600 - c, y estado b = 1600
        self.__sponge = Sponge(Keccak_f(1600).keccak, self.__pad101, 1600 - c, 1600)

    def __pad101(self, x, m):
        """
        Implementa el padding multi-rate pad10*1 de Keccak.

        Entrada:
        - x: tasa del algoritmo
        - m: longitud del mensaje en bits antes del padding

        Salida:
        - Lista de bits correspondiente al padding necesario para completar un múltiplo de x
        """
        assert(x > 0)
        assert(m >= 0)
        
        j = (-m - 2) % x
        return [1] + [0 for _ in range(j)] + [1]
    
    def keccak(self, N, d):
        """
        Aplica el algoritmo Keccak al mensaje N para generar una salida de d bits.

        Entrada:
        - N: lista de bits (mensaje de entrada)
        - d: longitud deseada de la salida en bits

        Salida:
        - Lista de bits con la salida del hash
        """
        self.__sponge.absorb(N + [0, 1])
        return self.__sponge.squeeze(d)


class SHA_3:
    """
    Implementación de las variantes SHA3-224, SHA3-256, SHA3-384 y SHA3-512 usando Keccak-f[1600].
    """

    def __init__(self):
        """
        Inicializa los objetos SHA3 para las distintas capacidades:
        - SHA3-224 (c = 448)
        - SHA3-256 (c = 512)
        - SHA3-384 (c = 768)
        - SHA3-512 (c = 1024)
        """
        self.__sha3_224_keccak = SHA_3_Keccak(448)
        self.__sha3_256_keccak = SHA_3_Keccak(512)
        self.__sha3_384_keccak = SHA_3_Keccak(768)
        self.__sha3_512_keccak = SHA_3_Keccak(1024)

    def sha_3_224(self, M):
        """
        Calcula SHA3-224 sobre el mensaje M dado como lista de bytes.
        Devuelve el digest como lista de bytes (28 bytes).
        """
        return BitsToBytes(self.__sha3_224_keccak.keccak(BytesToBits(M), 224))
    
    def sha_3_256(self, M):
        """
        Calcula SHA3-256 sobre el mensaje M dado como lista de bytes.
        Devuelve el digest como lista de bytes (32 bytes).
        """
        return BitsToBytes(self.__sha3_256_keccak.keccak(BytesToBits(M), 256))
    
    def sha_3_384(self, M):
        """
        Calcula SHA3-384 sobre el mensaje M dado como lista de bytes.
        Devuelve el digest como lista de bytes (48 bytes).
        """
        return BitsToBytes(self.__sha3_384_keccak.keccak(BytesToBits(M), 384))
    
    def sha_3_512(self, M):
        """
        Calcula SHA3-512 sobre el mensaje M dado como lista de bytes.
        Devuelve el digest como lista de bytes (64 bytes).
        """
        return BitsToBytes(self.__sha3_512_keccak.keccak(BytesToBits(M), 512))
    
class SHAKE_Keccak:
    """
    Implementación general de SHAKE usando Keccak-f[1600].
    """

    def __init__(self, c):
        """
        Inicializa el objeto SHAKE_Keccak con una capacidad dada.

        Entrada:
        - c: capacidad del algoritmo SHAKE (por ejemplo: 256, 512)
        """
        # Se define el objeto esponja con Keccak-f, padding pad10*1, tasa r = 1600 - c, y estado b = 1600
        self.__sponge = Sponge(Keccak_f(1600).keccak, self.__pad101, 1600 - c, 1600)

    def __pad101(self, x, m):
        """
        Implementa el padding multi-rate pad10*1 usado en SHAKE.

        Entrada:
        - x: tasa del algoritmo (en bits)
        - m: longitud del mensaje original (en bits)

        Salida:
        - Lista de bits con el padding adecuado
        """
        assert(x > 0)
        assert(m >= 0)

        j = (-m - 2) % x
        return [1] + [0 for _ in range(j)] + [1]

    def absorb(self, N):
        """
        Absorbe el mensaje de entrada N en la esponja, incluyendo los bits de dominio para SHAKE.

        Entrada:
        - N: lista de bits correspondiente al mensaje de entrada
        """
        # Se añaden los bits de dominio [1,1,1,1] para SHAKE antes de la absorción
        self.__sponge.absorb(N + [1, 1, 1, 1])

    def squeeze(self, d):
        """
        Extrae d bits de salida del estado esponjado.

        Entrada:
        - d: número de bits deseado en la salida

        Salida:
        - Lista de bits con la salida de longitud d
        """
        return self.__sponge.squeeze(d)


class SHAKE:
    """
    Implementación de las variantes SHAKE128 y SHAKE256 usando Keccak-f[1600].
    """

    def __init__(self):
        """
        Inicializa los objetos SHAKE128 y SHAKE256 con sus capacidades correspondientes.
        """
        self.__shake128_keccak = SHAKE_Keccak(256)
        self.__shake256_keccak = SHAKE_Keccak(512)

    def shake128(self, M, d):
        """
        Calcula SHAKE128 sobre el mensaje M con salida de d bits.

        Entrada:
        - M: lista de bytes (mensaje de entrada)
        - d: número de bits deseado en la salida

        Salida:
        - Lista de bytes de longitud d // 8
        """
        self.__shake128_keccak.absorb(BytesToBits(M))
        return BitsToBytes(self.__shake128_keccak.squeeze(d))

    def shake256(self, M, d):
        """
        Calcula SHAKE256 sobre el mensaje M con salida de d bits.

        Entrada:
        - M: lista de bytes (mensaje de entrada)
        - d: número de bits deseado en la salida

        Salida:
        - Lista de bytes de longitud d // 8
        """
        self.__shake256_keccak.absorb(BytesToBits(M))
        return BitsToBytes(self.__shake256_keccak.squeeze(d))
    
class XOF:
    """
    Clase XOF (eXtendable Output Function) basada en SHAKE128 con capacidad para absorber
    y extraer (squeeze) bits de longitud variable.

    Métodos:
    - absorb(N): absorbe la entrada N en forma de bits.
    - squeeze(l): extrae l bytes de salida pseudoaleatoria.
    """

    def __init__(self):
        """
        Inicializa un objeto SHAKE_Keccak con capacidad 256 (SHAKE128).
        """
        self.__shake128_keccak = SHAKE_Keccak(256)
        
    def absorb(self, N):
        """
        Absorbe la entrada N en el estado interno del SHAKE en forma de bits.

        Entrada:
        - N: lista de bytes
        """
        self.__shake128_keccak.absorb(BytesToBits(N))
    
    def squeeze(self, l):
        """
        Extrae l bytes de salida pseudoaleatoria a partir del estado interno.

        Entrada:
        - l: número entero, cantidad de bytes a extraer.

        Salida:
        - Lista de bytes de longitud l.
        """
        return BitsToBytes(self.__shake128_keccak.squeeze(8 * l))
    
def PRF(eta, s, b):
    """
    Función pseudoaleatoria determinista (PRF) parametrizada para el esquema.

    Genera una salida SHAKE256 con entrada la concatenación de s y un byte b,
    y longitud 8 * 64 * eta bits, con eta ∈ {2, 3}.

    Entrada:
    - eta: entero 2 o 3, parámetro del esquema.
    - s: lista de 32 bytes (clave/secreto).
    - b: entero entre 0 y 255 (byte).

    Salida:
    - Lista de bytes con la salida de la función PRF.
    """
    assert(eta == 2 or eta == 3)
    assert(len(s)  == 32)
    assert(0 <= b <= 255)
    
    return SHAKE().shake256(s + [b], 8 * 64 * eta)

def H(s):
    """
    Función hash H basada en SHA3-256.

    Entrada:
    - s: lista de bytes (mensaje).

    Salida:
    - Hash SHA3-256 de s, lista de bytes (32 bytes).
    """
    return SHA_3().sha_3_256(s)

def J(s):
    """
    Función hash extendida J basada en SHAKE256.

    Entrada:
    - s: lista de bytes (mensaje).

    Salida:
    - Hash SHAKE256 de longitud 256 bits (32 bytes).
    """
    return SHAKE().shake256(s, 8 * 32)

def G(c):
    """
    Función hash G basada en SHA3-512.

    Entrada:
    - c: lista de bytes (mensaje).

    Salida:
    - Tupla de dos elementos, cada uno con 32 bytes, que son las dos mitades
      de la salida SHA3-512 (64 bytes).
    """
    g = SHA_3().sha_3_512(c)
    return g[:32], g[32:]