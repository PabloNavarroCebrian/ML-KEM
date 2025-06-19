from keccak import G, PRF
from sampling import SampleNTT, SamplePolyCBD
from ntt import NTT, INTT, NTT_matrix_vector_multiply, SumNTTs, NTT_vector_vector_multiply, SubtractNTTs
from conversions import ByteEncode, ByteDecode, transpose, Compress, Decompress
from functools import reduce

class K_PKE:
    
    def __init__(self, k, eta1, eta2, du, dv):
        """
        Inicializa una instancia del esquema K-PKE.

        Entrada:
        - k (int): número de filas/columnas de la matriz A (dimensión del esquema).
        - eta1 (int): parámetro de ruido para las distribuciones de s y e.
        - eta2 (int): parámetro de ruido para las distribuciones de e1 y e2.
        - du (int): parámetro de compresión para el componente c1 del cifrado.
        - dv (int): parámetro de compresión para el componente c2 del cifrado.
        """
        self.__k = k
        self.__eta1 = eta1
        self.__eta2 = eta2
        self.__du = du
        self.__dv = dv
    
    def KeyGen(self, d):
        """
        Genera un par de claves (pública y secreta) para el esquema K-PKE.

        Entrada:
        - d (list[int]): semilla de entrada de tamaño 32 bytes.

        Salida:
        - ek_PKE (list[int]): clave pública comprimida.
        - dk_PKE (list[int]): clave secreta comprimida (en dominio NTT).
        """
        (rho, sigma) = G(d + [self.__k])  # Expansión determinista de la semilla en rho y sigma
        N = 0  # Contador para la función PRF
        
        # Construcción de la matriz pública A ∈ R_q^{k×k} en dominio NTT
        A = [[[0 for _ in range(256)] for _ in range(self.__k)] for _ in range(self.__k)]
        for i in range(self.__k):
            for j in range(self.__k):
                A[i][j] = SampleNTT(rho + [j, i])  # A[i][j] = XOF(rho || j || i)
        
        # Generación del vector secreto s ∈ R_q^k usando CBD con semilla sigma
        s = [[0 for _ in range(256)] for _ in range(self.__k)]
        for i in range(self.__k):
            s[i] = SamplePolyCBD(self.__eta1, PRF(self.__eta1, sigma, N))
            N = N + 1
        
        # Generación del vector de errores e ∈ R_q^k
        e = [[0 for _ in range(256)] for _ in range(self.__k)]
        for i in range(self.__k):
            e[i] = SamplePolyCBD(self.__eta1, PRF(self.__eta1, sigma, N))
            N = N + 1
        
        # Transformación NTT de s y e
        s_gorro = list(map(NTT, s))
        e_gorro = list(map(NTT, e))
        
        # Cálculo de t̂ = A·s_gorro + e_gorro
        t_gorro = list(map(SumNTTs, NTT_matrix_vector_multiply(A, s_gorro), e_gorro))
        
        # Codificación de la clave pública: incluye t_gorro y rho
        ek_PKE = []
        for i in range(self.__k):
            ek_PKE = ek_PKE + ByteEncode(12, t_gorro[i])
        ek_PKE = ek_PKE + rho
        
        # Codificación de la clave secreta: solo s_gorro
        dk_PKE = []
        for i in range(self.__k):
            dk_PKE = dk_PKE + ByteEncode(12, s_gorro[i])
            
        return ek_PKE, dk_PKE
    
    def Encrypt(self, ek_PKE, m, r):
        """
        Cifra un mensaje m utilizando la clave pública y una semilla aleatoria.

        Entrada:
        - ek_PKE (list[int]): clave pública.
        - m (list[int]): mensaje de 32 bytes a cifrar.
        - r (list[int]): semilla aleatoria para la generación de ruido.

        Salida:
        - c (list[int]): cifrado (c1 || c2).
        """
        N = 0  # Contador para PRF
        
        # Decodificación de t̂ a partir de ek_PKE
        t_gorro = []
        for i in range(self.__k):
            t_gorro.append(ByteDecode(12, ek_PKE[384 * i : 384 * (i + 1)]))
        rho = ek_PKE[384 * self.__k:]  # Extracción de la semilla rho
        
        # Reconstrucción de la matriz A a partir de rho
        A = [[[0 for _ in range(256)] for _ in range(self.__k)] for _ in range(self.__k)]
        for i in range(self.__k):
            for j in range(self.__k):
                A[i][j] = SampleNTT(rho + [j, i])
        
        # Generación del vector aleatorio y ∈ R_q^k
        y = [[0 for _ in range(256)] for _ in range(self.__k)]
        for i in range(self.__k):
            y[i] = SamplePolyCBD(self.__eta1, PRF(self.__eta1, r, N))
            N = N + 1
        
        # Generación del vector de errores e1 ∈ R_q^k
        e1 = [[0 for _ in range(256)] for _ in range(self.__k)]
        for i in range(self.__k):
            e1[i] = SamplePolyCBD(self.__eta2, PRF(self.__eta2, r, N))
            N = N + 1
        
        # Generación del error e2 ∈ R_q
        e2 = SamplePolyCBD(self.__eta2, PRF(self.__eta2, r, N))
        
        # Transformación NTT del vector y
        y_gorro = list(map(NTT, y))
        
        # Cálculo de u = INTT(Aᵗ·y_gorro) + e1
        u = list(map(SumNTTs, list(map(INTT, NTT_matrix_vector_multiply(transpose(A), y_gorro))), e1))
        
        # Transformación del mensaje m a mu (0 --> 0 y 1 --> floor(q/2))
        mu = [Decompress(1, x) for x in ByteDecode(1, m)]
        
        # Cálculo de v = INTT(t_gorro·_gorroy) + e2 + μ
        v = reduce(SumNTTs, [INTT(NTT_vector_vector_multiply(t_gorro, y_gorro)), e2, mu])
        
        # Codificación del componente c1: compresión de u
        c1 = []
        for i in range(self.__k):
            c1 = c1 + ByteEncode(self.__du, [Compress(self.__du, x) for x in u[i]])
        
        # Codificación del componente c2: compresión de v
        c2 = ByteEncode(self.__dv, [Compress(self.__dv, x) for x in v])
        
        return c1 + c2
    
    def Decrypt(self, dk_PKE, c):
        """
        Descifra un cifrado c utilizando la clave secreta.

        Entrada:
        - dk_PKE (list[int]): clave secreta.
        - c (list[int]): cifrado (c1 || c2).

        Salida:
        - m (list[int]): mensaje descifrado como lista de 32 bytes.
        """
        # Separación del cifrado en componentes c1 y c2
        c1 = c[:32 * self.__du * self.__k]
        c2 = c[32 * self.__du * self.__k:]
        
        # Reconstrucción de u' a partir de c1
        u_prime = []
        for i in range(self.__k):
            u_prime.append([Decompress(self.__du, x) for x in ByteDecode(self.__du, c1[32 * self.__du * i: 32 * self.__du * (i + 1)])])
        
        # Reconstrucción de v' a partir de c2
        v_prime = [Decompress(self.__dv, x) for x in ByteDecode(self.__dv, c2)]
        
        # Decodificación de s_gorro desde la clave secreta
        s_gorro = []
        for i in range(self.__k):
            s_gorro.append(ByteDecode(12, dk_PKE[384 * i: 384 * (i + 1)]))
        
        # Cálculo de w = v' - INTT(s_gorro·NTT(u'))
        w = SubtractNTTs(v_prime, INTT(NTT_vector_vector_multiply(s_gorro, list(map(NTT, u_prime)))))
        
        # Decodificación del mensaje final m
        m = ByteEncode(1, [Compress(1, x) for x in w])
        
        return m