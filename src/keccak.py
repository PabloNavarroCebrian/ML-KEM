import math
from conversions import BytesToBits, BitsToBytes, b2h

class Keccak_p:
    
    def __init__(self, b, nr):
        assert(b in [25, 50, 100, 200, 400, 800, 1600])
        self.__b = b
        self.__w = b // 25
        self.__l = int(math.log(self.__w, 2))
        self.__nr = nr
        
    def __theta(self, A):
        
        A_prime = [[[0 for _ in range (self.__w)] for _ in range(5)] for _ in range(5)]
        C = [[0 for _ in range (self.__w)] for _ in range(5)]
        D = [[0 for _ in range (self.__w)] for _ in range(5)]
        for x in range(5):
            for z in range(self.__w):
                C[x][z] = A[x][0][z] ^ A[x][1][z] ^ A[x][2][z] ^ A[x][3][z] ^ A[x][4][z]
                
        for x in range(5):
            for z in range(self.__w):
                D[x][z] = C[(x - 1) % 5][z] ^ C[(x + 1) % 5][(z - 1) % self.__w]
                
        for x in range(5):
            for y in range(5):
                for z in range(self.__w):
                    A_prime[x][y][z] = A[x][y][z] ^ D[x][z]
                    
        return A_prime

    def __rho(self, A):
        
        A_prime = [[[0 for _ in range (self.__w)] for _ in range(5)] for _ in range(5)]
        
        for z in range(self.__w):
            A_prime[0][0][z] = A[0][0][z]
            
        (x, y) = (1, 0)
        
        for t in range(24):
            for z in range(self.__w):
                A_prime[x][y][z] = A[x][y][(z - (t + 1)*(t + 2) // 2) % self.__w]
            (x, y) = (y, (2*x + 3*y) % 5)
                
        return A_prime

    def __pi(self, A):
        
        A_prime = [[[0 for _ in range (self.__w)] for _ in range(5)] for _ in range(5)]
        
        for x in range(5):
            for y in range(5):
                for z in range(self.__w):
                    A_prime[x][y][z] = A[(x + 3*y) % 5][x][z]
        
        return A_prime

    def __chi(self, A):
        
        A_prime = [[[0 for _ in range (self.__w)] for _ in range(5)] for _ in range(5)]
        
        for x in range(5):
            for y in range(5):
                for z in range(self.__w):
                    A_prime[x][y][z] = A[x][y][z] ^ ((A[(x + 1) % 5][y][z] ^ 1) & A[(x + 2) % 5][y][z])
                    
        return A_prime

    def __rc(self, t):
        
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
        
        A_prime = [[[A[x][y][z] for z in range (self.__w)] for y in range(5)] for x in range(5)]
        
        RC = [0 for _ in range(self.__w)]
        for j in range(self.__l + 1):
            RC[2 ** j - 1] = self.__rc(j + 7*ir)
                    
        for z in range(self.__w):
            A_prime[0][0][z] = A_prime[0][0][z] ^ RC[z]
            
        return A_prime

    def __Rnd(self, A, ir):
        return self.__iota(self.__chi(self.__pi(self.__rho(self.__theta(A)))), ir)

    def __string_to_state(self, S):
        A = [[[S[self.__w * (5 * y + x) + z] for z in range (self.__w)] for y in range(5)] for x in range(5)]
        return A

    def __state_to_string(self, A):
        
        S = []
        for j in range(5):
            for i in range(5):
                S = S + A[i][j]
                
        return S

    def keccak(self, S):
        
        A = self.__string_to_state(S)
        
        for ir in range(12 + 2*self.__l - self.__nr, 12 + 2*self.__l):
            A = self.__Rnd(A, ir)
            
        S_prime = self.__state_to_string(A)
        
        return S_prime
    
class Keccak_f:
    
    def __init__(self, b):
        self.__keccak = Keccak_p(b, 12 + 2*int(math.log(b // 25, 2)))
        
    def keccak(self, S):
        return self.__keccak.keccak(S)
        
class Sponge:
    def __init__(self, f, pad, r, b):
        self.__f = f  # Función de permutación
        self.__pad = pad  # Función de padding
        self.__r = r  # Tasa de absorción
        self.__b = b  # Tamaño del estado
        self.__c = b - r  # Capacidad
        self.__S = [0] * b  # Estado inicializado a ceros
        self.__pos = 0

    def absorb(self, N):
        """Absorbe los datos de entrada N en el estado del esponjado."""
        P = N + self.__pad(self.__r, len(N))  # Aplicar padding
        n = len(P) // self.__r  # Número de bloques

        for i in range(n):
            # Construir bloque con padding de capacidad
            P_padded = P[self.__r * i : self.__r * (i + 1)] + [0] * self.__c
            # XOR con el estado actual y aplicar la función f
            self.__S = self.__f([self.__S[j] ^ P_padded[j] for j in range(self.__b)])

    def squeeze(self, d):
        """Genera d bytes de salida en modo de goteo."""
        assert d >= 0
        
        Z = self.__S[self.__pos : self.__r]

        while d > len(Z):
            self.__S = self.__f(self.__S) # Aplicar f al estado
            Z += self.__S[:self.__r]  # Extraer más datos

        self.__pos = self.__r - (len(Z) - d)
        
        return Z[:d]
    
class SHA_3_Keccak:
    
    def __init__(self, c):
        self.__sponge = Sponge(Keccak_f(1600).keccak, self.__pad101, 1600 - c, 1600)
        
    def __pad101(self, x, m):
        assert(x > 0)
        assert(m >= 0)
        
        j = (-m - 2) % x
        return [1] + [0 for _ in range(j)] + [1]
    
    def keccak(self, N, d):
        self.__sponge.absorb(N + [0, 1])
        return self.__sponge.squeeze(d)

class SHA_3:
    def __init__(self):
        self.__sha3_224_keccak = SHA_3_Keccak(448)
        self.__sha3_256_keccak = SHA_3_Keccak(512)
        self.__sha3_384_keccak = SHA_3_Keccak(768)
        self.__sha3_512_keccak = SHA_3_Keccak(1024)
        
    def sha_3_224(self, M):
        return BitsToBytes(self.__sha3_224_keccak.keccak(BytesToBits(M), 224))
    
    def sha_3_256(self, M):
        return BitsToBytes(self.__sha3_256_keccak.keccak(BytesToBits(M), 256))
    
    def sha_3_384(self, M):
        return BitsToBytes(self.__sha3_384_keccak.keccak(BytesToBits(M), 384))
    
    def sha_3_512(self, M):
        return BitsToBytes(self.__sha3_512_keccak.keccak(BytesToBits(M), 512))

class SHAKE_Keccak:
    def __init__(self, c):
        self.__sponge = Sponge(Keccak_f(1600).keccak, self.__pad101, 1600 - c, 1600)
        
    def __pad101(self, x, m):
        assert(x > 0)
        assert(m >= 0)
        
        j = (-m - 2) % x
        return [1] + [0 for _ in range(j)] + [1]
        
    def absorb(self, N):
        self.__sponge.absorb(N + [1, 1, 1, 1])
        
    def squeeze(self, d):
        return self.__sponge.squeeze(d)

class SHAKE:
    
    def __init__(self):
        self.__shake128_keccak = SHAKE_Keccak(256)
        self.__shake256_keccak = SHAKE_Keccak(512)
        
    def shake128(self, M, d):
        self.__shake128_keccak.absorb(BytesToBits(M))
        return BitsToBytes(self.__shake128_keccak.squeeze(d))
        
    def shake256(self, M, d):
        self.__shake256_keccak.absorb(BytesToBits(M))
        return BitsToBytes(self.__shake256_keccak.squeeze(d))
    
class XOF:
    
    def __init__(self):
        self.__shake128_keccak = SHAKE_Keccak(256)
        
    def absorb(self, N):
        self.__shake128_keccak.absorb(BytesToBits(N))
    
    def squeeze(self, l):
        return BitsToBytes(self.__shake128_keccak.squeeze(8 * l))
    
def PRF(eta, s, b):
    assert(eta == 2 or eta == 3)
    assert(len(s)  == 32)
    assert(0 <= b <= 255)
    
    return SHAKE().shake256(s + [b], 8 * 64 * eta)

def H(s):
    return SHA_3().sha_3_256(s)

def J(s):
    return SHAKE().shake256(s, 8 * 32)

def G(c):
    g = SHA_3().sha_3_512(c)
    return g[:32], g[32:]

sha = SHA_3()
print("".join(b2h(BytesToBits(sha.sha_3_256(list("Hello, World!".encode("utf-8")))))))