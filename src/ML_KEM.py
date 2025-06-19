from K_PKE import K_PKE
from keccak import H, G, J
from os import urandom
from conversions import ByteDecode, ByteEncode, b2h, BytesToBits

class ML_KEM:
    
    def __init__(self, k, eta1, eta2, du, dv):
        self.__k = k
        self.__eta1 = eta1
        self.__eta2 = eta2
        self.__du = du
        self.__dv = dv
        self.__k_pke = K_PKE(self.__k, self.__eta1, self.__eta2, self.__du, self.__dv)
        
    def __KeyGen_internal(self, d, z):
        
        (ek_PKE, dk_PKE) = self.__k_pke.KeyGen(d)
        ek = ek_PKE
        dk = dk_PKE + ek + H(ek) + z
        
        return ek, dk
    
    def __Encaps_internal(self, ek, m):
        
        (K, r) = G(m + H(ek))
        c = self.__k_pke.Encrypt(ek, m, r)
        
        return K, c
        
    def __Decaps_internal(self, dk, c):
        
        dk_PKE = dk[:384 * self.__k]
        ek_PKE = dk[384 * self.__k : 768 * self.__k + 32]
        h = dk[768 * self.__k + 32 : 768 * self.__k + 64]
        z = dk[768 * self.__k + 64:]
        
        m_prime = self.__k_pke.Decrypt(dk_PKE, c)
        (K_prime, r_prime) = G(m_prime + h)
        K_barra = J(z + c)
        c_prime = self.__k_pke.Encrypt(ek_PKE, m_prime, r_prime)
        
        if c != c_prime:
            K_prime = K_barra
            
        return K_prime
    
    def KeyGen(self):
        
        d = list(urandom(32))
        z = list(urandom(32))
        
        (ek, dk) = self.__KeyGen_internal(d, z)
        
        return ek, dk
    
    def Encaps(self, ek):
        assert(len(ek) == (384 * self.__k + 32))
        assert(all([0 <= x <= 255 for x in ek]))
        assert([ek[384 * i : 384 * (i + 1)] for i in range(self.__k)] == [ByteEncode(12, ByteDecode(12, ek[384 * i : 384 * (i + 1)])) for i in range(self.__k)])
        
        m = list(urandom(32))
        
        (K, c) = self.__Encaps_internal(ek, m)
        
        return K, c
    
    def Decaps(self, dk, c):
        assert(len(c) == (32 * (self.__du * self.__k + self.__dv)))
        assert(all([0 <= x <= 255 for x in c]))
        assert(len(dk) == (768 * self.__k + 96))
        assert(H(dk[384 * self.__k : 768 * self.__k + 32]) == dk[768 * self.__k + 32 : 768 * self.__k + 64])
        
        K_prime = self.__Decaps_internal(dk, c)
        
        return K_prime
    
class ML_KEM_512:
    
    def __init__(self):
        self.__ml_kem = ML_KEM(2, 3, 2, 10, 4)
    
    def KeyGen(self):
        return self.__ml_kem.KeyGen()
    
    def Encaps(self, ek):
        return self.__ml_kem.Encaps(ek)
    
    def Decaps(self, dk, c):
        return self.__ml_kem.Decaps(dk, c)
    
class ML_KEM_768:
    
    def __init__(self):
        self.__ml_kem = ML_KEM(3, 2, 2, 10, 4)
    
    def KeyGen(self):
        return self.__ml_kem.KeyGen()
    
    def Encaps(self, ek):
        return self.__ml_kem.Encaps(ek)
    
    def Decaps(self, dk, c):
        return self.__ml_kem.Decaps(dk, c)
    
class ML_KEM_1024:
    
    def __init__(self):
        self.__ml_kem = ML_KEM(4, 2, 2, 11, 5)
    
    def KeyGen(self):
        return self.__ml_kem.KeyGen()
    
    def Encaps(self, ek):
        return self.__ml_kem.Encaps(ek)
    
    def Decaps(self, dk, c):
        return self.__ml_kem.Decaps(dk, c)
    

ml = ML_KEM_512()
(ek, dk) = ml.KeyGen()
(K, c) = ml.Encaps(ek)
K_prime = ml.Decaps(dk, c)

if K != K_prime:
    print("¡Fallo en la generación de la clave secreta compartida!\n")
else:
    print(f"Clave de encapsulado(ek) de tamaño {len(ek)} bytes:\n{"".join(b2h(BytesToBits(ek)))}\n")
    print(f"Clave de decapsulado(dk) de tamaño {len(dk)} bytes:\n{"".join(b2h(BytesToBits(dk)))}\n")
    print(f"Texto cifrado (c) de tamaño {len(c)} bytes:\n{"".join(b2h(BytesToBits(c)))}\n")
    print(f"Clave secreta compartida (K) de tamaño {len(K)} bytes:\n{"".join(b2h(BytesToBits(K)))}\n")