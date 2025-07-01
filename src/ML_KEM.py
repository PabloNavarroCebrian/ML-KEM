from K_PKE import K_PKE
from keccak import H, G, J
from os import urandom
from conversions import ByteDecode, ByteEncode, b2h, BytesToBits
    
class ML_KEM:
    
    def __init__(self, k, eta1, eta2, du, dv):
        """
        Inicializa una instancia del esquema ML-KEM con los parámetros dados.

        Entrada:
        - k: número de polinomios en el esquema (determina el nivel de seguridad)
        - eta1: parámetro de ruido para la generación de claves
        - eta2: parámetro de ruido para el cifrado
        - du, dv: parámetros de compresión de la cápsula

        Internamente, se instancia una versión correspondiente del esquema K-PKE.
        """
        self.__k = k
        self.__eta1 = eta1
        self.__eta2 = eta2
        self.__du = du
        self.__dv = dv
        self.__k_pke = K_PKE(self.__k, self.__eta1, self.__eta2, self.__du, self.__dv)
        
    def __KeyGen_internal(self, d, z):
        """
        Algoritmo interno de generación de claves.

        Entrada:
        - d: semilla aleatoria para generación determinista
        - z: cadena aleatoria utilizada para KDF alternativa en caso de fallo de descifrado

        Salida:
        - ek: clave pública
        - dk: clave privada extendida (incluye ek, H(ek), z)
        """
        (ek_PKE, dk_PKE) = self.__k_pke.KeyGen(d)
        ek = ek_PKE
        dk = dk_PKE + ek + H(ek) + z
        
        return ek, dk
    
    def __Encaps_internal(self, ek, m):
        """
        Algoritmo interno de encapsulación.

        Entrada:
        - ek: clave pública del receptor
        - m: mensaje aleatorio (preimagen de la clave)

        Salida:
        - K: clave simétrica derivada mediante función hash
        - c: cápsula (ciphertext) que encapsula el mensaje m
        """
        (K, r) = G(m + H(ek))
        c = self.__k_pke.Encrypt(ek, m, r)
        
        return K, c
        
    def __Decaps_internal(self, dk, c):
        """
        Algoritmo interno de desencapsulación.

        Entrada:
        - dk: clave privada extendida del receptor
        - c: cápsula recibida

        Salida:
        - K': clave simétrica recuperada (o clave alternativa si el descifrado falla)
        """
        # Se extraen las partes necesarias de la clave privada
        dk_PKE = dk[:384 * self.__k]
        ek_PKE = dk[384 * self.__k : 768 * self.__k + 32]
        h = dk[768 * self.__k + 32 : 768 * self.__k + 64]
        z = dk[768 * self.__k + 64:]
        
        # Se intenta recuperar el mensaje original
        m_prime = self.__k_pke.Decrypt(dk_PKE, c)
        (K_prime, r_prime) = G(m_prime + h)
        K_barra = J(z + c)  # Clave alternativa en caso de fallo
        c_prime = self.__k_pke.Encrypt(ek_PKE, m_prime, r_prime)
        
        # Se comprueba si el descifrado fue correcto
        if c != c_prime:
            K_prime = K_barra
            
        return K_prime
    
    def KeyGen(self):
        """
        Genera un par de claves pública y privada para ML-KEM.

        Salida:
        - ek: clave pública
        - dk: clave privada extendida
        """
        d = list(urandom(32))
        z = list(urandom(32))
        
        (ek, dk) = self.__KeyGen_internal(d, z)
        
        return ek, dk
    
    def Encaps(self, ek):
        """
        Realiza el algoritmo de encapsulación usando una clave pública.

        Entrada:
        - ek: clave pública del receptor

        Salida:
        - K: clave simétrica generada
        - c: cápsula correspondiente
        """
        assert(len(ek) == (384 * self.__k + 32))
        assert(all([0 <= x <= 255 for x in ek]))
        # Verifica que la clave pública es válida según el estándar
        assert([ek[384 * i : 384 * (i + 1)] for i in range(self.__k)] == [ByteEncode(12, ByteDecode(12, ek[384 * i : 384 * (i + 1)])) for i in range(self.__k)])
        
        m = list(urandom(32))  # Mensaje aleatorio que se encapsula
        
        (K, c) = self.__Encaps_internal(ek, m)
        
        return K, c
    
    def Decaps(self, dk, c):
        """
        Realiza el algoritmo de desencapsulación usando una clave privada.

        Entrada:
        - dk: clave privada extendida del receptor
        - c: cápsula recibida

        Salida:
        - K': clave simétrica recuperada
        """
        # Verificaciones de integridad sobre cápsula y clave
        assert(len(c) == (32 * (self.__du * self.__k + self.__dv)))
        assert(all([0 <= x <= 255 for x in c]))
        assert(len(dk) == (768 * self.__k + 96))
        assert(H(dk[384 * self.__k : 768 * self.__k + 32]) == dk[768 * self.__k + 32 : 768 * self.__k + 64])
        
        K_prime = self.__Decaps_internal(dk, c)
        
        return K_prime


class ML_KEM_512:
    
    def __init__(self):
        """
        Inicializa una instancia ML-KEM con parámetros correspondientes al nivel de seguridad 1 (512).
        """
        self.__ml_kem = ML_KEM(2, 3, 2, 10, 4)
    
    def KeyGen(self):
        """
        Ejecuta la generación de claves para ML-KEM-512.
        """
        return self.__ml_kem.KeyGen()
    
    def Encaps(self, ek):
        """
        Ejecuta la encapsulación con clave pública para ML-KEM-512.
        """
        return self.__ml_kem.Encaps(ek)
    
    def Decaps(self, dk, c):
        """
        Ejecuta la desencapsulación con clave privada para ML-KEM-512.
        """
        return self.__ml_kem.Decaps(dk, c)


class ML_KEM_768:
    
    def __init__(self):
        """
        Inicializa una instancia ML-KEM con parámetros correspondientes al nivel de seguridad 3 (768).
        """
        self.__ml_kem = ML_KEM(3, 2, 2, 10, 4)
    
    def KeyGen(self):
        """
        Ejecuta la generación de claves para ML-KEM-768.
        """
        return self.__ml_kem.KeyGen()
    
    def Encaps(self, ek):
        """
        Ejecuta la encapsulación con clave pública para ML-KEM-768.
        """
        return self.__ml_kem.Encaps(ek)
    
    def Decaps(self, dk, c):
        """
        Ejecuta la desencapsulación con clave privada para ML-KEM-768.
        """
        return self.__ml_kem.Decaps(dk, c)


class ML_KEM_1024:
    
    def __init__(self):
        """
        Inicializa una instancia ML-KEM con parámetros correspondientes al nivel de seguridad 5 (1024).
        """
        self.__ml_kem = ML_KEM(4, 2, 2, 11, 5)
    
    def KeyGen(self):
        """
        Ejecuta la generación de claves para ML-KEM-1024.
        """
        return self.__ml_kem.KeyGen()
    
    def Encaps(self, ek):
        """
        Ejecuta la encapsulación con clave pública para ML-KEM-1024.
        """
        return self.__ml_kem.Encaps(ek)
    
    def Decaps(self, dk, c):
        """
        Ejecuta la desencapsulación con clave privada para ML-KEM-1024.
        """
        return self.__ml_kem.Decaps(dk, c)

    
import time

ml = ML_KEM_768()

start = time.time()
(ek, dk) = ml.KeyGen()
end = time.time()
t1 = end - start
start = time.time()
(K, c) = ml.Encaps(ek)
end = time.time()
t2 = end - start
start = time.time()
K_prime = ml.Decaps(dk, c)
end = time.time()
t3 = end - start

if K != K_prime:
    print("¡Fallo en la generación de la clave secreta compartida!\n")
else:
    print(f"Clave de encapsulado(ek) de tamaño {len(ek)} bytes:\n{"".join(b2h(BytesToBits(ek)))}\n")
    print(f"Clave de desencapsulado(dk) de tamaño {len(dk)} bytes:\n{"".join(b2h(BytesToBits(dk)))}\n")
    print(f"Texto cifrado (c) de tamaño {len(c)} bytes:\n{"".join(b2h(BytesToBits(c)))}\n")
    print(f"Clave secreta compartida (K) de tamaño {len(K)} bytes:\n{"".join(b2h(BytesToBits(K)))}\n")
        
    print(f"Tiempo de Generación de Claves: {t1:.3f} segundos")
    print(f"Tiempo de Encapsulado: {t2:.3f} segundos")
    print(f"Tiempo de Desencapsulado: {t3:.3f} segundos")
    print(f"Tiempo total: {t1 + t2 + t3:.3f} segundos")