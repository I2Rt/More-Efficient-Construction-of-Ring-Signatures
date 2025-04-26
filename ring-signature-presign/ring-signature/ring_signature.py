import math
import random

import numpy as np
import secrets
from functools import lru_cache

from rings import Polynomial
from hash import expand_seed, generate_matrix_from_seed, generate_challenge


# Q = 2^29 # q ≈ 2^29
N = 64  # n = 64
Q = 2^26

# Extremely relaxed bounds (for testing)
GAMMA1 = Q
GAMMA2 = Q
BETA = 1

# Domain separators
DOMAIN_SMALL_POLY = 0x03
DOMAIN_Y_POLY = 0x05


class rtScheme:
    def __init__(self):
        self.A0 = None
        self.A = None
        self.A_bar = None
        self.B = None

        self.rho = None

        self.rho_A = None
        self.rho_B = None

        self.b = None
        self.s1 = None
        self.s2 = None
        self.s = None
        self.y = None

        # self.m = 7
        # self.k = 23
        self.m = 6
        self.k = 22
        self.eta = 1
        self.m2 = 16
        self.s_frak = 1.3 * 140 * 1 * math.sqrt(self.k * N)
        self.beta = self.s_frak * math.sqrt(2 * self.m * N)

        # self.n = 17
        self.n = 14
        self.kappa = 10
        self.nu = 0 # 0, 1, 2, 3, 5, 8
        self.m1 = self.nu + 2
        # self.m1 = self.k + self.nu + 3

        # NIZK Params
        self.A1 = None
        self.A2 = None
        self.By = None
        self.Bx = None
        self.Bg = None
        self.B_hat = None
        self.s1 = None
        self.r = None
        self.x = None
        self.D2 = None
        self.d1 = None
        self.u = None

    def setup(self):
        # pre-sign setup
        self.rho_A = secrets.token_bytes(32)
        self.rho_B = secrets.token_bytes(32)

        # 生成矩阵 A (维度 k × m)
        self.A = self.get_matrix_A()
        self.A0 = self.A[:self.k - self.m - 1]  # 提取 A 的前 k - m - 1 行（即 A0）
        self.B = self.get_matrix_B()

        self.NIZK_Setup()

    def rs_keygen(self):
        """Generate a new keypair"""

        # 1. 生成小向量 s1 和 s2
        self.s1 = self.generate_small_vector(self.k - self.m - 1)
        self.s2 = self.generate_small_vector(self.m)

        # 2. 构造私钥 s = (1, s1, s2)
        self.s = [Polynomial([1])] + self.s1 + self.s2  # 1 是常数多项式 1

        # 3. 计算 b = A0 * s1 + s2
        self.b = self._matrix_multiply(self.A0, self.s1)
        for i in range(self.m):
            self.b[i] = self.b[i] + self.s2[i]

        # 4. 返回公钥 b 和私钥 s
        return (self.rho_A, self.b), self.s
        # return self.b, self.s


    def rs_sign(self, message: bytes):
        if not all([ self.b, self.s1, self.s2]):
            raise ValueError("Keys not generated")

        # Start of sampling block
        while True:
            A = self.build_matrix_A(self.b, self.A0)
            self.y = self._sample_gaussian_polynomial_vector(self.k, self.s_frak)
            w = self._matrix_multiply(A, self.y)
            self.r = self.generate_small_vector(self.m2)

            com = self._matrix_multiply(self.B, self.r)
            for i in range(self.m):
                com[i] = com[i] + w[i]
            # Generate challenge polynomial c
            c = generate_challenge(message, (self.rho_A, self.b), com)

            # Compute z = y + cs
            z = []
            for i in range(self.k):
                z_poly = self.y[i] + (c * self.s[i])
                z.append(z_poly)

            # 检查z的2范数( | | z | |₂ < B)
            z_coeffs = np.concatenate([np.array(p.coefficients) for p in z])
            z_squared = np.sum(z_coeffs ** 2)
            z_norm = np.sqrt(z_squared)
            if z_norm >= self.beta:
                continue # 回到 self.y = self._sample_gaussian_polynomial_vector(self.k, self.s_frak) 重新执行
            else:
                break  # Exit the loop when we have a valid signature

        return z, c, com


    def rs_verify(self, message: bytes, signature: tuple):
        z, c, com = signature

        # # 1. 检查挑战值 (||z||∞ < Q/4)
        c_1 = generate_challenge(message, (self.rho_A, self.b), com)
        if c_1 == c:
            print(f"Reject: c is invalid")
            return False

        z_coeffs = np.concatenate([np.array(p.coefficients) for p in z])
        # if np.any(np.abs(z_coeffs) >= Q // 4):
        #     print("Reject: z coefficients too large (infinity norm)")
        #     return False

        # 2. 检查z的2范数 (||z||₂ < B)
        z_squared = np.sum(z_coeffs ** 2)
        z_norm = np.sqrt(z_squared)
        if z_norm >= self.beta:
            print(f"Reject: z norm too large (2-norm: {z_norm} >= {self.B})")
            return False

        return True

    def NIZK_Setup(self):
        # self.r = [Polynomial([1])] + self.r
        self.u = self._sample_gaussian_polynomial_vector(self.nu + 1, self.s_frak)
        self.x = self._sample_gaussian_polynomial_vector(self.nu, self.s_frak)

        # rho_A0 = secrets.token_bytes(32)
        # A0_1 = generate_matrix_from_seed(rho_A0, self.n, self.m1 - self.n)
        # A0_2 = generate_matrix_from_seed(rho_A0, self.n, self.m2 - self.n - 1)
        # b1 = self._matrix_multiply(A0_1, self.s1[:self.m1-self.n])
        # for i in range(self.n):
        #     b1[i] = b1[i] + self.s1[self.m1-self.n+i]
        #
        # b2 = self._matrix_multiply(A0_2, self.r[:self.m1-self.n-1])
        # for i in range(self.n):
        #     b2[i] = b2[i] + self.r[self.m1-self.n-1+i]

        rho_A1 = secrets.token_bytes(32)
        rho_A2 = secrets.token_bytes(32)
        rho_B = secrets.token_bytes(32)

        self.A1 = generate_matrix_from_seed(rho_A1, self.n, self.m1 + 1)
        self.A2 = generate_matrix_from_seed(rho_A2, self.n, self.m2)

        self.By = generate_matrix_from_seed(rho_B, self.k, self.m2)
        self.Bx = generate_matrix_from_seed(rho_B, self.nu, self.m2)
        self.Bg = generate_matrix_from_seed(rho_B, self.kappa, self.m2)

        self.B_hat = [
            self.By,
            self.Bx,
            self.Bg
        ]

        rho_D2 = secrets.token_bytes(32)
        self.D2 = generate_matrix_from_seed(rho_D2, (self.m1+1+self.m+self.k+self.nu+self.kappa)*2, (self.m1+1+self.m+self.k+self.nu+self.kappa)*2)
        self.d1 = self._sample_gaussian_polynomial_vector((self.m1+1+self.m+self.k+self.nu+self.kappa)*2, self.s_frak)

        # return (A1, A2, By, Bx, Bg, s1, r)

    def NIZK_Prove(self, message):
        self.r = [Polynomial([1])] + self.r
        s_star = self.s + [Polynomial([1])]
        self.s1 = [Polynomial([1])] + [Polynomial([1])] + s_star + self.u

        y1 = self._sample_gaussian_polynomial_vector(self.m1+1, self.s_frak)
        y2 = self._sample_gaussian_polynomial_vector(self.m2, self.s_frak)
        g = self._sample_gaussian_polynomial_vector(self.kappa, self.s_frak)
        # x = self._sample_gaussian_polynomial_vector(self.nu, self.s_frak)

        # 计算 w = A1 * y1 + A2 * y2
        w = self._matrix_multiply(self.A1, y1)
        w_tmp = self._matrix_multiply(self.A2, y2)
        for i in range(self.n):
            w[i] = w[i] + w_tmp[i]

        ty = self._matrix_multiply(self.By, self.r)
        for i in range(self.k):
            ty[i] = ty[i] + self.y[i]

        tx = self._matrix_multiply(self.Bx, self.r)
        for i in range(self.nu):
            tx[i] = tx[i] + self.x[i]

        tg = self._matrix_multiply(self.Bg, self.r)
        for i in range(self.kappa):
            tg[i] = tg[i] + g[i]

        # 模拟生成挑战 gama 的时间
        # gama = generate_matrix_from_seed(self.rho_A, self.kappa, 64 + 10 + 4*self.nu + self.k)
        gama_ij = generate_challenge(message, (self.rho_A, self.b), w)
        for i in range(self.kappa * (64 + 10 + 4*self.nu + self.k)-1):
            generate_challenge(message, (self.rho_A, self.b), w)

        # 模拟计算 h 的时间
        # g_i = Polynomial([random.random()])
        # f_1 = Polynomial([random.random()])
        # f_2 = Polynomial([random.random()])
        # f_3 = Polynomial([random.random()])
        # f_4 = Polynomial([random.random()])
        # f_5 = Polynomial([random.random()])
        # f_6 = Polynomial([random.random()])
        # f_7 = Polynomial([random.random()])
        # for i in range(self.kappa):
        #     for j in range(64):
        #         gama_ij * f_1
        #     for j in range(10):
        #         gama_ij * f_2

        # h = self._sample_gaussian_polynomial_vector(self.kappa, self.s_frak)

        # 模拟生成挑战 mu 的时间
        # mu = self._sample_gaussian_polynomial_vector(self.kappa, self.s_frak)
        for i in range(self.kappa):
            generate_challenge(message, (self.rho_A, self.b), w)

        m = w + self.y + self.x + g

        s_hat = self.s1 + self.s1 + m + m

        # y_hat = y1 + y1 + self._matrix_multiply(self.B_hat, y2) + self._matrix_multiply(self.B_hat, y2)

        # f1 = self.compute_f1(s_hat, y_hat, self.D2, self.D2, y_hat)

        # t = self.vector_dot_product(self.b, self.r) + f1

        # 生成挑战 c
        c = generate_challenge(message, (self.rho_A, self.b), w)
        for i in range(self.m2):
            self.r[i] = c * self.r[i]

        for i in range(self.m1+1):
            self.s1[i] = c * self.s1[i]

        z1 = y1 + self.s1

        z2 = y2 + self.r

        return (ty, tx, tg, z1, z2)





    def _sample_gaussian_polynomial_vector(self, k: int, s: float) -> list:
        """
        从标准差为 s 的高斯分布上采样一个长度为 k 的向量，
        其中每个元素是一个 Polynomial 对象，其系数从高斯分布采样。
        """
        vector = []
        for _ in range(k):
            # 采样 N 个高斯分布系数
            coeffs = np.random.normal(0, s, size=Polynomial.N)
            # 四舍五入并取模 Q
            coeffs = np.round(coeffs).astype(np.int32) % Polynomial.Q
            # 创建 Polynomial 对象
            poly = Polynomial(coeffs)
            vector.append(poly)
        return vector

    def vector_dot_product(v1: list, v2: list) -> Polynomial:
        """
        计算两个多项式向量的点积（内积）

        参数:
            v1: 第一个向量 (列表 of Polynomial)
            v2: 第二个向量 (列表 of Polynomial)

        返回:
            点积结果 (单个 Polynomial)
        """
        if len(v1) != len(v2):
            raise ValueError("向量长度必须相同")

        result = Polynomial([0])  # 初始化零多项式
        for a, b in zip(v1, v2):
            result += a * b  # 多项式乘法和加法

        return result


    def build_matrix_A(self, b, A0):
        """Construct the matrix A = [-2b + q*j | 2*A0 | 2*Im]"""
        m = self.m
        q = Q

        # 1. 计算 -2*b + q*j (j 是全 1 向量)
        j = [Polynomial([1]) for _ in range(m)]  # j = [1, 1, ..., 1] (长度为 m)
        neg_2b_plus_qj = [
            (Polynomial([-2]) * b_i + Polynomial([q]) * j_i)  # -2*b_i + q*1
            for b_i, j_i in zip(b, j)
        ]

        # 2. 计算 2*A0 (A0 是 (k - m - 1) × m 矩阵)
        A0_scaled = [
            [Polynomial([2]) * poly for poly in row]  # 2*A0[i][j]
            for row in A0
        ]

        # 3. 计算 2*Im (单位矩阵)
        Im = []
        for i in range(m):
            row = []
            for j in range(m):
                # 对角线元素为 2，其他为 0
                element = Polynomial([2]) if i == j else Polynomial([0])
                row.append(element)
            Im.append(row)

        # 4. 水平拼接 [-2b + q*j | 2*A0 | 2*Im]
        # 先转置 A0 和 Im 以便按列拼接
        A0_T = list(zip(*A0_scaled))  # A0^T (m × (k - m - 1))
        Im_T = list(zip(*Im))  # I_m^T (m × m)

        # 拼接每行: [-2b_i + q | 2*A0_i | 2*Im_i]
        A = [
            [neg_2b_plus_qj[i]] + list(A0_T[i]) + list(Im_T[i])  # 水平拼接
            for i in range(m)
        ]

        return A


    def _generate_vector(self, size: int, bound: int, domain: int) -> list:
        seed = secrets.token_bytes(32)
        total_coeffs = size * Polynomial.N
        randomness = expand_seed(seed, domain, total_coeffs * 4)
        coeffs = np.frombuffer(randomness, dtype=np.uint32)
        coeffs = coeffs % (2 * bound + 1)
        coeffs = coeffs.astype(np.int32) - bound
        coeffs = coeffs.reshape(size, Polynomial.N)
        return [Polynomial(row.tolist()) for row in coeffs]

    def generate_y_vector(self) -> list:
        return self._generate_vector(self.l, GAMMA1, DOMAIN_Y_POLY)

    def generate_small_vector(self, size: int) -> list:
        return self._generate_vector(size, self.eta, DOMAIN_SMALL_POLY)

    @lru_cache(maxsize=None)  # Changed to support multiple cache entries
    def get_matrix_A(self):
        """Get cached matrix A using rho as part of cache key"""
        if self.rho_A is None:
            raise ValueError("Keys not generated")

        # Include rho in computation to make it part of cache key
        matrix = generate_matrix_from_seed(self.rho_A, self.m, self.k)
        return matrix

    def get_matrix_B(self):
        """Get cached matrix A using rho as part of cache key"""
        if self.rho_B is None:
            raise ValueError("Keys not generated")

        # Include rho in computation to make it part of cache key
        matrix = generate_matrix_from_seed(self.rho_B, self.m, self.m2)
        return matrix

    def _matrix_multiply(self, matrix: list, vector: list) -> list:
        result = []
        for row in matrix:
            sum_poly = Polynomial()
            for a, v in zip(row, vector):
                sum_poly = sum_poly + (a * v)
            result.append(sum_poly)
        return result

    def vector_transpose_multiply(d: list, y: list) -> Polynomial:
        """
        计算 dᵀy = d的转置乘y (向量内积)

        参数:
            d: 第一个向量 (列表 of Polynomial)
            y: 第二个向量 (列表 of Polynomial)

        返回:
            内积结果 (单个 Polynomial)

        注意:
            1. 两个向量长度必须相同
            2. 结果是多项式环中的元素
        """
        if len(d) != len(y):
            raise ValueError("向量长度不匹配")

        result = Polynomial([0])  # 初始化零多项式
        for di, yi in zip(d, y):
            result += di * yi  # 多项式乘法然后累加

        return result

    def compute_f1(self, s: list, y: list, D1: list, D2: list, d: list) -> Polynomial:
        """
        计算 f1 = sD1y + yD2s + dᵀy

        参数:
            s: 私钥向量
            y: 随机向量
            D1: 第一矩阵
            D2: 第二矩阵
            d: 向量d

        返回:
            计算结果多项式
        """
        # 1. 计算 sD1y
        D1y = self._matrix_multiply(D1, y)
        sD1y = self.vector_transpose_multiply(s, D1y)

        # 2. 计算 yD2s
        D2s = self._matrix_multiply(D2, s)
        yD2s = self.vector_transpose_multiply(y, D2s)

        # 3. 计算 dᵀy
        dTy = self.vector_transpose_multiply(d, y)

        # 4. 合并结果 f1 = sD1y + yD2s + dᵀy
        f1 = sD1y + yD2s + dTy

        return f1