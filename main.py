import struct
from scipy.sparse.linalg import SuperLU
import skimage

import time


def main():
    a = 0.0
    for i in range(100):
        a += test()
        print()
    print(a / 10)


def test():
    N = 10
    print(f"試行回数{N}回")
    start = time.time()
    for i in range(N):
        model = PWNCB("./resources/pwncb/bunZipper34834.pwnb")
        model.length += 1
    print("バイナリ")
    print(f"経過時間:{time.time() - start}秒")
    binary = time.time() - start

    # print(len(model.v), len(model.vn), len(model.h), model.length, model.n_length)

    start = time.time()
    for i in range(N):
        model2 = Obj("./resources/obj/bunZipper34834.obj")
        model2.length += 1
    print("obj")
    print(f"経過時間:{time.time() - start}秒")
    obj = time.time() - start

    # print(len(model2.v), len(model2.vn), len(model2.h), model2.length, model2.n_length)

    return obj / binary


def read(file_name):
    with open(file=file_name, mode="r", encoding="UTF-8") as f:
        while line := f.readline():
            yield line


def mk_obj(data: "PWNCB", file_path):
    with open(file_path, "w", encoding="UTF-8") as f:
        f.write("# https://github.com/Cyberdog90/iis\n")
        for vertex, _ in zip(data.v, range(data.length)):
            f.write(f"v {vertex.x} {vertex.y} {vertex.z}\n")
        for normal, _ in zip(data.vn, range(data.length)):
            f.write(f"vn {normal.x} {normal.y} {normal.z}\n")


def calc(model):
    for i in range(model.length):
        pass


class Vec3f:
    x: float
    y: float
    z: float
    _counter: int

    def __init__(self, x: float, y: float, z: float):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    # |-----Iterator Methods-----|
    def __iter__(self):
        self._counter = 0
        return self

    def __next__(self):
        self._counter += 1
        if self._counter > 3:
            raise StopIteration
        return self._get()[self._counter - 1]

    # |-----Comparison Methods-----|

    # "=="operator
    def __eq__(self, other: "Vec3f") -> bool:
        return self.x == other.x and self.y == other.y and self.z == other.z

    # "!="operator
    def __ne__(self, other):
        return self.x != other.x or self.y != other.y or self.z != other.z

    # "<"operator
    def __lt__(self, other):
        return self.x < other.x and self.y < other.y and self.z < other.z

    # "<="operator
    def __le__(self, other):
        return self.x <= other.x and self.y <= other.y and self.z <= other.z

    # ">"operator
    def __gt__(self, other):
        return self.x > other.x and self.y > other.y and self.z > other.z

    # ">="operator
    def __ge__(self, other):
        return self.x >= other.x and self.y >= other.y and self.z >= other.z

    # |-----Arithmetic Methods-----|
    # "+(single)"operator
    def __pos__(self) -> "Vec3f":
        return Vec3f(self.x, self.y, self.z)

    # "-(single)"operator
    def __neg__(self) -> "Vec3f":
        return Vec3f(-self.x, -self.y, -self.z)

    def __add__(self, other: "Vec3f") -> "Vec3f":
        return Vec3f(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: "Vec3f") -> "Vec3f":
        return Vec3f(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):
        return Vec3f(self.x * other.x, self.y * other.y, self.z * other.z)

    def __truediv__(self, other):
        return Vec3f(self.x / other.x, self.y / other.y, self.z / other.z)

    def __floordiv__(self, other):
        return Vec3f(self.x // other.x, self.y // other.y, self.z // other.z)

    def __pow__(self, power, modulo=None):
        if modulo is None:
            return Vec3f(self.x ** power, self.y ** power, self.z ** power)
        else:
            return Vec3f(self.x ** power % modulo, self.y ** power % modulo, self.z ** power % modulo)

    def __abs__(self):
        return Vec3f(abs(self.x), abs(self.y), abs(self.z))

    def __len__(self):
        return 3

    def increment(self) -> None:
        self.x += 1
        self.y += 1
        self.z += 1

    def decrement(self) -> None:
        self.x -= 1
        self.y -= 1
        self.z -= 1

    def _get(self) -> list:
        return [self.x, self.y, self.z]


class Obj:
    length: int  # 物体表面の節点数
    n_length: int  # 補助節点を含めた全ての節点数
    v: list
    vn: list
    h: list

    def __init__(self, file_name):
        self.v = []
        self.vn = []
        for line in read(file_name=file_name):
            if (p := line.split())[0] == "#":
                continue

            elif p[0] == "v":
                x, y, z = map(float, p[1:4])
                self.v.append(Vec3f(x, y, z))
            elif p[0] == "vn":
                x, y, z = map(float, p[1:4])
                self.vn.append(Vec3f(x, y, z))

            else:
                continue

        self.length = len(self.v)
        self.h = [0] * self.length
        beta = -10 ** -2
        for i in range(0, self.length, 2):
            self.v.append(self.v[i] + self.vn[i] * Vec3f(beta, beta, beta))
            self.h.append(beta)

        self.n_length = len(self.v)


class PWNCB:
    # PWNCB File Specification
    # Little Endian

    # |----File Length----|
    #      (4Bytes:Int)
    # |-----Vertex X0-----||-----Vertex Y0-----||-----Vertex Z0-----|
    #    (8Bytes:Double)      (8Bytes:Double)      (8Bytes:Double)
    #           ︙                   ︙                   ︙
    # |-----Vertex Xn-----||-----Vertex Yn-----||-----Vertex Zn-----|
    #    (8Bytes:Double)      (8Bytes:Double)      (8Bytes:Double)
    # |-----Normal X0-----||-----Normal Y0-----||-----Normal Z0-----|
    #    (8Bytes:Double)      (8Bytes:Double)      (8Bytes:Double)
    #           ︙                   ︙                   ︙
    # |-----Normal Xn-----||-----Normal Yn-----||-----Normal Zn-----|
    #    (8Bytes:Double)      (8Bytes:Double)      (8Bytes:Double)

    length: int  # 物体表面の節点数
    n_length: int  # 補助節点を含めた全ての節点数
    v: list
    vn: list
    h: list

    def __init__(self, file_name):
        self.v = []
        self.vn = []

        with open(file=file_name, mode="rb") as f:
            data = f.read()

        # バイナリの長さを読み取り
        self.length, *_ = struct.unpack_from("<i", data, 0)
        self.length = int(self.length)

        # 頂点座標の読み取り
        for i in range(self.length):
            x, y, z, *_ = struct.unpack_from("<ddd", data, 4 + i * 8 * 3)
            self.v.append(Vec3f(x, y, z))

        # 頂点法線の読み取り
        for i in range(self.length):
            x, y, z, *_ = struct.unpack_from("<ddd", data, 4 + (self.length * 8 * 3) + i * 8 * 3)
            self.vn.append(Vec3f(x, y, z))

        # 補助接点の生成
        self.h = [0] * self.length
        beta = -10 ** -2
        for i in range(0, self.length, 2):
            self.v.append(self.v[i] + self.vn[i] * Vec3f(beta, beta, beta))
            self.h.append(beta)

        self.n_length = len(self.v)


if __name__ == "__main__":
    main()
