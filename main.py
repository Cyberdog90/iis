import struct
from typing import Union
import numpy as np


def main():
    model = Obj("./resources/obj/tri-pyramid.obj")
    print(f"接点数{model.length}")
    print(f"補助接点を含めた節点数:{len(model.v)}")
    print(f"法線数:{len(model.vn)}")

    for i in range(model.n_length):
        if i >= model.length:
            print(f"補助接点:{i + 1 - model.length} -> x:{model.v[i].x}, y:{model.v[i].y}, z:{model.v[i].z}")
        else:
            print(f"接点:{i + 1} -> x:{model.v[i].x}, y:{model.v[i].y}, z:{model.v[i].z}")

    for i in range(len(model.h)):
        if i >= model.length:
            print(f"補助接点{i + 1 - model.length}の物体表面からの距離 -> {model.h[i]}")
        else:
            print(f"接点{i + 1}の物体表面からの距離 -> {model.h[i]}")

    left, right = calc(model=model)

    with open("./resources/data/left.csv", "w", encoding="UTF-8") as f:
        for i in left:
            f.write(", ".join(list(map(str, i))) + "\n")
    with open("./resources/data/right.csv", "w", encoding="UTF-8") as f:
        for i in right:
            f.write(str(i) + "\n")

    ans = np.linalg.solve(np.array(left), np.array(right)).tolist()

    with open("./resources/data/lu.csv", "w", encoding="UTF-8") as f:
        for i in ans:
            f.write(str(i) + "\n")

    lambda_v, alpha_v = ans[:-4], ans[-4:]

    for i in range(model.n_length):
        fun = func(model.v[i], model, lambda_v, alpha_v)
        print(f"f(x_{i}) -> {fun}")


def read(file_name):
    with open(file=file_name, mode="r", encoding="UTF-8") as f:
        while line := f.readline():
            yield line


def mk_obj(data: "PWNCB", file_path: str) -> None:
    with open(file_path, "w", encoding="UTF-8") as f:
        f.write("# https://github.com/Cyberdog90/iis\n")
        for vertex, _ in zip(data.v, range(data.length)):
            f.write(f"v {vertex.x} {vertex.y} {vertex.z}\n")
        for normal, _ in zip(data.vn, range(data.length)):
            f.write(f"vn {normal.x} {normal.y} {normal.z}\n")


def calc(model: Union["PWNCB", "Obj"]) -> tuple:
    left = []
    for y in range(model.n_length + 4):
        tmp = []
        if y == model.n_length:
            left.append([1.0] * model.n_length + [0.0, 0.0, 0.0, 0.0])
        elif y == model.n_length + 1:
            left.append([i.x for i in model.v] + [0.0, 0.0, 0.0, 0.0])
        elif y == model.n_length + 2:
            left.append([i.y for i in model.v] + [0.0, 0.0, 0.0, 0.0])
        elif y == model.n_length + 3:
            left.append([i.z for i in model.v] + [0.0, 0.0, 0.0, 0.0])
        else:
            for x in range(model.n_length + 4):
                if x == model.n_length:
                    tmp.append(1.0)
                elif x == model.n_length + 1:
                    tmp.append(model.v[y].x)
                elif x == model.n_length + 2:
                    tmp.append(model.v[y].y)
                elif x == model.n_length + 3:
                    tmp.append(model.v[y].z)
                else:
                    tmp.append(model.v[y].distance(model.v[x]) ** 3)
            left.append(tmp)
    right = model.h + [0.0, 0.0, 0.0, 0.0]
    return left, right


def func(x: "Vec3f", model, lambda_v, a):
    def p():
        return a[0] + a[1] * x.x + a[2] * x.y + a[3] * x.z

    c = 0
    for i in range(model.n_length):
        c += lambda_v[i] * x.distance(model.v[i]) ** 3
    return c + p()


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

    def __mul__(self, other: "Vec3f") -> "Vec3f":
        return Vec3f(self.x * other.x, self.y * other.y, self.z * other.z)

    def __truediv__(self, other: "Vec3f") -> "Vec3f":
        return Vec3f(self.x / other.x, self.y / other.y, self.z / other.z)

    def __floordiv__(self, other: "Vec3f") -> "Vec3f":
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

    def distance(self, arg: "Vec3f") -> float:
        return ((self.x - arg.x) ** 2 + (self.y - arg.y) ** 2 + (self.z - arg.z) ** 2) ** .5

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
        self.h = [0.0] * self.length
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
        self.h = [0.0] * self.length
        beta = -10 ** -2
        for i in range(0, self.length, 2):
            self.v.append(self.v[i] + self.vn[i] * Vec3f(beta, beta, beta))
            self.h.append(beta)

        self.n_length = len(self.v)


if __name__ == "__main__":
    main()
