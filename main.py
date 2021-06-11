import struct
import scipy
import skimage


def main():
    model = PWNCB("./resources/pwncb/bunZipper34834.pwnb")
    calc(model)


def read(file_name):
    with open(file=file_name, mode="r", encoding="UTF-8") as f:
        while line := f.readline():
            yield line


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
    v: list
    vn: list

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

    length: int
    v: list
    vn: list
    support_points: list

    def __init__(self, file_name):
        self.v = []
        self.vn = []
        self.support_points = []

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
        for i in range(self.length // 2):
            self.support_points.append(self.v[i])


if __name__ == "__main__":
    main()
