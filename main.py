import struct
from os import path
from typing import Union, Optional

import mcubes
import numpy as np

from vector import Vec3f


def main():
    model = File3D("./resources/obj/cube2.obj")


class File3D:
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
    length: int = -1  # 物体表面の節点数
    n_length: int = -1  # 補助節点を含めた全ての節点数
    v: list = []  # 節点
    vn: list = []  # 節点法線
    h: list = []  # 物体表面からの節点の距離

    _precision: int = 1
    file_type: str = ""
    position_min: Vec3f = Vec3f(2 ** 32, 2 ** 32, 2 ** 32)
    position_max: Vec3f = Vec3f(-2 ** 32, -2 ** 32, -2 ** 32)
    _v_min = 2 ** 32
    _v_max = -2 ** 32
    _vn_min = 2 ** 32
    _vn_max = -2 ** 32

    def __init__(self, file_name: str, precision: int = 1) -> None:
        self._precision = precision
        _, self.file_type = path.splitext(file_name)

        if ".pwnb" == self.file_type or ".pwncb" == self.file_type:
            with open(file=file_name, mode="rb") as f:
                data = f.read()
            self.length, *_ = struct.unpack_from("<i", data, 0)
            self.length = int(self.length)
            for i in range(self.length * 2):
                x, y, z, *_ = struct.unpack_from("<ddd", data, 4 + i * 8 * 3)
                if i < self.length:
                    self._v_min = min(x, y, z, self._v_min)
                    self._v_max = max(x, y, z, self._v_max)
                    self._set_v(x=x, y=y, z=z)
                else:
                    self._vn_min = min(x, y, z, self._vn_min)
                    self._vn_max = max(x, y, z, self._vn_max)
                    self._set_vn(x=x, y=y, z=z)
        else:
            with open(file=file_name, mode="r", encoding="UTF-8") as f:
                for line in f.readlines():
                    if (p := line.split())[0] == "v":
                        x, y, z = map(float, p[1:4])
                        self._v_min = min(x, y, z, self._v_min)
                        self._v_max = max(x, y, z, self._v_max)
                        self._set_v(x=x, y=y, z=z)
                    elif p[0] == "vn":
                        x, y, z = map(float, p[1:4])
                        self._vn_min = min(x, y, z, self._vn_min)
                        self._vn_max = max(x, y, z, self._vn_max)
                        self._set_vn(x=x, y=y, z=z)
                    else:
                        continue
            self.length = len(self.v)

        if precision != 1:
            self._precision_f()
        self.h = [0.0] * self.length
        beta = -10 ** -2
        for i in range(0, self.length, 2):
            self.v.append(self.v[i] + self.vn[i] * Vec3f(beta, beta, beta))
            self.h.append(beta)

        self.n_length = len(self.v)
        self._normalize()

    def _set_v(self, x, y, z, index: Optional[int] = None) -> None:
        if index is None:
            self.v.append(Vec3f(x=x, y=y, z=z))
        else:
            self.v[index] = Vec3f(x=x, y=y, z=z)

    def _set_vn(self, x, y, z, index: Optional[int] = None) -> None:
        if index is None:
            self.vn.append(Vec3f(x=x, y=y, z=z))
        else:
            self.vn[index] = Vec3f(x=x, y=y, z=z)

    def _precision_f(self):
        lv = []
        lvn = []
        for i in range(self.length):
            if i % self._precision != 0:
                continue
            lv.append(self.v[i])
            lvn.append(self.vn[i])
        self.length = len(lv)
        self.v = lv
        self.vn = lvn

    def _normalize(self):
        for i in range(self.length):
            x, y, z = map(lambda w: (w - 0.5) * 2,
                          [((j - self._v_min) / (self._v_max - self._v_min)) for j in self.v[i]])
            self.position_min.x = min(x, self.position_min.x)
            self.position_min.y = min(y, self.position_min.y)
            self.position_min.z = min(z, self.position_min.z)
            self._set_v(x=x, y=y, z=z, index=i)

        for i in range(self.length):
            x, y, z = map(lambda w: (w - 0.5) * 2,
                          [((j - self._vn_min) / (self._vn_max - self._vn_min)) for j in self.vn[i]])
            self._set_vn(x=x, y=y, z=z, index=i)

    def mk_obj(self, file_path: str) -> None:
        if self.file_type == ".obj":
            return
        with open(file_path, "w", encoding="UTF-8") as f:
            f.write("# https://github.com/Cyberdog90/iis\n")
            for vertex, _ in zip(self.v, range(self.length)):
                f.write(f"v {vertex.x} {vertex.y} {vertex.z}\n")
            for normal, _ in zip(self.vn, range(self.length)):
                f.write(f"vn {normal.x} {normal.y} {normal.z}\n")


class Itoh:
    model: Union["PWNCB", "Obj"]
    left: list
    right: list
    lambda_vec: list
    alpha_vec: list

    def __init__(self, filename: str):
        # モデルの読み込み
        self.model = File3D(file_name=filename, precision=1)
        print(self.model.n_length)

        self.calc()
        self.lu_decomposition()
        for i in range(self.model.n_length):
            self.func(x=self.model.v[i])
        vertices, triangles = mcubes.marching_cubes_func((-1, -1, -1), (1, 1, 1), 50, 50, 50, self.f_func, 0)
        mcubes.export_obj(vertices, triangles, "./resources/data/bunny.obj")

    def calc(self):
        self.left = []
        for y in range(self.model.n_length + 4):
            tmp = []
            if y == self.model.n_length:
                self.left.append([1.0] * self.model.n_length + [0.0, 0.0, 0.0, 0.0])
            elif y == self.model.n_length + 1:
                self.left.append([i.x for i in self.model.v] + [0.0, 0.0, 0.0, 0.0])
            elif y == self.model.n_length + 2:
                self.left.append([i.y for i in self.model.v] + [0.0, 0.0, 0.0, 0.0])
            elif y == self.model.n_length + 3:
                self.left.append([i.z for i in self.model.v] + [0.0, 0.0, 0.0, 0.0])
            else:
                for x in range(self.model.n_length + 4):
                    if x == self.model.n_length:
                        tmp.append(1.0)
                    elif x == self.model.n_length + 1:
                        tmp.append(self.model.v[y].x)
                    elif x == self.model.n_length + 2:
                        tmp.append(self.model.v[y].y)
                    elif x == self.model.n_length + 3:
                        tmp.append(self.model.v[y].z)
                    else:
                        tmp.append(self.model.v[y].distance(self.model.v[x]) ** 3)
                self.left.append(tmp)
        self.right = self.model.h + [0.0, 0.0, 0.0, 0.0]

    def lu_decomposition(self):
        ans = np.linalg.solve(np.array(self.left), np.array(self.right)).tolist()
        self.lambda_vec, self.alpha_vec = ans[:-4], ans[-4:]

    def f_func(self, x, y, z) -> float:
        return self.func(Vec3f(x, y, z))

    def func(self, x: "Vec3f") -> float:
        def p():
            return self.alpha_vec[0] + self.alpha_vec[1] * x.x + self.alpha_vec[2] * x.y + self.alpha_vec[3] * x.z

        c = 0
        for i in range(self.model.n_length):
            c += self.lambda_vec[i] * x.distance(self.model.v[i]) ** 3
        return c + p()


class MarchingCubes:
    pass


if __name__ == "__main__":
    main()
