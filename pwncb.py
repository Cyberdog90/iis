import struct

from data_class import Vec3f

p = []
n = []

with open("./resources/bunZipper34834.pwnb", "rb") as f:
    data = f.read()

node = int(struct.unpack_from("<i", data, 0)[0])
print(node)

for i in range(node * 3 + 1):
    x, y, z = struct.unpack_from("<ddd", data, 4 + (12 * i))
    p.append(Vec3f(x, y, z))

for i, j in zip(p, range(3)):
    print(i.get())


class PWNCB:
    _file_name: str
    node: int
    p: list
    normal: list

    def __init__(self, file_name):
        self._file_name = file_name
        self.read()

    def read(self):
        i = 0
        with open(self._file_name, "rb") as f:
            while b4b := f.read(4):
                if i == 0:
                    print(b4b.decode())
                    exit(1)

