class Vec3f:
    x: float
    y: float
    z: float
    _counter: int

    def __init__(self, x: float, y: float, z: float):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    # Special Methods
    def __iter__(self):
        self._counter = 0
        return self

    def __next__(self):
        self._counter += 1
        if self._counter > 3:
            raise StopIteration
        return self._get()[self._counter - 1]

    def __eq__(self, other: "Vec3f") -> bool:
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __ne__(self, other):
        return self.x != other.x or self.y != other.y or self.z != other.z

    def __lt__(self, other):
        return self.x < other.x and self.y < other.y and self.z < other.z

    def __le__(self, other):
        return self.x <= other.x and self.y <= other.y and self.z <= other.z

    def __gt__(self, other):
        return self.x > other.x and self.y > other.y and self.z > other.z

    def __ge__(self, other):
        return self.x >= other.x and self.y >= other.y and self.z >= other.z

    def __pos__(self) -> "Vec3f":
        return Vec3f(self.x, self.y, self.z)

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

    # Util Methods
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

    def min(self):
        return min(self._get())

    def max(self):
        return max(self._get())

    def mix(self):
        return min(self._get()), max(self._get())
