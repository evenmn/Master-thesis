import numpy as np

class Shapes:
    def __init__(self, r):
        self.r = r
        
    def getArea(self):
        raise NotImplementedError("Class {} has no instance 'getArea'."
                                  .format(self.__class__.__name__))
        
        
class Circle(Shapes):
    def getArea(self):
        return np.pi*self.r**2
        
    def getCircumference(self):
        return 2*np.pi*self.r
        
    def getExtent(self, x, y):
        if np.linalg.norm([x,y]) < self.r:
            return True
        return False
        
class Square(Shapes):
    def getArea(self):
        return 4*self.r**2
    
    def getCircumference(self):
        return 8*self.r
        
    def getExtent(self, x, y):
        if abs(x) < self.r and abs(y) < self.r:
            return True
        return False
        
class Triangle(Shapes):
    def getCircumference(self):
        return 6*self.r
        

circ = Circle(1)
squr = Square(1)

for i in range(1,10):
    M = 10**i
    A_square = 0
    A_circle = 0
    for m in range(int(M)):
        x = np.random.uniform(-1,1)
        y = np.random.uniform(-1,1)
        if circ.getExtent(x,y):
            A_circle += 1
        if squr.getExtent(x,y):
            A_square += 1
    print("Pi: ", 4 * A_circle/A_square)
