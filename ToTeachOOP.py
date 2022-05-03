class MyClass():
    i = 1
    def __init__(self):
        self.i = 2
     

a = MyClass()
print(a.i)
print(MyClass.i)
