class MyClass():
    i = 123
    def __init__(self):
        self.i = 345
     

a = MyClass()
print(a.i)
print(MyClass.i)
