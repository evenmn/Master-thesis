class Animal:
    def __init__(self, animal, name, race, age):
        self.animal = animal
        self.name = name
        self.race = race
        self.age = age
		
    def __call__(self):
        return "%s is a %s that's %d years old and of the race %s"%(self.name, self.animal, self.age, self.race)

Alma = Animal("cat", "Alma", "Ragdoll", 4)

print(Alma())
