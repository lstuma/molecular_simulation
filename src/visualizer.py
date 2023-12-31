import pyglet as pg
from random import randrange
from logger import log


elements = None
window = None
callbacks = []

class Window(pg.window.Window):
    def __init__(self):
        super(Window, self).__init__(1080, 720, caption="Molecular Simulator")
        
    def on_draw(self):
        global elements
        self.clear()
        elements.draw()
    
def init():
    global elements
    elements = pg.graphics.Batch()
    # Create window for... well the window tf you think it is..?
    global window
    window = Window()
        
class Circle:
    def __init__(self, pos=(400, 300), radius=10, color=(75,255,125)):
        global elements
        self.shape = pg.shapes.Circle(pos[0], pos[1], radius, color=color, batch=elements)
        
def color(color: str):
    match color:
        case 'magenta': return (randrange(150, 255), randrange(10, 100), randrange(150, 255))
        case 'red': return (randrange(200, 255), randrange(10, 50), randrange(10, 50))
        case 'yellow': return (randrange(200, 255), randrange(150, 200), randrange(10, 50))
        case 'green': return (randrange(10, 100), randrange(150, 255), randrange(10, 50))
        case 'blue': return (randrange(10, 70), randrange(10, 100), randrange(150, 255))
        case 'orange': return (randrange(150, 255), randrange(50, 100), randrange(10, 50))
        case 'rainbow': return (randrange(50, 255), randrange(50, 255), randrange(50, 255))

def run(fixed_callback=None, interval:float=0.2, slow_fixed_callback=None, slow_interval=1):
    global callbacks
    callbacks.append(fixed_callback)
    callbacks.append(slow_fixed_callback)
    if fixed_callback:
        pg.clock.schedule_interval(fixed_callback, interval)
    if slow_fixed_callback:
        pg.clock.schedule_interval(slow_fixed_callback, slow_interval)
    pg.app.run()


if __name__ == '__main__':
    init()
    circles = []
    for i in range(0,50):
        circles.append(Circle(pos=(randrange(5,1075), randrange(5,695))))
    run()
