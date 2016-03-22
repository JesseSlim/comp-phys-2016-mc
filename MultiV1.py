from multiprocessing import Process, Queue
import random
import time

def f(q, name):
    print('a')
    for i in range(1,1000):
        time.sleep(0.001)
    
starttime = time.clock()

if __name__ == '__main__':
    q = Queue()
    a = Process(target=f, args=(q,'a', ))
    b = Process(target=f, args=(q,'b', ))
    c = Process(target=f, args=(q,'c', ))
    d = Process(target=f, args=(q,'d', ))
    a.start();
    b.start();
    c.start();
    d.start();
    a.join();
    b.join();
    c.join();
    d.join();
    print(time.clock()-starttime)