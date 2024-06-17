_skip_list=['pyqt','swig']
#_skip_list=['pyqt']
f = open("./requirements.txt",'r')
req=f.readlines()
f.close()

req=[n.strip() for n in req]
for r in req[:]:
    for s in _skip_list[:]:

        if s in r:
            req.remove(r)

f = open("./requirements.txt",'w')
for r in req:
    print(r,file=f)
f.close()