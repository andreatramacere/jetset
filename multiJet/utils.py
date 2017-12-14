
def commands(obj):
    cmd_list= dir(obj)
    for item in cmd_list:
        if item[0]!='_':
            print item