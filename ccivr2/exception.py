class CcivrFormatException(Exception):
    '''
    tell insufficiency of columns in the input data
    '''
    def __init__(self,list,msg=""):
        num = len(list)
        msg = f'error: missing {num} required column{(lambda x:"s" if x >= 2 else "")(num)}: {str(list)[1:-1]}'
        super().__init__(msg)

def triggerException(items):

    required = {'id','Chr','Start','End','Strand'}

    if not required <= set(items):
        absent = required - set(items)
        raise CcivrFormatException(absent)

def item_check(items):
    try:
        triggerException(items)
    except CcivrFormatException as e:
        print(e)
        exit()