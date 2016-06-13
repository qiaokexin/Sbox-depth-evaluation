try:
    import cPickle as pickle
except ImportError:
    import pickle

def find_expression(sbox):
    myfile=open('all_expr_hw8_within_depth4.5.txt','rb')
    all_expr=pickle.load(myfile)
    myfile.close()
    key=list(all_expr.keys())
    d=0
    flag=0
    print('Given S-box:')
    print(sbox)
    for i in range(4):
        bl=''
        for k in range(16):
            bl=bl+str((sbox[k]>>(3-i))& 0x1)
        if bl in all_expr:
            if i==0:
                print('a\' = ')
            if i==1:
                print('b\' = ')
            if i==2:
                print('c\' = ')
            if i==3:
                print('d\' = ')
            print(all_expr[bl])
            if all_expr[bl][1]>d:
                d=all_expr[bl][1]
        else:
            print('No appear!')
            flag=1    
    if flag:
        print('Can not detect depth!')
    else:
        print('Depth is '+str(d)+'.')
    return
##example: test depth of Sb0 in Midori:
find_expression([0xc,0xa,0xd,0x3,0xe,0xb,0xf,0x7,0x8,0x9,0x1,0x5,0x0,0x2,0x4,0x6])
