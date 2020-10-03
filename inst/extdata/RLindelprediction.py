import Lindel, os, sys
from Lindel.Predictor import *
import pickle as pkl


def write_file(seq,array,freq):
    sequences,frequency,indels = [],[],[]
    ss = 13
    sequences.append(seq[0:30] + ' | '+ seq[30:60])
    frequency.append('0')
    indels.append('')
    for i in range(len(array)):
        pt = array[i][0]
        try:
            idx1,dl = map(int,pt.split('+'))
            idx1 += ss+17
            idx2 = idx1 + dl
            cs = ss+17
            if idx1 < cs:
                if idx2>=cs:
                    s = seq[0:idx1]+'-'*(cs-idx1) + ' ' + '|' + ' ' + '-'*(idx2-cs)+seq[idx2:]
                else:
                    s = seq[0:idx1]+'-'*(idx2-idx1) + seq[idx2:cs]  + ' ' + '|' + ' ' + seq[cs:]
            elif idx1 > cs:
                s = seq[0:cs]+' ' + '|' + ' '+ seq[cs:idx1]+'-'*int(dl)+seq[idx2:]
            else:
                s = seq[0:idx1]+ ' ' + '|' + ' ' +'-'*int(dl)+seq[idx2:]
            indels.append('D' + str(dl) + '  ' +str(idx1-30))
        except ValueError:
            idx1 = int(pt.split('+')[0])
            if pt!='3':
                bp = pt.split('+')[1]
                il = str(idx1)
                indels.append('I' +il +'+' + bp)
            else:
                bp ='X' # label any insertion >= 3bp as X
                il = '>=3'
                indels.append('I3' + '+' + bp)
            s = seq[0:ss+17]+' '+bp+' '*(2-len(bp))+seq[ss+17:]
        sequences.append(s)
        frequency.append("{0:.8f}".format(freq[pt]*100))
    res = []
    for s,f,i in zip(sequences,frequency,indels):
        # res.append(s+'\t'+f + '\t'+i +'\n')
        res.append([s,f,i])
    return(res)


def rlp(inputseq):
    weights = pkl.load(open(os.path.join(Lindel.__path__[0], "Model_weights.pkl"),'rb'))
    prerequesites = pkl.load(open(os.path.join(Lindel.__path__[0],'model_prereq.pkl'),'rb'))
    seq = inputseq.upper() #input your sequence here
    try:
        y_hat, fs = gen_prediction(seq,weights,prerequesites)
        #filename += '_fs=' + str(round(fs,3))+'.txt'
        rev_index = prerequesites[1]
        pred_freq = {}
        for i in range(len(y_hat)):
            if y_hat[i]!=0:
                pred_freq[rev_index[i]] = y_hat[i]
        pred_sorted = sorted(pred_freq.items(), key=lambda kv: kv[1],reverse=True)
        res = write_file(seq,pred_sorted,pred_freq)
        return(res)
    except ValueError:
        return ('Error: No PAM sequence is identified.Please check your sequence and try again')
