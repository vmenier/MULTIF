from pylab import *
import cPickle
import sys
from shaowuML import prepTrain, prepTest, trainAdaboost, testEnsemble

def read_fann_format(filename):
    with open(filename, "r") as f:
        content = f.readlines()
    header = content[0]
    features = content[1::2]
    beta = content[2::2]
    assert(len(beta) == len(features))
    n = len(beta)
    
    beta_array = np.zeros(n)
    p1 = np.zeros_like(beta_array)
    p2 = np.zeros_like(beta_array)
    p3 = np.zeros_like(beta_array)
    
    for i in range(n):
        beta_array[i] = float(beta[i])
        fl = features[i].replace("\n", "").split(",")
        p1[i] = float(fl[0])
        p2[i] = float(fl[1])
        p3[i] = float(fl[2])
    return p1, p2, p3, beta_array

# read training data
p1, p2, p3, beta = read_fann_format("train.dat")

# read features file
pid, p1_t, p2_t, p3_t = np.loadtxt("sample_features.dat", unpack=True)

config = {}
config['eps'] = 1e-6
config['verbose'] = False
config['cv_folds'] = 2
config['normalize_scheme'] = ['log-normal', 'log-normal','sqrt-normal']

config['adaboost'] = {'max_depth':16, 'loss':'square', 'lr':0.05, 'n_est':1500 }

##################### prepare data #####################################
train_feature = np.array([p1,p2,p3]).T
train_target = np.array(beta)
test_feature = np.array([p1_t, p2_t, p3_t]).T

##################### training & testing data ##########################
df_train_feature, normalize_dict, config, Y_train = prepTrain(train_feature, train_target, config)

ensemble_model = trainAdaboost( df_train_feature, Y_train, config )
del ensemble_model
del config
del normalize_dict
with open('config.pkl', 'rb') as fid:
    config = cPickle.load(fid)

with open('normalize_dict.pkl', 'rb') as fid:
    normalize_dict = cPickle.load(fid)

with open('ensemble_model.pkl', 'rb') as fid:
    ensemble_model = cPickle.load(fid)

df_test_feature = prepTest(test_feature, config, normalize_dict)
test_target = testEnsemble(ensemble_model, df_test_feature, config)

print test_target.shape
#show()


f = open("beta_for_su2.dat", "w")
for i in range(len(pid)):
    f.write("%i"%pid[i] + "\t%.12e\n"%(test_target[i]))
f.close()
#
#figure()
#plot(beta, test_target, "r.")
#plot(beta, beta, "g-", lw=2)
#savefig("tmp.png")

#figure()
#plot(p1, p2, 'r.')
#xlabel('p1')
#ylabel('p2')
#savefig("p1_p2.png")

#figure()
#plot(p1, p3, 'r.')
#xlabel('p1')
#ylabel('p3')
#savefig("p1_p3.png")

#figure()
#plot(p2, p3, 'r.')
#xlabel('p2')
#ylabel('p3')
#savefig("p2_p3.png")

#show()
