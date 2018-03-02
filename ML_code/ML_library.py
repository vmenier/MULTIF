# automatic data normalization machine learning library

# author:
# Shaowu Pan

# date:
# 02/11/2017, Michigan


# packages included
import pandas as pd
import numpy as np
from pandas.tools.plotting import scatter_matrix
import cPickle


# Note:
# - input features, target (1 column) must be numpy array columns, each sample is one row

# ----- content of functions-------

# 1. prepTrain
# 2. prepTest
# 3. trainAdaboost
# 4. testEnsemble

def prepTrain(feature, target, config):

# input list
#
# ----------
# 
# feature, target, config

# output list
# 
# ----------
# 
# df_new, normalize_dict

# 0. input configuration of training
    eps = config.get('eps',1e-6);
    normalize_scheme = config.get('normalize_scheme',None)

# 1. preprocessing on target
    df_beta = pd.DataFrame(target, columns=['beta'])
    df_beta_new = pd.DataFrame()
    mean_beta = df_beta['beta'].mean()
    std_beta = df_beta['beta'].std()
    
    ## show target distribution before normal    
    if config['verbose'] == True: df_beta.hist(column='beta', xlabelsize=12, xrot = 45, figsize=(14,8),bins = 1000) 

    # centralize target
    df_beta_new['beta'] = df_beta['beta'].apply(lambda x: (x - mean_beta)/std_beta)
    print 'std of target', std_beta
    print 'mean of target', mean_beta

    ## show target distribution after normalization
    if config['verbose'] == True: df_beta_new.hist(column='beta', xlabelsize=12, xrot = 45, figsize=(14,8),bins = 1000)

    config['target'] = {'mean': mean_beta, 'std':std_beta}

# 2. check normalized target distribution
    print df_beta_new.describe()

# 3. preprocessing on features
    X = feature
    Y = np.vstack((df_beta_new['beta'].values))
# X = np.vstack((p1,p2,p3))
# X = X.T
# Y = np.vstack((df_beta_new['beta'].values))

# # 3.1 random permutation the X,Y
    random_index = np.random.permutation(X.shape[0])
    X = X[random_index]
    Y = Y[random_index]

    nfeature=X.shape[1]
    col_names = [ 'p'+str(i+1) for i in range(nfeature) ]
    config['col_names']=col_names;

    # col_names = ['p1','p2','p3']
    df = pd.DataFrame(X,columns=col_names)

# # 3.1.1 show the histogram of features
    if config['verbose'] == True: df.hist(column=col_names, xlabelsize=12, xrot = 45, figsize=(14,8),bins = 1000)

# # 3.1.2 show the scatter matrix of features
    if config['verbose'] == True: scatter_matrix(df,alpha=0.05, figsize=(16,8),diagonal='kde')

# # 3.1.3 show the statistics of feature space
    print df.describe()

# # 3.2 normalize on feature space
    df_new = pd.DataFrame()
    normalize_dict = {}
    for i in range(nfeature):
        featureName = col_names[i]
        featureNormalizeScheme = normalize_scheme[i]

        if featureNormalizeScheme == 'normal':
            df_new[featureName] = df[featureName]
        elif featureNormalizeScheme == 'log-normal':
            df_new[featureName] = df[featureName].apply(lambda x: np.log10(x+eps))
        elif featureNormalizeScheme == 'sqrt-normal':
            df_new[featureName] = df[featureName].apply(lambda x: np.sqrt(x+eps))
        
        mean_p = df_new[featureName].mean()
        std_p  = df_new[featureName].std()
        normalize_dict[featureName] = {'mean':mean_p, 'std':std_p }

        df_new[featureName] = df_new[featureName].apply(lambda x: (x-mean_p)/std_p)

# # 3.3 display the normalized feature space for scatter matrix
    if config['verbose'] == True: scatter_matrix(df_new,alpha=0.05, figsize=(16,8),diagonal='kde')

# # 3.4 show the statistics for normalized dataset
    print df_new.describe()

# # 3.4.1 add target mean and std
    normalize_dict['target'] = { 'mean':mean_beta, 'std':std_beta }

    with open('normalize_dict.pkl', 'wb') as fid:
        cPickle.dump(normalize_dict, fid, protocol=cPickle.HIGHEST_PROTOCOL)

# # 3.5 return normalized data
    return df_new, normalize_dict, config, Y


def trainAdaboost(df_new,Y,config):
    import os 
    from sklearn.ensemble import AdaBoostRegressor
    from sklearn.cross_validation import KFold, train_test_split
    from sklearn.metrics import r2_score
    from sklearn.tree import DecisionTreeRegressor
    from matplotlib import pyplot as plt

    df_x_train = df_new[config['col_names']]
    df_y_train = Y;

    x=df_x_train.values;
    y=Y;
    kf = KFold(y.shape[0], n_folds=config['cv_folds'], shuffle=True)
    tmp_count = 0
# # consider an stack-average-ensemble modeling
    ensemble_model = []
    idd=0
    fh = open('error','w')
    fh.write('adaboost'+'\n')
    fh.write('square loss function'+'\n')
    for train_index, test_index in kf:
    
        model = AdaBoostRegressor(DecisionTreeRegressor(max_depth=config['adaboost']['max_depth'],max_features=len(config['col_names'])), n_estimators=config['adaboost']['n_est'], learning_rate=config['adaboost']['lr'], loss=config['adaboost']['loss']).fit(x[train_index],y[train_index]);
        ensemble_model.append(model)
        predictions = model.predict(x[test_index])
        actual = y[test_index]
    
        _predictions = model.predict(x[train_index])
        _actual = y[train_index]
        print 'cross validated R2, training =', r2_score(_actual, _predictions)
        print 'cross validated R2 = ', r2_score(actual, predictions)
    
        # write into file

        fh.write('---\n')
        fh.write('CV R2 : training = ' + str(r2_score(_actual, _predictions))+'\n')
        fh.write('CV R2 : validation = ' + str(r2_score(actual, predictions))+'\n')
    
        idd = idd + 1
        if not os.path.exists('cv_data'):
            os.makedirs('cv_data')
        np.savez_compressed('./cv_data/tp_ta_vp_va_'+str(idd),_predictions, _actual,predictions,actual)
    
        fig = plt.figure(figsize=(8,8))
        ax1 = fig.add_subplot(1,1,1)
        ax1.scatter(predictions,actual,c='r',marker='o',edgecolors='black')
        plt.axis('equal')
    
        ax1.plot(actual,actual)
        plt.xlabel(r'$\beta_{prediction}$')
        plt.ylabel(r'$\beta_{truth}$')
        ax1.set_xlim([-2,4])
        ax1.set_ylim([-2,4])
        np.savetxt("tmp_scatter_%i.dat"%tmp_count, np.c_[predictions, actual])
        tmp_count += 1
    fh.close()

    with open('config.pkl', 'wb') as fid:
        cPickle.dump(config, fid, protocol=cPickle.HIGHEST_PROTOCOL)

    with open('ensemble_model.pkl', 'wb') as fid:
        cPickle.dump(ensemble_model, fid, protocol=cPickle.HIGHEST_PROTOCOL)

    return ensemble_model




def prepTest(featureTest, config, normalize_dict):

    normalize_scheme = config['normalize_scheme'] 
    eps = config['eps']
 
    # normalize testing data
    # test features need to be numpy ndarray

    df_test = pd.DataFrame()
    df_test_new = pd.DataFrame()

    col_names = config['col_names'];
    nfeature = featureTest.shape[1]

    for i in range(nfeature): 
        featureName = col_names[i]
        featureNormalizeScheme = normalize_scheme[i]
 
        df_test[featureName] = featureTest[:,i]

        if featureNormalizeScheme == 'normal':
            df_test_new[featureName] = df_test[featureName]
        elif featureNormalizeScheme == 'log-normal':
            df_test_new[featureName] = df_test[featureName].apply(lambda x: np.log10(x+eps))
        elif featureNormalizeScheme == 'sqrt-normal':
            df_test_new[featureName] = df_test[featureName].apply(lambda x: np.sqrt(x+eps))

        mean_p = normalize_dict[featureName]['mean']
        std_p  = normalize_dict[featureName]['std']
        df_test_new[featureName] = df_test_new[featureName].apply(lambda x: (x-mean_p)/std_p)

    # # 1. report statistical information for testing features
    if False:
        df_test.describe()
        df_test_new.describe()
        # # 2. report PDF for testing features
        scatter_matrix(df_test_new,alpha=0.05, figsize=(16,8),diagonal='kde')
        # # 3. return normalized feature
    return df_test_new







def testEnsemble(ensemble_model, df_test_new, config):

    #filtering = config['filtering']

    # filtering?
    if False:
        
        # loading filtering variables
        gid_test = filtering['gid_test']
        x_test = filtering['x_test']
        p3_test = filtering['p3_test']
        d_max = filtering['d_max']
        x_max = filtering['x_max']
        x_min = filtering['x_min']
        beta_min = filtering['beta_min']
        beta_max = filtering['beta_max']
        beta_test = filtering['beta_test']    
    
        # loading testing data to test + test
        x = df_test_new.values
        mean_beta = config['target']['mean']
        std_beta = config['target']['std']
        
        ind = 0
        tmp = 0
        for model in ensemble_model:
            tmp = tmp + model.predict(x)
        beta_test = tmp/len(ensemble_model)
        ones_ = np.ones_like(gid_test)
        # transform beta_test back to physical scale
        beta_test = beta_test*std_beta + mean_beta

        n = len(beta_test)
        # filtering the final result
        for i in range(n):
            if p3_test[i] < d_max and x_test[i] > x_min and x_test[i] < x_max:
                # filter the beta (target)
                beta_test[i] = min(max(beta_test[i], beta_min), beta_max)
                ind = ind + 1
        else:
            beta_test[i] = 1.0

        print 'the remaining point after filtering: ', ind

    else:
        
        # loading testing data to test + test
        x = df_test_new.values
        n = x.shape[0]
        mean_beta = config['target']['mean']
        std_beta = config['target']['std']
        #print 'total testing data = ',n
        ind = 0
        tmp = 0
        for model in ensemble_model:
            tmp = tmp + model.predict(x)
        beta_test = tmp/len(ensemble_model)
        # ones_ = np.ones_like(gid_test)
        # transform beta_test back to physical scale
        beta_test = beta_test*std_beta + mean_beta





    ## saving test prediction data depends on how you want to do

    # savetxt("alpha_test.dat", np.c_[gid_test, x_test, y_test, z_test, beta_test, ones_, ones_], fmt="%i %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f")
    # savetxt("alpha_test_wfiltering.dat", np.c_[gid_test, x_test, y_test, z_test, ones_, ones_, beta_test], fmt="%i %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f")
    # savetxt("alpha_test_wfiltering.dat", np.c_[gid_test, x_test, y_test, z_test, beta_test, ones_, ones_], fmt="%i %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f")


    return beta_test















