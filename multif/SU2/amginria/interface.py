import os

def amg_call(config):
    
    cmd = '';
    cmd = "amg -in %s -sol %s -p 2 \
         -c %f -hgrad %.2f -hmin %le -hmax %le -out %s \
        -itp  %s  -nordg " \
        % (config['mesh_in'], config['sol_in'],  \
        config['size'],  config['hgrad'], config['hmin'], config['hmax'], \
        config['mesh_out'], config['itp_sol_in']);
        
    if config['adap_source'] != "":
        cmd += ' -source %s ' % config['adap_source']
    
    if config['adap_back'] != "":
        cmd += ' -back %s ' % config['adap_back']
    
    cmd += ' > %s' % config['amg_log']
    os.system(cmd);