import logging

def addFileHandler(fn):
    ll = logging.getLogger("xmmpy")
    fh = logging.FileHandler(fn, 'w')
    fh.formatter = logging.Formatter(fmt="%(asctime)s [%(name)s - %(levelname)s] %(filename)s -> %(funcName)s :  %(message)s")
    ll.addHandler(fh)

def setup_logging(config):
    """
    """
    if isinstance(config, dict):
       conf = config 
    elif isinstance(config, str):
        conf = read_config(config)
    else:
        conf = {}
    
    if "LOGGING" in conf and "LOG_FILE" in conf["LOGGING"]:
        log_file = config["LOGGING"]["LOG_FILE"]
    else:
        log_file = "xmmpy.log"
        
    #print("lof", log_file)    
    #logger = logging.getLogger()

    logger = logging.getLogger("xmmpy")
    
    sh = logging.StreamHandler()
    sh.formatter = logging.Formatter(fmt="%(asctime)s [%(name)s - %(levelname)s] %(filename)s -> %(funcName)s :  %(message)s")
    logger.addHandler(sh)
    
    # Set logging level to the logger
    logger.setLevel(logging.DEBUG) # <-- THIS!
    logger.propagate = False
    #logger.info(80*"=")
    return logger
    
setup_logging(None)
