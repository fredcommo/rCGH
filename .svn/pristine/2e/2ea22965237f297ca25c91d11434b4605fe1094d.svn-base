
#############################
## CLASS DEFINITION
#############################
setClass(
    'rCGH',
    representation(
        info = 'character',
        cnSet = 'data.frame',
        param = 'list',
        segTable = 'data.frame'
        ),
    prototype=prototype(
        info = NULL,
        cnSet = data.frame(),
        param = list(),
        segTable = data.frame()
        )
)

setClass('rCGH-Agilent', contains = 'rCGH')
setClass('rCGH-SNP6', contains = 'rCGH')
setClass('rCGH-cytoScan', contains = 'rCGH')


#############################
## SHOW METHOD FOR THIS CLASS
#############################
setMethod(
    'show',
    signature = 'rCGH',
    definition = function(object){
        d <- dim(object@cnSet)
        message('\nInstance of class ', class(object))
        message()
        message('Dataset with ', d[1], ' probes and ', d[2], ' columns.')
        message('Array information:')
        message()
        print(getInfo(object))
        message()
    }
)
