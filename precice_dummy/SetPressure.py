myTargetData = 0
##def performAction (time,sourcedata, targetData):
def performAction(time,targetData):
    global myTargetData
    print targetData, type(targetData)
##  mySourceData = sourceData # store (reference to) sourceData for later use
    myTargetData = targetData # store (reference to) targetData for later use
    for i in range(targetData.size):
        targetData[i] = 100002
        i = i+1


