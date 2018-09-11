


class region(object):

    def __init__(self, ID):
        self.ID = int(ID)
        self.dens=0
        self.diffPer= 0
        self.state=0
        self.numNbs=0
        self.nodelist=[]
        self.diffCells=0
        self.pos=[]

    def changeValue(self,value,dens):
        self.diffPer=value
        self.dens=int(dens)
        self.diffCells=value*dens

    def savePos(self,posx,posy):
        self.pos=[posx,posy]
        
    def add_nodes(self,node_IDs):
        self.nodelist=node_IDs
    def change_state(self,state):
        self.state=state

    def __hash__(self):
        return hash(self.ID)

    def __repr__(self):
        x=round(self.diffPer,2)
        return repr(self.ID)#+"\n"+repr(self.dens)
##        return repr(self.ID)+" , "+repr(x) #+"\n"+repr(self.dens)
##        return repr(x)
