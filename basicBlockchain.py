# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 22:50:51 2018

@author: GAR
"""
import hashlib
import datetime as date

class myBlock:
    def __init__ (self,index,proof,data,timestamp='',previousHash=''):
        self.index = index;
        self.timestamp = timestamp; 
        self.proof = proof;
        self.data = data;
        self.previousHash = previousHash;
        self.hash = self.calculateHash();                                     

    def calculateHash(self):
        return hashlib.sha256(str(self.index) + str(self.timestamp) + str(self.proof) + str(self.data) + str(self.previousHash)).hexdigest()
    
    def mineBlock(self,difficulty):        
        while self.hash[0:difficulty] != '0' * difficulty:
            self.proof += 1;
            self.hash=self.calculateHash();
        print("Block mined: " + self.hash);    
        

    
class myBlockchain:
        def __init__ (self):
            self.chain = [self.createGenesisBlock()]
            self.difficulty=4;
        
        def createGenesisBlock(self):
            return myBlock(0,0,"Genesis data",date.datetime.now(),"0")
        
        def getLastBlock(self):
            return self.chain[len(self.chain)-1]
        
        def addBlock(self,newBlock):
            newBlock.previousHash = self.getLastBlock().hash
            newBlock.timestamp = date.datetime.now();
            newBlock.mineBlock(self.difficulty);
            newBlock.hash = newBlock.calculateHash()
            self.chain.append(newBlock)
            
        def isChainValid(self):
            for cont in range(1,len(self.chain)):

                currentBlock = self.chain[cont]
                previousBlock = self.chain[cont-1]

                if currentBlock.hash != currentBlock.calculateHash():
                    return False
                
                if currentBlock.previousHash != previousBlock.hash:
                    return False
            return True
                

garChain=myBlockchain()

for attr, value in garChain.chain[0].__dict__.iteritems():
    print ("{" + str(attr) + ": " +  str(value) + "}")
print("}")

print("Mining block...")
garChain.addBlock( myBlock(1,1,{'from': 'GAR', 'to': 'GAR', 'concept': 'car','value': 10} ) )
print("Mining block...")
garChain.addBlock( myBlock(2,2,{'from': 'GAR_1', 'to': 'GAR_2', 'concept': 'house','value': 20} ) )

for cont in range(len(garChain.chain)):
    print("{")
    for attr, value in garChain.chain[cont].__dict__.iteritems():
        print ("{" + str(attr) + ": " +  str(value) + "}")
    print("}")
    
print ("The current chain is: " + str(garChain.isChainValid()))
garChain.chain[1].data={'concept': 'car','value': 10};
garChain.chain[1].hash=garChain.chain[1].calculateHash()
print ("The current chain is: " + str(garChain.isChainValid()))