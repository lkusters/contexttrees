# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 11:42:43 2016

@author: Lieneke Kusters
"""

ALPHABET = "ACGT"
import numpy as np
import warnings
class contexttree:
    """ parent class for all the tree objects, the class object can be
    initialized with a depth and an input string for which the counts are stored
    """
    
    def __init__(self,depth,sequence=None):
        if not isinstance(depth,int):
            raise ValueError("Invalid maximum depth",depth)
        
        self._symbolcounts = None
        self._initialcontext = None
        self._sequencelength = 0
        self._maximumdepth = depth  
        if sequence == None:
            return None
        self._symbolcounts = self.__countsymbols(sequence)
        return None
        
    def __str__(self):
        return "contexttree of depth {0}, initcontext {1}".format(self._maximumdepth,self._initialcontext )+\
        "\n symbolcounts: \n {0}".format(str(self._symbolcounts))
    
    def __verifyinputsequence(self,sequence):
        """ verify that the sequence has only symbols in the alphabet
        """
        # alphabet in sequence
        no_valid_symbols = 0 # number of valid symbols
        for symbol in ALPHABET:
            no_valid_symbols += sequence.count(symbol)
        if not no_valid_symbols == len(sequence):
            # invalid
            raise ValueError(
                "Sequence has values that are not in alphabet: "+\
                "{0} valid of total {1} symbols".format(str(no_valid_symbols), \
                str(len(sequence))),ALPHABET
                )
                
    def __countsymbols(self,sequence):
        """ Count the symbol occurences for a context tree, for the contexts at 
            maximum depth and return dict()
    
        Keyword arguments:
        sequence:   (str) The sequence for which we count the symbolcounts
        
        Returns:
        counts: (dict): keys are occuring contexts (tuple), counts are 
                            symbol counts for symbols of alphabet given context
        """
        
        # first verify if input is valid
        self.__verifyinputsequence(sequence)
        if len(sequence)<=self._maximumdepth:
            # we need a sequence of at least length > self._maximumdepth
            warnings.warn("sequence length {0}, is too short, return None".\
                format(str(len(sequence)))
                )
            return None
    
        # Now prepare the data
        if self._symbolcounts == None:
            counts = dict()
        else:
            counts = self._symbolcounts
        
        keys = dict(zip(ALPHABET,range(len(ALPHABET)))) # Conversion table
        sequence = sequence.upper() # upper case
        
        # Special case, tree of depth 0
        if self._maximumdepth == 0:
            counts[''] = [sequence.count(ALPHABET[index]) for index in range(4)]
            
            return counts
        
        initcontext = ''
        for i in range(self._maximumdepth):
            initcontext+=sequence[i]   
        sequence = sequence[self._maximumdepth:]
        initcontext = initcontext[::-1]
        
        # we start with initial context initialcontext
        context = initcontext
        # now each next state is just a shift
        for symbol in sequence:
            if context in counts:
                counts[context][keys[symbol]] += 1
            else:
                counts[context] = [0 for sym in range(len(ALPHABET))]
                counts[context][keys[symbol]] += 1
            context = symbol+context[:-1]
        self._initialcontext  = initcontext
        self._sequencelength += len(sequence)
        return counts
    
    def updatesymbolcounts(self,sequence):
        """ Count the symbol occurences for a context tree, for the contexts at 
            maximum depth and return dict()
    
        Keyword arguments:
        sequence:   (str) The sequence for which we count the symbolcounts
        
        Returns:
        counts: (dict): keys are occuring contexts (tuple), counts are 
                            symbol counts for symbols of alphabet given context
        """
        
#        if self._symbolcounts == None:
#            raise ValueError(
#                "symbolcounts cannot be updated as there are none")
        self._symbolcounts = self.__countsymbols(sequence)
    def combinecontexttrees(self,tree):
        """ add the input tree to this tree
        """
    
        if not tree._maximumdepth == self._maximumdepth:
            raise ValueError("cannot combine the trees as their depth is "+\
            "not matching ",self._maximumdepth,tree._maximumdepth)
            
        counts = self._symbolcounts
        for key, val in tree._symbolcounts.items():
            if key in counts:
                counts[key] = [a+b for (a,b) in zip(counts[key],val)]
            else:
                counts[key] = val
        self._sequencelength += tree._sequencelength
    def copycontexttree(self):
        """ Make a copy of this tree and return it
        """
        
        tree = contexttree(self._maximumdepth)
        tree._symbolcounts = dict(self._symbolcounts)
        tree._sequencelength = self._sequencelength
        tree._initialcontext = self._initialcontext
        return tree
          
class fulltree(contexttree):
    """ inherits from contexttree, adds functionality of rates 
    probabilities and rates calculation
    """
    
    def __verifytreedephts(self,tree):
        """ verify that the input tree has the same depth as the source
        """
        if not tree._maximumdepth == self._maximumdepth:
            raise ValueError("cannot combine the trees as their depth is "+\
            "not matching ",self._maximumdepth,tree._maximumdepth)
    def __verifyfulltree(self,tree):
        """ verify that the compared tree is a fulltree of correct depth"""
        self.__verifytreedephts(tree)
        if not(type(tree)==type(self)):
            raise ValueError("cannot combine the trees as the tree type is "+\
            "not matching ",type(tree),type(self))
        
    def __counts2logprobs(self,counts):
        """ convert symbolcounts to probabilities
        """
        denum = np.log2(sum(counts)+len(ALPHABET)/2)
        logprobs = [np.log2(c+1/2) - denum for c in counts]
        return logprobs
    
    def getprobs(self):
        """ Compute the corresponding probabilities of the symbols in log2-space
        given the counts and return rself
        """
        rself = 0
        symbollogprobs = dict()
        for key,counts in self._symbolcounts.items():
            logprobs = self.__counts2logprobs(counts)
            symbollogprobs[key] = logprobs
            rself -= sum([a*b for a,b in zip(counts,logprobs)])
        return symbollogprobs,rself/self._sequencelength
    def getratetree(self,tree):
        """ estimate the achievable compression rate when applying this
        tree model for compression of the sequence corresponding to the input
        tree
        """
        self.__verifytreedephts(tree)
        symbolprobs,_ = self.getprobs()
        
        rate = 0
        for key,val in tree._symbolcounts.items():
            if key in symbolprobs:
                rate -= sum([a*b for a,b in zip(val,symbolprobs[key])])
            else:
                rate -= sum([a*-2 for a in val])
        return rate/tree._sequencelength
    def getratesequence(self,sequence):
        """ apply the model to a sequence and return the list of corresponding
        rates corresponding to each symbol in the sequence
        """
        
        raise ValueError(
                "Sorry this functionality has not been implemented yet")
    def getdivergence(self,tree):
        """ This function returns the estimated KL-Divergence of the probability
        distribution of the own tree in comparison to other tree
        
        We use D(q_z||p_x) ~ lim_{n-> inf} 1/n log2(q_z(Z)/p_x(Z))
        
        Here q_z is the probability distribution of this tree and p_x 
        corresponds to the other tree
        """
        
        
        self.__verifyfulltree(tree)
        _,rself = self.getprobs()
        rother = tree.getratetree(self)
        
        divergence = rother-rself
        return divergence
        
    def getdistance(self,tree):
        """ Use divergence as a distance metric, by estimating it in
        both directions """
        
        div1 = self.getdivergence(tree)
        div2 = tree.getdivergence(self)
        
        print(div1)
        print(div2)
        return (div1+div2)/2
    def copycontexttree(self):
        """ Make a copy of this tree and return it
        """
        
        tree = fulltree(self._maximumdepth)
        tree._symbolcounts = dict(self._symbolcounts)
        tree._sequencelength = self._sequencelength
        tree._initialcontext = self._initialcontext
        return tree
        
        
       
            

