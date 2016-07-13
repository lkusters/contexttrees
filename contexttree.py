# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 11:42:43 2016

@author: Lieneke Kusters
"""

ALPHABET = "ACGT"
import numpy as np
import warnings
import copy
class contexttree:
    """ parent class for all the tree objects, the class object can be
    initialized with a depth and an input string for which the counts are stored
    
    usually this is not called externally
    instead we use fulltree or maptree
    """
    
    def __init__(self,depth,sequence=None):
        if not isinstance(depth,int):
            raise ValueError("Invalid maximum depth",depth)
        
        self._initialcontext = []
        self._symbolcounts = None
        self._sequencelength = 0
        self._maximumdepth = depth 
        self._rself = None # achievable compression rate of full source tree
        
        if not sequence == None: # None should only occur when self.copy is used
            self.__countsymbols(sequence)
            
        
    def __str__(self):
        return "tree: {2} of depth {0}, initcontext {1}".format(self._maximumdepth,self._initialcontext,type(self) )+\
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
                
    def __verifytreedephts(self,tree):
        """ verify that the input tree has the same depth as the source
        """
        if not tree._maximumdepth == self._maximumdepth:
            raise ValueError("trees cannot interact, as their depth is "+\
            "not matching ",self._maximumdepth,tree._maximumdepth)
    
    def _verifysametype(self,tree):
        """ verify that both trees are same type and depth"""
        self.__verifytreedephts(tree)
        if not(type(tree)==type(self)):
            raise ValueError("both trees should be of same type",type(self),type(tree))
                
    def _counts2logprobs(self,counts):
        """ convert symbolcounts to probabilities
        """
        denum = np.log2(sum(counts)+len(ALPHABET)/2)
        logprobs = [np.log2(c+1/2) - denum for c in counts]
        return logprobs
        
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
        self._initialcontext += [initcontext]
        self._sequencelength += len(sequence)
        
        self._symbolcounts = counts
        self._rself = None # just became invalid in case it was set before
    
    def updatesymbolcounts(self,sequence):
        """ update symbol counts: call __countsymbols()
    
        Keyword arguments:
        sequence:   (str) The sequence for which we count the symbolcounts
        
        """
        
        if self._symbolcounts == None:
            warnings.warn("cannot update, since no symbols were counted "+\
            "yet, initializing with current input sequence instead")
        
        self.__countsymbols(sequence)
        
    def combine(self,tree):
        """ add the input tree to this tree
        """
    
        self.__verifytreedephts(tree)
        
        if self._symbolcounts == None and tree._symbolcounts == None:
            warnings.warn("combining with an empty tree")
            # do nothing
        elif tree._symbolcounts == None:
            warnings.warn("combining with an empty tree")
            # do nothing
        elif self._symbolcounts == None:
            warnings.warn("combining with an empty tree")
            # copy tree to self
            for attr in vars(tree):
                setattr(self, attr, getattr(tree, attr))
        else:
            for key, val in tree._symbolcounts.items():
                if key in self._symbolcounts:
                    self._symbolcounts[key] = [a+b for (a,b) in zip(self._symbolcounts[key],val)]
                else:
                    self._symbolcounts[key] = val
            self._sequencelength += tree._sequencelength
            self._initialcontext += [tree._initialcontext]
            self._rself = None # just became invalid
        
    def getcopy(self):
        """ Make a copy of this tree and return it
        """
        
        tree = contexttree(self._maximumdepth)
        for attr in vars(self):
            setattr(tree, attr, copy.deepcopy(getattr(self, attr)))
    
        return tree
                
          
class fulltree(contexttree):
    """ inherits from contexttree, adds functionality of rates 
    probabilities and rates calculation
    """
        
    
    def __getprobs(self):
        """ Compute the corresponding probabilities of the symbols in log2-space
        given the counts and return rself
        """
        rself = 0
        
        symbollogprobs = dict()
        for key,counts in self._symbolcounts.items():
            logprobs = self._counts2logprobs(counts)
            symbollogprobs[key] = logprobs
            rself -= sum([a*b for a,b in zip(counts,logprobs)])
        self._rself = rself/self._sequencelength
        return symbollogprobs
      
    def getrself(self):
        """ return rself or calculate rself if not calculated yet"""
        
        if self._rself == None:
            self.__getprobs()
            
        return self._rself
        
    def getratetree(self,tree):
        """ estimate the achievable compression rate when applying this
        tree model for compression of the sequence corresponding to the input
        contexttree
        """
        
        self._verifysametype(tree)
        symbolprobs = self.__getprobs()
        
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
         = Rother - Rself  (for sequence Z)
        
        Here 'tree' corresponds to the (model of the) input sequence Z
        Here q_z is the probability distribution of this tree and p_x 
        corresponds to the other tree
        """
        
        rother = self.getratetree(tree)
        rself = tree.getrself()
        
        divergence = rother-rself
        return divergence
        
    def getdistance(self,tree):
        """ Use divergence as a distance metric, by estimating it in
        both directions """
        
        div1 = self.getdivergence(tree)
        div2 = tree.getdivergence(self)
        
        return (div1+div2)/2
    def getcopy(self):
        """ Make a copy of this tree and return it
        """
        
        tree = fulltree(self._maximumdepth)
        for attr in vars(self):
            setattr(tree, attr, copy.deepcopy(getattr(self, attr)))
    
        return tree
    
""" =========================================== """
class maptree(contexttree):        
    """ this one is initialized by a contexttree (counts) """
      
    def getcountsdepth(self,depth):
        """ get counts corresponding to certain depth"""
        tree = self.__getcountsallnodes()
        tree2 = dict()
        for key,val in tree.items():
            if len(key)==depth:
                tree2[key] = [v for v in val]
        
        return tree2
        
    def getrself(self):
        """ return rself or calculate rself if not calculated yet"""
        
        if self._rself == None:
            self.__getprobs()
            
        return self._rself
    def getleafs(self):
        """ return the labels of the leafs"""
        
        if self._rself == None:
            """ not yet initialized"""
            self.__setleafs()
            
        return self._leafs
        
    def __getprobs(self):
        """ Compute the corresponding probabilities of the symbols in log2-space
        given the counts and leafs and return logprobs and counts in the leafs
        """
        if self._rself == None:
            """ not yet initialized"""
            self.__setleafs()
            
        rself = 0
        newleafs = [] # it is possible that we have some leafs that
                            # do actually not occur in the sequence
        
        symbollogprobs = dict()
        symbolcountsleafs = dict()
        symbolcounts = self.__getcountsallnodes()
        for key in self._leafs:
            if key in symbolcounts:
                counts = symbolcounts[key]
                logprobs = self._counts2logprobs(counts)
                symbollogprobs[key] = logprobs
                symbolcountsleafs[key] = [count for count in counts]
                rself -= sum([a*b for a,b in zip(counts,logprobs)])
                newleafs += [key]
                
        self._rself = rself/self._sequencelength
        self._leafs = newleafs
        return symbollogprobs,symbolcountsleafs
        
    def __setleafs(self):
        
        """ Construct the MAP model corresponding to the given symbol counts at 
        maximum depth of the tree, and also calculate achievable 
        compression rate on self (given the model)
        
        """
        # Construct the model corresponding to Counts
        allcounts = self.__getcountsallnodes() # counts at various depths
        treepe = self.__getpe(allcounts)
        self.__getpm(treepe)
        
    
    def __getcountsallnodes(self):
        """ Given the symbol counts at maximum depth, recover the counts at 
        decreasing depths until and including the root of the tree
        """
            
        tree = dict(self._symbolcounts)
        counts_previous_depth = dict(self._symbolcounts)
        
        for depth in reversed(range(1,self._maximumdepth+1)):
            counts_current_depth = dict()
            for key in counts_previous_depth:
                # calc predecessor:
                newkey = key[:-1]
                if newkey in counts_current_depth:
                    counts_current_depth[newkey] =  [sum(x) for \
                        x in zip(counts_current_depth[newkey],
                        counts_previous_depth[key])
                        ]
                else:
                    counts_current_depth[newkey] = counts_previous_depth[key]
                    
            counts_previous_depth = counts_current_depth
            tree.update(counts_current_depth)
        return tree
    def __getpe(self,tree):
        """ Given the symbol counts at various depths, calculate the memoryless
        probabilities (in log2-space) of the corresponding sequences 
        using the KT estimator
        
        KT-estimate is defined as:
        Pe := Prod_{foreach symbol in ALPHABET}( (symbol counts-1/2)! ) / ..
           ( ( total symbol counts - len(ALPHABET)/2 )! ) 
        
        Keyword arguments:
        tree:   (dict) : keys are occuring contexts (str), counts are 
                            symbol counts for symbols of alphabet given context
        
        Returns:
        tree_pe: (dict): keys are occuring contexts (str), values are the 
            memoryless probabilities of the sequence corresponding to this context
            / we define log_2(0) = 0
        """
        
        treepe = dict()
        for context,vals in tree.items():
            lengthsubseq = sum(vals)
            if lengthsubseq > 0:
                # KT - estimator
                denum = np.log2(np.fromiter(range(lengthsubseq), float)+len(ALPHABET)/2).sum()
                numer = 0
                for x in vals:
                    if x>0:
                        numer+=np.log2(np.fromiter(range(1,x+1), float)-1/2).sum()
                treepe[context] = numer-denum
            else:
                treepe[context] = 0            
        return treepe
    
    def __getpm(self,treepe):
        """ Given the the memoryless probabilities (in log2-space) of the
        sequences corresponding to the different contexts, find the maximum a
        posteriori probability for each context (= node)
        
        maximum a posteriori probability is defined as:
        Pm := max( alpha* Pe ; (1-alpha)*Prod_{children} Pe)
        where alpha = (len(ALPHABET)-1)/len(ALPHABET)
        
        Keyword arguments:
        treepe:   (dict) : keys are occuring contexts (str), values are the 
            memoryless probabilities (in log2-space) of the sequences 
            corresponding to the different contexts
        
        Calculates:
        treepm: (dict) keys are occuring contexts (str), values are the 
            maximum a posteriori probabilities of the sequence corresponding to 
            this context
            / we defined log_2(0) = 0
        
        Returns:
        leafs: (list): tuples that correspond to the selected leafs
        """
        # define the penalties
        alpha = np.log2(len(ALPHABET)-1)-np.log2(len(ALPHABET))
        alpha_inv = -np.log2(len(ALPHABET))
        
        branches = [] # 
        tree_pm = dict()
        # from max depth to less depth
        for ctxt in sorted(treepe, key=len, reverse=True):
            if len(ctxt)== self._maximumdepth :
                # this is initialization
                tree_pm[ctxt] = treepe[ctxt]
            else:
                # calculate sum of children Pm
                child_sum = 0
                for symbol in ALPHABET:
                    child = ctxt+symbol
                    children = []
                    if child in tree_pm:
                        child_sum += tree_pm[child]
                        children.append(child)
                # compare the values
                if (alpha + treepe[ctxt]) >= (alpha_inv+child_sum):   
                    tree_pm[ctxt] = alpha + treepe[ctxt]
                    # this may be a leaf
                else:
                    tree_pm[ctxt] = alpha_inv+child_sum
                    # this may be a branch
                    branches.append(ctxt)
        # Now start at root and select leaf if in leafs
        self._leafs = self.__findleafs(branches)
        
    def __findleafs(self,branches,key=''):
        """ This recursive formula returns the leafs of the tree, given a list of
        branches. 
        
        Returns key if key is not a branch, else returns the result of self on the 
        children of key.
        
        Keyword arguments:
        branches:   (list) these are possible branches in the tree. (nodes that 
        need memory) 
        key: (default None), this parameter is used to recursively increase the 
        depth of the processed node
        
        Returns:
        leafs: (list): tuples that correspond to the selected leafs
        """
        leafs = []
        if key in branches:
            # then look at children
            for symbol in ALPHABET:
                leafs += self.__findleafs(branches,key+symbol)
        else:
            leafs.append(key)
        return leafs
       
            

