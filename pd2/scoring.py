"""Collection of classes that can score associations."""


class AvgMaxScore:
    """Compute average and max values for numbers.
    
    To use, initialize, then add score components with add(xxx).
    To extract max and averages, use properties max_score and avg_score.
    
    Caveat: adding values of 0 is ignored.
    For example, add(4) and add(0) gives avg_score of 4, not 2.
    """
    
    def __init__(self):
        self.sum = 0
        self.max_score = 0
        self.counter = 0
        
    def add(self, score, ignore_zero=True):        
        """Associate a numerical score to this instance.        
        Function performs eager tracking of maximal score.
        
        score -- number, new score to track
        ignore_zero -- logical, when true adding a zero score does 
            not change the counters.
        """
        
        if score == 0 and ignore_zero:
            return
        
        self.sum += score
        self.counter += 1    
        self.max_score = max(score, self.max_score)
        
    @property
    def avg_score(self):
        """Gives an average score."""

        if self.counter == 0:
            return 0
        
        return self.sum / self.counter


class PhenodigmScoring:
    """Computation of phenodigm scores based on phenotype lists."""
    
    def __init__(self, oomap):
        """Setup including a cross-ontology map.
        
        This object will determine how the scoring will be interpreted. 
        For example, when oomap is structured as oomap[HPs][MPs],
        the score functions can compute associations between
        mouse models (phenos2; 2nd index in oomap) and 
        disease (phenos1; 1st index in oomap).
        """
        self.oomap = oomap

    def max_scoring_phenotypes(self, phenos1, phenos2):
        """Produce an explanation for what contributes to non-zero score
        between phenos1 and phenos2."""
        
        result = []
        score = 0            
        
        # similar calculation to max_score, but keep track of what 
        # terms give the max association
        for p1 in phenos1:                        
            for p2 in phenos2:
                oldscore = score
                # note here p2 and p1 are flipped!
                score = self.valp1p2(score, p2, p1)
                if score > oldscore:
                    result = [p1, p2]
        if score == 0:
            return "", ""
        return result[0], result[1]            

    def all_scoring_phenotypes(self, phenos1, phenos2):
        """Produce an explanation for all phenotypes that contribute
        to a score between phenos1 and phenos2."""
        
        # set up subsets of phenos1 and phenos2 
        result1, result2 = set(), set()

        for p2 in phenos2:
            for p1 in phenos1:            
                p1p2 = self.oomap.get(p2, p1)
                if p1p2:
                    result1.add(p1)
                    result2.add(p2)                
            
        return ",".join(result1), ",".join(result2)
            
    def valp1p2(self, val, p1, p2):
        """Helper function, max between val0 and p1/p2 score."""
        
        p1p2 = self.oomap.get(p1, p2)
        if p1p2:
            val = max(val, p1p2.score)
        return val

    def max_score(self, phenos1, phenos2):
        """Get a max score based on two phenotype lists.
        
        phenos1 -- list of phenotypes from the first dimension in oomap
        phenos2 -- list of phenotypes from the second dimension in oomap
        """
        
        result = 0
        # Brute force approach to get max across pairs phenos1/phenos2
        for p1 in phenos1:                        
            for p2 in phenos2:
                result = self.valp1p2(result, p1, p2)
        return result

    def avg_score(self, phenos1, phenos2):
        """compute the phenodigm avgScore."""
                                
        p12len = len(phenos1)+len(phenos2)        
        if p12len == 0:
            return 0                
                    
        scoresum = 0                    
        for p1 in phenos1:
            # elegant approach is to re-use max_score
            # scoresum += self.max_score([p1], phenos2)
            # but a more efficient way (benefit using cProfile)            
            temp = 0
            for p2 in phenos2:
                temp = self.valp1p2(temp, p1, p2)                         
            scoresum += temp
                        
        for p2 in phenos2:
            # elegant approach re-uses max_score
            # scoresum += self.max_score(phenos1, [p2])
            # but a more efficient way (benefit using cProfile)
            temp = 0
            for p1 in phenos1:
                temp = self.valp1p2(temp, p1, p2)            
            scoresum += temp            
                                
        return scoresum / p12len

    def calc_score(self, phenos2, phenos1, ideal1, scorefun):
        """Generic perc score that uses a specified scoring function.
        
        Note order or phenos2 and phenos1 in the function declaration. 
        
        phenos2 - list of phenotypes, e.g. associated with a model 
        phenos1 - list of phenotypes, e.g. associated with a disease
        ideal1 - list of ideal phenotypes, e.g. for a disease
        scorefun - a scoring function, e.g. self.ps.avg_score 
                        
        Returns at 2-tuple (value in [0,100], value in [0, inf])
        The first value is the normalized score, the second is unnormalized    
        """
        
        a = scorefun(phenos1, phenos2)
        if a == 0:
            return 0, 0
        b = scorefun(ideal1, phenos2)
        if b == 0:
            return 100, a
        return min(100, 100*a/b), a

    def score_avg_perc_raw(self, phenos2, phenos1, ideal1):
        """Obtain 'avgPerc' and 'avg' scores for model/disease combination."""
        return self.calc_score(phenos2, phenos1, ideal1, 
                               scorefun=self.avg_score)

    def score_max_perc_raw(self, phenos2, phenos1, ideal1):
        """Obtain 'maxPerc' and 'max' scores for model/disease combination."""        
        return self.calc_score(phenos2, phenos1, ideal1, 
                               scorefun=self.max_score)

    def common_phenos(self, phenos2, phenos1):
        """Obtain a vector with common phenotypes."""
        
        p12len = len(phenos1)+len(phenos2)        
        if p12len == 0:
            return []                
                    
        common = set()                    
        for p1 in phenos1:
            tempmax, labmax = 0, ""
            for p2 in phenos2:
                p1p2 = self.oomap.get(p1, p2)
                if p1p2:
                    if p1p2.score > tempmax:
                        tempmax = p1p2.score
                        labmax = p1+"_"+p2
            if labmax != "":
                common.add(labmax)
                        
        for p2 in phenos2:
            tempmax, labmax = 0, ""
            for p1 in phenos1:
                p1p2 = self.oomap.get(p1, p2)
                if p1p2:
                    if p1p2.score > tempmax:
                        tempmax = p1p2.score
                        labmax = p1+"_"+p2
            if labmax != "":
                common.add(labmax)            
                                
        return list(common)

    def score(self, phenos2, phenos1, ideal1):
        """Obtain a set of 'phenodigm' scores for phenotypes lists.
         
        e.g. comparing a model and a disease.
        
        phenos2 -- list of phenotypes, e.g. MPs for a mouse model
        phenos1 -- list of phenotypes, e.g. HPs for a disease
        ideal1 - list of ideal phenotypes, e.g. HPs for an ideal disease
        that would be a good fit for the mouse model
        
        returns -- an array of four scores, using avg_perc, max_perc.         
        """
        
        # shortcut and definition, empty phenotypes give null score
        if len(phenos1) == 0 and len(phenos2) == 0:
            return 0, 0, 0, 0

        a = self.score_avg_perc_raw(phenos2, phenos1, ideal1)
        b = self.score_max_perc_raw(phenos2, phenos1, ideal1)                    

        return [a[0], a[1], b[0], b[1]]                                    

