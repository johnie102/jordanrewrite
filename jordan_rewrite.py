import random

# Convenience Methods
def dict_add(d,k,v):
    '''Adds values v to d[k] with a default of d[k]=0'''
    if k in d: d[k]+= v
    else: d[k] = v
def dict_remove_zeroes(d):
    '''Remove all keys of d if the value is very close to zero'''
    for k,v in list(d.items()):
        if abs(v)<0.000001:
            del d[k]
def dict_to_expr(d):
    '''takes in a dictionary where keys are strings,
    and values are numbers and prints it as a sum.
    e.g.: {"a": 1, "b": -1, "c": 2} -> a - b + 2c'''
    if not d:
        return "0"
    s = ""
    sort = sorted(d.keys())
    for term in sort:
        v = d[term]
        if abs(v-round(v))<0.00001:
            v = round(v)
        if v==1: s+= "+ " + term
        elif v==-1: s+= "- " + term
        else: s+= " {0:+} {1}".format(v, term)
    s = s.strip()
    if s[0] == "+": return s[1:].strip()
    return s


class Base(object):
    '''Base class implementing some
    boilerplate arithmetic functions'''
    def copy(self):
        return self
    def __add__(self, other):
        r = self.copy()
        r += other
        return r
    def __neg__(self):
        r = self.copy()
        return -1*r
    def __sub__(self, other):
        r = self + (-other)
        return r
    def __pow__(self, n):
        if not isinstance(n,int) or n<=0:
            raise Exception("Can only raise to a positive int")
        r = self
        for i in range(n-1):
            r = r*self
        return r
    def __hash__(self):
        return hash(str(self))
    def __eq__(self, other):
        return str(self)==str(other)
    def __repr__(self):
        return str(self)
    def to_latex(self):
        return str(self)
    

class JordanMonomialBase(Base):
    '''Base Class for representing the content of
    Jordan multiplication operators, e.g. 'ab' in T_{ab}'''
    def __mul__(self, other):
        return JMProduct(self, other)

    def __rmul__(self, other):
        if isinstance(other, (int,float)):
            return other * Words(self)
        if isinstance(other, JordanMonomialBase):
            return JMProduct(self, other)
        raise NotImplementedError

    def __add__(self, other):
        if isinstance(other, Words):
            other.add(1, (self,))
            return other
        return Words(self) + Words(other)

    def is_normalised(self):
        return True

class JMSingle(JordanMonomialBase):
    '''Used for representing T_a, and T_b, e.g. where
    the content of the operator is a single variable'''
    def __init__(self, variable_name):
        self._variable_name = variable_name

    def __str__(self):
        return self._variable_name

    def to_special(self):
        d = {self._variable_name: 1}
        return d

class JMProduct(JordanMonomialBase):
    '''Used for representing T_{LR} where L and R are
    other instances of JordanMonomialBase'''
    def __init__(self, L, R):
        if isinstance(L, str):
            L = JMSingle(L)
        if isinstance(R, str):
            R = JMSingle(R)

        self._L = L
        self._R = R

    def __str__(self):
        ltext = None
        rtext = None
        if isinstance(self._L, JMSingle):
            ltext = str(self._L)
        else:
            ltext = "({})".format(str(self._L))
        if isinstance(self._R, JMSingle):
            rtext = str(self._R)
        else:
            rtext = "({})".format(str(self._R))
        if ltext<rtext:
            return ltext+rtext
        else:
            return rtext+ltext

    def to_latex(self):
        s = str(self)
        if s=="aa": return "a^2"
        if s=="bb": return "b^2"
        return s

    def to_special(self):
        '''L=a, R=b, then returns 1/2ab + 1/2ba (or its dict equivalent)'''
        l = self._L.to_special()
        r = self._R.to_special()
        result = {}
        for a,t in l.items():
            for b,s in r.items(): #a and b are strings, a+b is str concat
                dict_add(result, a+b, t*s*0.5)
                dict_add(result, b+a, s*t*0.5)
        dict_remove_zeroes(result)
        return result

    def is_normalised(self):
        '''Wether the normalisation equation can be applied'''
        return (isinstance(self._L, JMSingle) and isinstance(self._R,JMSingle))

    def normalise_step(self):
        '''The product operator is of the form T_{a(bc)}. Apply the normalisation rewrite and return the result'''
        if not isinstance(self._R,JMProduct):
            L = self._R
            R = self._L
        else:
            L = self._L
            R = self._R
        
        a = L
        b = R._L
        c = R._R
        result = W(a,b*c) + W(b,a*c) + W(c,a*b) - W(b,a,c) - W(c,a,b)
        return result

a = JMSingle("a")  # T_a
b = JMSingle("b")  # T_b
ab = a*b  # T_{ab}
aa = a*a  # T_{a^2}
bb = b*b  # T_{b^2}

class Words(Base):
    '''Class that can represent linear combinations of words
    of JM terms like 2T_aT_{b(ab)} - 3T_b.
    A synonym of this class is defined below as W=Words,
    because the name is used many times.'''
    def __init__(self, *term):
        self.terms = {}
        self.is_unit = True #whether the class represent the identity, e.g. whether it is empty or not
        if term:
            self.terms[Words.normalise_word(term)] = 1
            self.is_unit=False

    def __len__(self):
        return len(self.terms)

    def copy(self):
        r = Words()
        r.terms = {Words.normalise_word(k):v for k,v in self.terms.copy().items()}
        r.is_unit = self.is_unit
        return r

    @staticmethod
    def normalise_word(word):
        ''''makes sure that T_{a^2} terms always appear before T_a terms (and the same for b)'''
        r = list(word)
        while True:
            for i in range(len(word)-1):
                if ((r[i] == a and r[i+1] == aa) or
                    (r[i] == b and r[i+1] == bb)):
                    tmp = r[i]
                    r[i] = r[i+1]
                    r[i+1] = tmp
                    break
            else: #not broken out of loop, so no more rewrites
                break
        return tuple(r)
        
    def __str__(self):
        d = {}
        for k in self.terms:
            d[self._format_word(k)] = self.terms[k]
        return dict_to_expr(d)

    def _format_word(self, word):
        return "[" + ','.join([str(item) for item in word]) + "]"

    def to_latex(self):
        '''Formats the expression suitable for LaTeX output'''
        d = {}
        for k in self.terms:
            d[self._format_word_latex(k)] = self.terms[k]
        return dict_to_expr(d)

    def _format_word_latex(self,word):
        return "".join(["L_{}".format(("{"+item.to_latex()+"}") if len(item.to_latex())!=1
                                      else item.to_latex()) for item in word])

    
    def add(self, scalar, word):
        '''Adds scalar amount of the specified word to the expression'''
        w = Words.normalise_word(word)
        dict_add(self.terms,w,scalar)
        if abs(self.terms[w])<0.0001:
            del self.terms[w]
    
    def __iadd__(self, other):
        if other.is_unit:
            raise Exception("adding unit Words to another Words")
        for word, scalar in other.terms.items():
            self.add(scalar, word)
        return self

    def __mul__(self, other):
        if isinstance(other, JordanMonomialBase):
            other = Words(other)
        if other.is_unit: return self
        if self.is_unit: return other
        result = Words()
        result.is_unit = False
        for w1, s1 in self.terms.items():
            for w2, s2 in other.terms.items():
                result.add(s1*s2, w1+w2)
        return result

    def scalar_mult(self, val):
        for term in self.terms: self.terms[term] *= val

    def __rmul__(self,other):
        if isinstance(other, (int,float)):
            r = self.copy()
            r.scalar_mult(other)
            return r
        raise NotImplementedError

    def normalise(self):
        '''Applies the normalisation equation until it can no longer be applied'''
        amount = 0
        while True:
            word, i = self._find_normalisable_word()
            if word==None:
                print("Normalisations applied: " + str(amount))
                return
            reduced = word[i].normalise_step()
            amount += 1

            val = self.terms[word]
            reduced.scalar_mult(val)
            del self.terms[word]
            newterm = Words(*word[:i]) * reduced * Words(*word[i+1:])
            self += newterm

    def _find_normalisable_word(self):
        for word in self.terms:
            for i in range(len(word)):
                if not word[i].is_normalised():
                    return (word, i)
        return None, None
    

    def get_random_rewrite(self):
        possibilities = []
        t = list(self.terms.keys())
        random.shuffle(t)
        random.shuffle(rewrites)
        for term in t:
            for A,B,C in rewrites:
                if A in term and B in term:
                    #find all the matches of (A,B) in term
                    matches = [i for i in range(len(term)-1) if (A,B)==term[i:i+2]]
                    if matches:
                        return (term, random.sample(matches,1)[0], C)           
        return None

    def do_rewrite_step(self, reduc):
        term,index,rewrite = reduc
        new = W(*term[:index]) * (self.terms[term]*rewrite) * W(*term[index+2:])
        del self.terms[term]
        self += new

    def do_random_rewrites(self, N=100, list_rewrites=False, silent=False):
        '''Does a specified amount of random rewrites. Returning the expression
        with the minimal amount of terms it has found in its random path.
        If list_rewrites is True it also outputs the rewrites it has done'''
        minterms = len(self.terms)
        best_so_far = self.copy()
        path = []
        bestpath = []
        for n in range(1,N):
            reduc = self.get_random_rewrite()
            if not reduc: #no reduction possible, so we are done
                if list_rewrites: return self, bestpath
                else: return self
            self.do_rewrite_step(reduc)
            path.append(reduc)
            if len(self.terms)< minterms:
                minterms = len(self.terms)
                if not silent: print("New minimal amount of terms: " + str(minterms))
                best_so_far = self.copy()
                bestpath = path.copy()
            if n%1000 == 0:
                if not silent: print("At iteration " + str(n))

        if list_rewrites: return best_so_far,bestpath
        else: return best_so_far
        
    def iterative_reduce(self, iterations=100):
        '''For a specified amount of iterations, do 100 random reductions and remember
         the shortest expression. Then repeat the process with this expression'''
        r = self.copy()
        shortest = len(self.terms)
        for i in range(iterations):
            r = r.do_random_rewrites(100)
            if len(r) < shortest:
                shortest = len(r)
                if len(r)==0: break
        return r

    def greedy_reduce(self):
        '''Tries 500 random reductions, and remembers the best one.
        It keeps doing this to reduce the expression as 'greedily' as possible.
        If it gets stuck, it tries to do 2 reductions at once, and then 3, etc.
        It returns the set of reductions it has found.
        '''
        path = []
        r = self.copy()
        r.normalise()
        temperature = 1
        while len(r.terms):
            reductions = []
            best = r
            length = len(r.terms)
            for i in range(500):
                r2 = r.copy()
                r3, reduc = r2.do_random_rewrites(temperature,True,True)
                if len(r3.terms)<length:
                    reductions = reduc
                    best = r3
                    length = len(r3.terms)
            if reductions:
                print("found new step. {!s} terms left".format(length))
                r = best
                path.extend(reductions)
                temperature = 1
            else:
                temperature += 1
            
        return path

W = Words

def Q(a,b=None):
    if b: return W(a,b) + W(b,a) - W(a*b)
    return 2*W(a,a) - W(a*a)

# list of tuples (A,B,C) with the understanding that AB
# may be interchanged with BA + C or BA with AB - C
rewrite_base = [
    (b, ab, 0.5*(W(bb,a) - W(a,bb))),
    (a, ab, 0.5*(W(aa,b) - W(b,aa))),
    (a, bb, 2*(W(ab,b) - W(b,ab))),
    (b, aa, 2*(W(ab,a) - W(a,ab))),
    (ab,aa, -2*(W(b,a,a,a)-W(a,a,a,b))-2*(W(aa,a,b)-W(b,aa,a))-2*(W(a,a,b,a)-W(a,b,a,a))),
    (ab,bb, -2*(W(a,b,b,b)-W(b,b,b,a))-2*(W(bb,b,a)-W(a,bb,b))-2*(W(b,b,a,b)-W(b,a,b,b))),
    (aa,bb, 4*(W(b,a,b,a) - W(a,b,a,b) + W(a,ab,b) - W(b,ab,a))),
    ]

# add reverse rewrites
r = []
for A,B,C in rewrite_base:
    r.append((B,A, -C))

rewrite_base.extend(r)
del r
rewrites = []
for A,B,C in rewrite_base:
    rewrites.append((A, B, W(B,A) + C))

def verify_rewrites():
    for A,B, l in rewrite_base:
        h = W(A,B)-W(B,A) - l
        print(to_special(h))
    

def to_special(l,name="c"):
    if not isinstance(l, W):
        l = W(l)
        
    result = {}
    for word in l.terms:
        w = list(reversed(word))
        r = {name:1}
        for t in w:
            r2 = {}
            for a, s1 in t.to_special().items():
                for k,v in r.items():
                    dict_add(r2,a+k,0.5*v*s1)
                    dict_add(r2,k+a,0.5*v*s1)
            r = r2
        for k,v in r.items():
            dict_add(result, k, v*l.terms[word])
            
    dict_remove_zeroes(result)
    return result


def fundamental_ident():
    l = l = Q(aa*b)+4*Q(a*ab) - 4*Q(a*ab,aa*b) # Q_{Q_a b}
    r = Q(a)*Q(b)*Q(a)
    h = l-r
    print(to_special(h))  # Verify h should indeed be reducable to zero
    return h

def equation7():
    l = 4*Q(a,b)**2
    r = 4*Q(ab) + 2*Q(aa,bb) - Q(a)*Q(b) - Q(b)*Q(a)
    h = l-r
    print(to_special(h))
    return h
def equation8(): # Equation (8)
    l = Q(a)*Q(b)*Q(a) + Q(aa)*Q(b) + Q(b)*Q(aa)\
        + 4*(Q(a,b)*Q(a,b)*Q(a)+Q(a,b)*Q(a)*Q(a,b) + Q(a)*Q(a,b)*Q(a,b))
    r = Q(aa*b)+4*Q(ab*a) + 4*Q(aa*b,ab*a) + 2*Q(aa*a,bb*a) + 4*Q(aa*a,ab*b)
    h = l-r
    print(to_special(h))
    return h
def theorem(): # Reducing this to zero proves the fundamental equality
    l = (2*(Q(aa,bb)*Q(a) + Q(a)*Q(aa,bb)) 
         + 4*(Q(ab)*Q(a) + Q(a)*Q(ab) + Q(a,b)*Q(a)*Q(a,b)))
    r = 2*Q(aa*b) + 8*Q(ab*a) + 2*Q(aa*a,bb*a) + 4*Q(aa*a,ab*b)
    h = l-r
    print(to_special(h))
    return h

def theorem2(): # Reducing this to zero proves the fundamental equality - but with less manual work.
    l = (Q(aa)*Q(b) + Q(b)*Q(aa) + 4*(Q(a,b)**2*Q(a)
        + Q(a,b)*Q(a)*Q(a,b) + Q(a)*Q(a,b)**2))
    r = 8*Q(aa*b,a*ab) + 2*Q(aa*a,bb*a) + 4*Q(aa*a,ab*b)
    h = l-r
    print(to_special(h))
    return h
