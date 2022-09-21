##############################################################################
############## To add in sage/combinat/knutson_tao_puzzles.py   ##############
############## (either as methods of PuzzleFilling or as        ##############
############## standalone functions.)                           ##############
##############################################################################

from sage.combinat.knutson_tao_puzzles import PuzzleFilling, H_grassmannian_pieces

# Method for equality on PuzzleFilling
# def __eq__(self, other):
    # return self._squares == other._squares

def _NE_to_S_path(puzzle, coord):
    r"""
    Returns the content row of the LR skew tableau corresponding to ``coord``.
        
    This method traces out a path from ``coord`` to its "mirror" coordinate.
    If ``coord`` specifies the `i`-th 1 from the top on the north-east border
    of the puzzle, then its mirror is the `i`-th 1 from the west on the south
    border. The algorithm records the content numbers of the traced horizontal
    rhombi. See [Purbhoo07]_.
        
    .. WARNING::
        
        This method only works for the classical cohomology puzzle pieces.
        
    REFERENCES:
        
    .. [Purbhoo07] K. Purbhoo, Puzzles, Tableaux and Mosaics, :arxiv:`0705.1184`
        
    INPUT:
        
    - ``coord`` -- north-east boundary position counting from the top whose label is 1
        
    OUTPUT:
        
    - a list of numbers giving the content of one row in the LR tableau
        
    TESTS::
        
        sage: ps = KnutsonTaoPuzzleSolver("H")
        sage: solns = ps('0101011','0101101')
        sage: puzzle = solns[4]
        sage: puzzle[(1,1)]
        0/1\10
        sage: _NE_to_S_path(puzzle, 2)
        [2]
        sage: _NE_to_S_path(puzzle, 3)
        Traceback (most recent call last):
        ...
        AssertionError: the coordinate needs to be a coordinate of a 1 on the north-east boundary
        
        sage: _NE_to_S_path(puzzle, 4)
        [1, 1]
        
    """
    assert puzzle._ne_labels[coord-1] == '1', "the coordinate needs to be a coordinate of a 1 on the north-east boundary"
    i = coord; j = puzzle._n
    
    k = 0
    moving = "west"
    LR_list=[]
    while j != 0:
        if moving == "west":
            if j-i == 0:
                current_piece = puzzle[(i,j)]
            else:
                current_piece = puzzle[(i,j)].north_piece()
            current_labels = current_piece.border()
            if current_labels == ('1', '1', '1'):
                moving = "south"
                LR_list = [i+1 for i in LR_list]
            elif current_labels == ('0', '10', '1'):
                i = i-1
                j = j-1
                LR_list.insert(0, k)
        elif moving == "south":
            if j-i < 1:
                break
            current_piece = puzzle[(i,j)].south_piece()
            current_labels = current_piece.border()
            if current_labels == ('1', '1', '1'):
                moving = "west"
            j = j-1
        assert j > 0, "something went wrong and the path escaped the puzzle"        
    return LR_list
    
def KTpuzzle_to_LRtableau(puzzle):
    r"""
    Creates a skew Littlewood--Richardson tableau from a puzzle.
    
    We follow the bijection given in [Purbhoo07]_. A similar but different
    bijection is given in [Vakil03]_.
    
    .. WARNING::
        
        This method only works for the classical cohomology puzzle pieces.
        
    REFERENCES:
    
    .. [Purbhoo07] K. Purbhoo, Puzzles, Tableaux and Mosaics, :arxiv:`0705.1184`
        
    .. [Vakil03] R. Vakil, A geometric Littlewood-Richardson rule, :arxiv:`0302294`
        
    INPUT:
    
    - ``puzzle`` -- a PuzzleFilling
    
    OUTPUT:
    
    - a LR skew tableau
    
    EXAMPLES::
        
        sage: ps = KnutsonTaoPuzzleSolver("H")
        sage: solns = ps('010101','010101')
        sage: sorted([KTpuzzle_to_LRtableau(puzzle) for puzzle in solns])
        [[[None, None, None], [1, 1], [2]],
         [[None, None, 1], [None, 1], [2]],
         [[None, None, 1], [None, 2], [1]],
         [[None, 1, 1], [None, 2], [None]]]
        sage: # Example from Purbhoo07 (runs slowly)
        ....: ps = KnutsonTaoPuzzleSolver('H')
        ....: solns = ps('00000010001000010100', '00000000100010010100')
        ....: puzzle = solns[168]
        ....: tab = KTpuzzle_to_LRtableau(puzzle); tab.pp()
          .  .  .  .  .  .  .  .  .  .  1  1  1  1
          .  .  .  .  .  .  .  .  1  1  2  2  2
          .  .  .  .  .  1  1  2  2  3  3
          .  1  1  2  2  3  4  4
    """
    south = puzzle.south_labels()
    ne = puzzle._ne_labels
    k = sum(int(i) for i in ne)
    lam = abacus_to_partition(south) + [0]*k
    tab = []
    for i in range(puzzle._n):
        row = []
        if ne[i] == '1':
            k -= 1
            row += [None] * lam[k]
            row += _NE_to_S_path(puzzle, i+1)
            tab.insert(0, row)
    while [] in tab:
        tab.remove([])
    return SkewTableau(tab)
    

def LRtableau_to_KTpuzzle(tableau, size=None):
    r"""
    Takes a Knutson--Tao puzzle and returns a skew LR tableau.
    
    This is the inverse of :func:`.KTpuzzle_to_LRtableau`.
    
    INPUT:
        
    - ``tableau`` -- a Littlewood--Richardson SkewTableau
    
    - ``size`` -- the size of the output Knutson--Tao puzzle (optional)
    
    TESTS::
    
        sage: # Example from Purbhoo07
        sage: tab = SkewTableau([
        ....: [None]*10 + [1]*4,
        ....: [None]*8 + [1]*2 + [2]*3,
        ....: [None]*5 + [1]*2 + [2]*2 + [3]*2,
        ....: [None] + [1]*2 + [2]*2 + [3] + [4]*2])
        sage: puzzle = LRtableau_to_KTpuzzle(tab, 20)
        sage: puzzle[(5,10)]
        1/\0  0\/1
        sage: ''.join(puzzle.south_labels())
        '01000010001001000000'
        sage: # puzzle.plot() # not tested
    
    EXAMPLES::
    
        sage: ps = KnutsonTaoPuzzleSolver("H")
        sage: puzzle = ps('01010','01001')[0]
        sage: tab = KTpuzzle_to_LRtableau(puzzle); tab
        [[None, 1, 1], [2]]
        sage: puzzle2 = LRtableau_to_KTpuzzle(tab)
        sage: puzzle._squares == puzzle2._squares # To change for puzzle == puzzle2 after implementing eq
        True
    """
    # Extract border of puzzle from tableau
    
    tab_w = tableau.weight()
    tab_out = tableau.outer_shape()
    tab_inn = tableau.inner_shape()
    L = len(tableau)
    
    if size == None:
        size = max(tab_w[0] + L, tab_out[0] + L)
    assert size >= max(tab_w[0] + L, tab_out[0] + L), "the puzzle size inputed is too small"
    
    n = size
    
    lam = tab_w + [0]*(L-len(tab_w))
    lam = partition_to_abacus(lam, n)[::-1]
    mu = [(n - L) - r for r in tab_out][::-1]
    mu = partition_to_abacus(mu, n)[::-1]
    nu = tab_inn; nu = nu + [0]*(L-len(tab_inn))
    nu = partition_to_abacus(nu, n)
    
    # Initialize puzzle
    
    puzzle = PuzzleFilling(lam, mu)
    
    # Find the locations of 1-triangles (blue by default in the plot)
    
    chosenCols = [i for i in [1..n] if nu[i-1]=='1']
    chosenRows = [i for i in [1..n] if lam[i-1]=='1']
    delta_blue_positions = []
    nabla_blue_positions = []
    for col in range(L):
        propagationRow = []
        propagationCol = []
        for row in range(col+1):        
            i = chosenRows[row]
            j = chosenCols[col-row]
            
            delta_blue_positions.append((j, i))
            
            k = tableau[L-col+row-1].count(row+1)
            nabla_blue_positions.append((j+k, i+k+1))
            
            propagationCol.insert(0, j+k)
            propagationRow.insert(0, i+k+1)        
        chosenRows = sorted(propagationRow) + chosenRows[col+1:]
        chosenCols = sorted(propagationCol) + chosenCols[col+1:] 
    
    # Create a dictionary of boundaries
    # (not all boundaries are in the dictionary)
    
    D = {(i,j) : {} for i in [1..n] for j in [i..n]}
    for (i,j) in delta_blue_positions:
        D[(i,j)]['north_west'] = '1'
        D[(i,j)]['north_east'] = '1'
        (a, b) = (i, j)
        while ((a-1, b) not in nabla_blue_positions and a>1):
            a -= 1
            D[(a,b)]['north_east'] = '0'
            D[(a,b)]['north_west'] = '1'
            D[(a,b)]['south_east'] = '1'
            D[(a,b)]['south_west'] = '0'
        (a, b) = (i, j)
        while ((a, b) not in nabla_blue_positions and b>a):
            b -= 1
            D[(a,b)]['north_west'] = '0'
            D[(a,b)]['north_east'] = '10'
    for (i,j) in nabla_blue_positions:
        if (i,j) in D.keys():
            D[(i,j)]['south_west'] = '1'
            D[(i,j)]['south_east'] = '1'
        (a, b) = (i, j)
        while ((a, b-1) not in delta_blue_positions and a>0):
            a -= 1; b -= 1
            D[(a,b)]['south_east'] = '10'
            D[(a,b)]['south_west'] = '1'
            D[(a+1,b)]['north_east'] = '1'
            D[(a+1,b)]['north_west'] = '10'
    for i in [1..n]:
        D[(1,i)]['north_west'] = lam[i-1]
        D[(i,i)]['south'] = nu[i-1]
        D[(i,n)]['north_east'] = mu[i-1]
    
    # Now fill piece by piece
    # (like KnutsonTaoPuzzleSolver._fill_puzzle_by_pieces
    #  but with the extra constraints given by the above)
        
    all_pieces = H_grassmannian_pieces()
    dirs = ('north_east', 'north_west', 'south_east', 'south_west')
    rhombi = sorted(all_pieces.rhombus_pieces(), key=lambda p : ''.join(p[dir] for dir in dirs))    
        # the sorting guarantees that the 0/\0 0\/0 piece is always
        # inserted if possible.
    dirs = ('north_east', 'north_west', 'south')
    triangles = sorted(all_pieces.delta_pieces(), key=lambda p : ''.join(p[dir] for dir in dirs))
    
    for (i,j) in D.keys():
        candidates = []
        if i == j:
            pieces = triangles
            dirs = ('north_east', 'north_west', 'south')
        else:
            pieces = rhombi
            dirs = ('north_east', 'north_west', 'south_east', 'south_west')
        # Check what pieces have the desired boundaries
        for piece in pieces:
            if all(piece[dir] == D[(i,j)][dir] for dir in dirs if dir in D[(i,j)].keys()):
                candidates.append(piece)
        # Place the smallest possible piece (the one with most 0s)
        puzzle._squares[(i,j)] = candidates[0]
        
    return puzzle
    
##############################################################################
############## To add in sage/combinat/partitions.py            ##############
##############################################################################
    
def abacus_to_partition(abacus):
    r"""
    Returns a partition from an abacus.
    
    An abacus is a function `w : \mathbb{Z} \to \{0,1\}` such that
    `w(n) = 0` for `n \ll 0` and `w(n) = 1` for `n \gg 0`. It is usually
    represented with an infinite tuple, e.g. `(...,0,0,0,1,0,0,1,0,1,1,...)`.
    Here, we use finite tuples, and assume everything to the left is a 0
    and everything to the right is a 1.
    
    An abacus determines a partition, via the following interpretation:
    a 1 in the abacus encodes a vertical line, and a 0 encodes a
    horizontal one. Reading from left to right, the abacus spells the
    outline of the Young diagram.
            
    INPUT:
            
    - ``abacus`` -- a tuple, a list or a string of 1's and 0's
            
    OUTPUT:
    
    - a list of weakly decreasing integers specifying the corresponding partition shape
            
    EXAMPLES::
            
        sage: abacus_to_partition('010101')
        [3, 2, 1]
        sage: abacus_to_partition([1,1,0,0,0])
        []
        sage: abacus_to_partition((0,0,0,1,1))
        [3, 3]
        sage: abacus_to_partition(('1','1','1','0','0','0'))
        []
        sage: abacus_to_partition(['0','0','0','1','1','1'])
        [3, 3, 3]
    """
    part = []
    n = len(abacus)
    k = 0
    for i in range(0,n):
        if abacus[i] == '1' or abacus[i] == 1:
            k=k+1
            part.insert(0, i+1-k)
        elif abacus[i] != '0' and abacus[i] != 0:
            raise ValueError('an abacus should be a tuple, list or string of 0s and 1s')
    return Partition(part)
    
def partition_to_abacus(lam, size = 0):
    r"""
    Return an abacus from a partition.
    
    This is the inverse to :func:`.abacus_to_partition`. Additionaly, if ``size``
    is given, the abacus will be of lenght at least ``size`` (by padding with
    0s on the right). The number of 1s in the abacus will be ``len(lam)``.
    
    INPUT:
    
    - ``lam`` -- Partition or decreasing list of nonnegative integers
    
    - ``size`` -- Integer (optional. Default: 0)
        
    OUTPUT:
    
    - a string of 0s and 1s.
        
    EXAMPLES::
    
        sage: partition_to_abacus([3,2,1])
        '010101'
        sage: partition_to_abacus([3,2,1], 10)
        '0101010000'
        sage: partition_to_abacus([3,3])
        '00011'
        sage: partition_to_abacus([0,0], 6)
        '110000'
        sage: partition_to_abacus([0,0])
        '11'
        sage: partition_to_abacus([2,2], 4)
        '0011'
        sage: partition_to_abacus([2,2], 5)
        '00110'
        sage: partition_to_abacus([2,2,0], 5)
        '10011'
  
    """
    L = len(lam)
    lam = lam+[0]
    abacus = [1] * L
    for i in range(L):
        abacus = abacus[:(L-i-1)] + [0]*(lam[i]-lam[i+1]) + abacus[(L-i-1):]
    return ''.join(map(str, abacus + [0]*(size-len(abacus))))  
