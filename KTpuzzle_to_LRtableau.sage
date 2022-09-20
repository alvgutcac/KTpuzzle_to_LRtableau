##############################################################################
############## To add in sage/combinat/knutson_tao_puzzles.pyc  ##############
############## (either as methods of PuzzleFilling or as        ##############
############## standalone functions.)                           ##############
##############################################################################

def _NE_to_S_path(puzzle, coord):
    """
    Returns the content row of the LR skew tableau corresponding to ``coord``.
        
    This method traces out a path from ``coord`` to its "mirror" coordinate.
    If ``coord`` specifies the `i`-th 1 from the top on the north-east border
    of the puzzle, then its mirror is the `i`-th 1 from the west on the south
    border. The algorithm records the content numbers of the traced horizontal
    rhombi. See [Purbhoo07]_.
        
    .. WARNING::
        
        This method only works for the classical cohomology puzzle pieces.
        
    REFERENCES:
        
    .. [Purbhoo07] K. Purbhoo, Puzzles, Tableaux and Mosaics, {{{:arXiv:`0705.1184`}}}
        
    INPUT:
        
    - ``coord`` -- north-east boundary position counting from the top whose label is 1
        
    OUTPUT:
    
    - a list of numbers giving the content of one row in the LR tableau
        
    TESTS::
        
        sage: ps = KnutsonTaoPuzzleSolver("H")
        sage: solns = ps('0101011','0101101')
        sage: puzzle = solns[4]
        sage: puzzle._ne_labels
        ('0', '1', '0', '1', '1', '0', '1')
        sage: _NE_to_S_path(puzzle, 2)
        [2]
        sage: _NE_to_S_path(puzzle, 3)
        Traceback (most recent call last):
        ...
        AssertionError: the coordinate needs to be a coordinate of a 1 on the north east boundary
        
        sage: _NE_to_S_path(puzzle, 4)
        [1, 1]
    """
    assert puzzle._ne_labels[coord-1] == '1', "the coordinate needs to be a coordinate of a 1 on the north east boundary"
    i = coord; j = puzzle._n
    
    k = 0
    moving = "west"
    LR_list=[]
    while i != 0:
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
            if j-i <= 1:
                break
            current_piece = puzzle[(i,j)].south_piece()
            current_labels = current_piece.border()
            if current_labels == ('1', '1', '1'):
                moving = "west"
            j = j-1
        assert j > 0, "something went wrong and the path escaped the puzzle"        
    return LR_list
    
def KTpuzzle_to_LRtableau(puzzle):
    """
    Creates a skew Littlewood--Richardson tableau from a puzzle.
    
    We follow the bijection given in [Purbhoo07]_. A similar but different
    bijection is given in [Vakil03]_.
    
    INPUT:
    
    - ``puzzle`` -- a PuzzleFilling
    
    OUTPUT:
    
    - a LR skew tableau
    
    .. WARNING::
        
        This method only works for the classical cohomology puzzle pieces.
        
    REFERENCES:
    
    .. [Purbhoo07] K. Purbhoo, Puzzles, Tableaux and Mosaics, {{{:arXiv:`0705.1184`}}}
        
    .. [Vakil03] R. Vakil, A geometric Littlewood-Richardson rule, {{{:arXiv:`0302294`}}}
    
    EXAMPLES::
        
        sage: ps = KnutsonTaoPuzzleSolver("H")
        sage: solns = ps('010101','010101')
        sage: sorted([KTpuzzle_to_LRtableau(puzzle) for puzzle in solns])
        [[[None, None, None], [1, 1], [2]],
         [[None, None, 1], [None], [2]],
         [[None, None, 1], [None, 2], [1]],
         [[None, 1, 1], [None, 1], [None]]]

    TESTS::
    
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
    return SkewTableau(tab)
    
##############################################################################
############## To add in sage/combinat/partitions.py            ##############
##############################################################################
    
def abacus_to_partition(abacus):
    """
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
