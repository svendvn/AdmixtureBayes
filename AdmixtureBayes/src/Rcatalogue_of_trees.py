tree_clean={'s1':['s1s2',None, None, 0.1,None],
      's2':['s1s2', None, None, 0.1,None],
      's1s2':['r',None, None, 0.2,None],
      's3':['r',None, None, 0.2, None]}

tree_one_admixture={'s1':['s1b',None, None, 0.1,None],
      's1b':['s1s2','s3b',0.2, 0.1,0.2],
      's2':['s1s2', None, None, 0.1,None],
      's1s2':['r',None, None, 0.2,None],
      's3b':['r',None, None, 0.2, None],
      's3':['s3b',None,None,0.2,None]}

tree_two_admixture={'s1':['s1b',None, None, 0.1,None],
      's1c':['s1s2','s3b', 0.4,0.05,0.1],
      's1b':['s1c','s3a',0.2, 0.05,0.2],
      's2':['s1s2', None, None, 0.1,None],
      's1s2':['r',None, None, 0.2,None],
      's3b':['r',None, None, 0.2, None],
      's3':['s3a',None,None,0.1,None],
      's3a':['s3b', None,None,0.1,None]
      }

tree_two_admixture_cross={'s1':['s1b',None, None, 0.1,None],
      's1c':['s1s2','s3a', 0.4,0.05,0.1],
      's1b':['s1c',None,None, 0.05,None],
      's2':['s1s2', None, None, 0.1,None],
      's1s2':['r',None, None, 0.2,None],
      's3b':['r','s1b', 0.4, 0.2, 0.2],
      's3':['s3a',None,None,0.1,None],
      's3a':['s3b', None,None,0.1,None]
      }

tree_illegal={'s1':['s1b',None, None, 0.1,None],
      's1c':['s1s2','s3a', 0.4,0.05,0.1],
      's1b':['s1c','s3b',0.2, 0.05,0.2],
      's2':['s1s2', None, None, 0.1,None],
      's1s2':['r',None, None, 0.2,None],
      's3b':['r',None, None, 0.2, None],
      's3':['s3a',None,None,0.1,None],
      's3a':['s3b', None,None,0.1,None]
      }

tree_on_the_border={'s1':['c',None, None, 0.1,None],
      's2':['a',None, None,0.05,None],
      's3':['b',None,None, 0.3,None],
      'a':['b','d', 0.5,0.2,0.1],
      'c':['r','d',0.5,0.1,0.1],
      'd':['e',None,None,0.05,None],
      'b':['e',None,None,0.02,None],
      'e':['r',None,None,0.05,None]}

tree_on_the_border2={'s1':['d',None, None, 0.1,None],
      's2':['a',None, None,0.05,None],
      's3':['e',None,None, 0.3,None],
      's4':['b',None,None, 0.3,None],
      'a':['b','c', 0.5,0.2,0.1],
      'c':['e','d',0.5,0.1,0.1],
      'b':['f',None,None,0.05,None],
      'f':['r',None,None,0.02,None],
      'e':['f',None,None,0.05,None],
      'd':['r',None,None,0.05,None]}

tree_on_the_border2_with_children={'s1':['d',None, None, 0.1,None,None,None],
      's2':['a',None, None,0.05,None,None,None],
      's3':['e',None,None, 0.3,None,None,None],
      's4':['b',None,None, 0.3,None,None,None],
      'a':['b','c', 0.5,0.2,0.1,'s2',None],
      'c':['e','d',0.5,0.1,0.1,'a',None],
      'b':['f',None,None,0.05,None,'s4','a'],
      'f':['r',None,None,0.02,None,'b','e'],
      'e':['f',None,None,0.05,None,'c','s3'],
      'd':['r',None,None,0.05,None,'s1','c']}


tree_admix_to_child={
    's1':['r',None,None, 0.1,None],
    's2':['s2a',None,None,0.1,None],
    's3':['s3s2',None,None,0.1,None],
    's2a':['s3s2','s3s2a', 0.5,0.1,0.13],
    's3s2':['s3s2a',None,None,0.1,None],
    's3s2a':['r',None,None,0.01]
    }

tree_good={'s1':['d',None, None, 0.1,None,None,None],
      's2':['a',None, None,0.05,None,None,None],
      's3':['e',None,None, 0.3,None,None,None],
      's4':['b',None,None, 0.3,None,None,None],
      'a':['b','c', 0.5,0.2,0.1,'s2',None],
      'c':['e','d',0.5,0.1,0.1,'a',None],
      'b':['f',None,None,0.05,None,'s4','a'],
      'f':['r',None,None,0.02,None,'b','e'],
      'e':['f',None,None,0.05,None,'c','s3'],
      'd':['r',None,None,0.05,None,'s1','c']}

tree_without_consensus={'s1':['d',None, None, 0.1,None,None,None],
      's2':['a',None, None,0.05,None,None,None],
      's3':['e',None,None, 0.3,None,None,None],
      's4':['a',None,None, 0.3,None,None,None],
      'a':['b','c', 0.5,0.2,0.1,'s2',None],
      'c':['e','d',0.5,0.1,0.1,'a',None],
      'b':['f',None,None,0.05,None,'s4','a'],
      'f':['r',None,None,0.02,None,'b','e'],
      'e':['f',None,None,0.05,None,'c','s3'],
      'd':['r',None,None,0.05,None,'s1','c']}

tree_with_self_connection={'s1':['r',None, None, 0.1,None,None,None],
      's2':['a',None, None,0.05,None,None,None],
      's3':['f',None,None, 0.3,None,None,None],
      's4':['b',None,None, 0.3,None,None,None],
      'a':['b','c', 0.5,0.2,0.1,'s2',None],
      'c':['c',None, 0.5,0.1,0.1,'a','c'],
      'b':['f',None,None,0.05,None,'s4','a'],
      'f':['r',None,None,0.02,None,'b','s3']}

tree_with_pseudo_node={'s1':['d',None, None, 0.1,None,None,None],
      's2':['a',None, None,0.05,None,None,None],
      's3':['e',None,None, 0.3,None,None,None],
      's4':['b',None,None, 0.3,None,None,None],
      'a':['b','c', 0.5,0.2,0.1,'s2',None],
      'c':['e',None,None,0.1,None,'a',None],
      'b':['f',None,None,0.05,None,'s4','a'],
      'f':['r',None,None,0.02,None,'b','e'],
      'e':['f',None,None,0.05,None,'c','s3'],
      'd':['r',None,None,0.05,None,'s1',None]}

tree_with_doppel_band={'s1':['d',None, None, 0.1,None,None,None],
      's2':['a',None, None,0.05,None,None,None],
      's3':['d',None,None, 0.3,None,None,None],
      's4':['b',None,None, 0.3,None,None,None],
      'a':['b','c', 0.5,0.2,0.1,'s2',None],
      'c':['e','e',0.5,0.1,0.1,'a',None],
      'b':['f',None,None,0.05,None,'s4','a'],
      'f':['r',None,None,0.02,None,'b','e'],
      'e':['f',None,None,0.05,None,'c','c'],
      'd':['r',None,None,0.05,None,'s1','s3']}

tree_with_negative_bl={'s1':['d',None, None, 0.1,None,None,None],
      's2':['a',None, None,0.05,None,None,None],
      's3':['e',None,None, 0.3,None,None,None],
      's4':['b',None,None, 0.3,None,None,None],
      'a':['b','c', 0.5,-0.2,0.1,'s2',None],
      'c':['e','d',0.5,0.1,0.1,'a',None],
      'b':['f',None,None,0.05,None,'s4','a'],
      'f':['r',None,None,0.02,None,'b','e'],
      'e':['f',None,None,0.05,None,'c','s3'],
      'd':['r',None,None,0.05,None,'s1','c']}
        
tree_with_illegal_alpha={'s1':['d',None, None, 0.1,None,None,None],
      's2':['a',None, None,0.05,None,None,None],
      's3':['e',None,None, 0.3,None,None,None],
      's4':['b',None,None, 0.3,None,None,None],
      'a':['b','c', 1.5,0.2,0.1,'s2',None],
      'c':['e','d',0.5,0.1,0.1,'a',None],
      'b':['f',None,None,0.05,None,'s4','a'],
      'f':['r',None,None,0.02,None,'b','e'],
      'e':['f',None,None,0.05,None,'c','s3'],
      'd':['r',None,None,0.05,None,'s1','c']}