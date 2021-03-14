def Model_$className( m_list, e_list, prn=$prnDefault):
    T = victoria.AuxMatrix(names="$names", prn=prn)
$connections
    T.End()
    # Split in Erlang series of length m
    T.SplitErlang( e_list, m_list)
    return T