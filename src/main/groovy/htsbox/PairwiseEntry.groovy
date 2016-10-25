package htsbox

/**
 * Created by nathandunn on 10/24/16.
 */
class PairwiseEntry {

//    In the output, each line gives QName, QLen, QStart, QEnd, Strand, RName, RLen, RStart, REnd, PerBaseDivergence, MapQ and semicolon-delimited misc information.

    String qName
    int qLen
    int qStart
    int qEnd
    int strand
    String rName
    int rLen
    int rStart
    int rEnd
    int perBaseDivergence
    Short mapQ
// and semicolon-delimited misc information.


}
