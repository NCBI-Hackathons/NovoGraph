package htsbox

/**
 * Created by nathandunn on 10/24/16.
 */
class PairwiseAlignmentSet {

    List<PairwiseEntry> entries = new ArrayList<>()

    def addEntry(PairwiseEntry htsBoxEntry) {
        entries.add(htsBoxEntry)
    }

}
