package vinamax

// FMM globals
var (
	root  Cell      // roots the entire FMM tree
	level [][]*Cell // for each level of the FMM tree: all cells on that level. Root = level 0

	// statistics:
	totalPartners int
	totalNear     int
	totalCells    int
)
