int search_Value( Value *query, DArray_Value *vals, double r, DArray_Value *res )
{
	int num;
	ValueIterator *it, *end;
	it = DArray_Value_start(vals);
	end = DArray_Value_end(vals);
	for ( num = 0; it != end; ++it ) {
		if(*it==query)
			continue;
		if ( (
			((query)->(x) - ((*it)->val)->(x))*((query)->(x) - ((*it)->val)->(x))
			+ ((query)->(y) - ((*it)->val)->(y))*((query)->(y) - ((*it)->val)->(y))
		) < r ) {
			DArray_Value_append(res, *it);
			++num; 
		}
	} 
	return num;
}

int search_Item( Value *query, DArray_Item *vals, double r, DArray_Item *res )
{
	int num;
	ItemIterator *it, *end;
	it = DArray_Item_start(vals);
	end = DArray_Item_end(vals);
	for ( num = 0; it != end; ++it ) {
		if ( (
			((query)->(x) - ((*it))->(x))*((query)->(x) - ((*it))->(x))
			+ ((query)->(y) - ((*it))->(y))*((query)->(y) - ((*it))->(y))
		) < r ) {
			DArray_Item_append(res, *it);
			++num;
		}
	}
	return num;
}
