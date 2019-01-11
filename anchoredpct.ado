*! Anchored ranking
* Joao Pedro Azevedo - World Bank Group
* Paul Corral - World Bank Group 
cap prog drop anchoredpct
program define anchoredpct, eclass
	version 11.2
	syntax varlist(min=1 numeric) [if] [in] [aw pw fw] , anchor(varlist numeric max=1) [nq(integer 100) GENerate(string) VARprefix(string) replace decrease]
	
//Mark the estimation sample
marksample touse

	local vars `varlist' 
	local vars: list vars - anchor
	
	
	
	if ("`varprefix'" == "") & ("`replace'"=="") & ("`generate'"==""){
		display as error "You must specify a variable prefix to store results, a new variable to store results (generate()), or the replace option"
		error 119
	}
	local valid = 0
	if ("`varprefix'" != "") local valid = `valid' +1
	if ("`replace'"!="")     local valid = `valid' +1
	if ("`generate'"!="")    local valid = `valid' +1
	
	
	if (`valid'!=1){
		display as error "You must specify only one option, varprefix, generate or replace"
		error 119
	}
	
	if "`replace'"!=""{
		local reemplaza `anchor' `vars'
	}
	
	if "`varprefix'"!=""{
		local reemplaza `varprefix'_`anchor'
		foreach x of local vars{
			local reemplaza `reemplaza' `varprefix'_`x'
		}
		foreach x of local reemplaza{
			qui: gen `x' = .
			qui: lab var `x' "Anchored percentile rank to `anchor'"
		}
	}
	
	if ("`generate'"!=""){
		local num = 1
		foreach x of varlist `anchor' `vars'{
			qui: gen `generate'_`num' = . 
			lab var `generate'_`num' "Anchored percentile rank `x' to `anchor'"
			local reemplaza `reemplaza' `generate'_`num'
			
			local num = `num'+1
		}
	}
	
	if ("`decrease'" == "")    local decrease = 0
	else                       local decrease = 1
	
	

	
	//Weights
	local wvar : word 2 of `exp'
	if "`wvar'"==""{
	tempvar peso
	gen `peso'=1
	local wvar `peso'
	}
	
	mata: nq   = strtoreal(st_local("nq"))
	mata: x    = st_data(.,tokens("`anchor' `vars'"),"`touse'")
	mata: w    = st_data(.,tokens("`wvar'"),"`touse'")
	mata: aa   = _fpctile2(x[.,1], nq, w)
	mata: st_matrix("percs", aa)
	mata: x    = _anchor(x,aa,nq)
	mata: st_store(.,st_varindex(tokens(st_local("reemplaza"))),"`touse'",x)
	
	mat coln percs = "qtile_ptile" "max" "min"
	ereturn matrix qminmax = percs
		
end

mata
// Create matrix containing percentiles - X is a vector


function _fpctile2(real colvector X, real scalar nq, |real colvector w) {
	if (args()==2) w = J(rows(X),1,1)
	if (rows(X) < nq) {
		_error(3200, "Number of bins is more than the number of observations")
		exit(error(3200))       
	}
	data = runningsum(J(rows(X),1,1)), X, w
	_sort(data,(2,1))
	nq0 = quadsum(data[.,3])/nq 
	q = trunc((quadrunningsum(data[.,3]):/nq0):-0.0000000000001):+1 	
	data = data, q
	
	info= panelsetup(data,cols(data))
	
	minmax=J(rows(info),3,.)
	
	ri = rows(info)
	
	for(i=1; i<=ri; i++){
	if (i==1) mm=0
		minmax[i,.] =data[info[i,1],cols(data)],colmax(data[|info[i,1],2 \ info[i,2],2|]),mm 
	mm=minmax[i,2]
	}
		
	return(minmax)
} //Output is: quantile, max val of quantile, max value of previous quantile

function _anchor(real matrix x, real matrix aa, real scalar nq){
	x2=J(rows(x),cols(x),0)
	decrease = strtoreal(st_local("decrease"))

	ra = rows(aa)
	rx2 = rows(x2)
	cx2 = cols(x2)
	
	
	for(i=1;i<=ra;i++){
		if (i==1)           x2= x2+((x:<=aa[i,2]):*(x:>=aa[i,3]):*J(rx2,cx2,aa[i,1]))
		if ((i>1) & (i<ra)) x2=	x2+((x:<=aa[i,2]):*(x:>aa[i,3]):*J(rx2,cx2,aa[i,1]))
		if (i==ra)          x2= x2+((x:>aa[i,3]):*J(rx2,cx2,aa[i,1]))
	}
	
	x2=x2:/nq
	if (decrease==1)   x2 = 1:-x2
	return(x2)
}



end

	
	



