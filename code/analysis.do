*==============================================================================
* analysis.do
* Automation and Inequality within Cities: Evidence from Industrial Robots in China
*
* Runs every regression reported in the paper, from the files in data/.
* All results are written to output/.
*
* How to run: set Stata's working directory to the repository root, then
*     do code/analysis.do
* Requires Stata 16 or later and the user-written packages
*     ivreg2, ranktest, reghdfe, ivreghdfe, ftools, estout
* (install with: ssc install <name>).
*==============================================================================
clear all
set more off
global ROOT "`c(pwd)'"
global DATA "$ROOT/data"
global OUT  "$ROOT/output"
cap mkdir "$OUT"
cap log close
log using "$OUT/analysis.log", text replace
graph set window fontface "Times New Roman"

*------------------------------------------------------------------------------
* 0. Benchmark cross-city 2SLS, 2010-2020 (Table 3)
*------------------------------------------------------------------------------
use "$DATA/init.dta", clear
egen cid = group(id)
global X for_invest_gdp sec_indus_gdp thd_indus_gdp indus_etp_asst for_dir_invest

eststo clear
eststo o1: ivreg2 gini  $X (x = bartik_eu bartik_oth) i.year, cluster(cid)
* years that enter the regressions with controls (the controls end in 2016)
tab year if e(sample)
eststo o2: ivreg2 theil $X (x = bartik_eu bartik_oth) i.year, cluster(cid)
eststo o3: ivreg2 gini  (x = bartik_eu bartik_oth) i.year, cluster(cid)
eststo o4: ivreg2 theil (x = bartik_eu bartik_oth) i.year, cluster(cid)
esttab o1 o2 o3 o4 using "$OUT/T3_benchmark.csv", replace keep(x) ///
    b(%7.4f) se(%7.4f) star(* 0.1 ** 0.05 *** 0.01) stats(N N_clust widstat jp) csv

*------------------------------------------------------------------------------
* 1. Build the city panel, 1992-2023
*------------------------------------------------------------------------------
* 1a College density (college students per 10,000 residents), base value = 2009-2010 mean
use "$DATA/city_covariates.dta", clear
keep city_id year coll_grad
replace coll_grad = . if year == 2013     // reporting break: mean jumps from about 170 to 268
tempfile cg
save `cg'
keep if inlist(year, 2009, 2010)
collapse (mean) edu0 = coll_grad, by(city_id)
tempfile edu0
save `edu0'

* 1b Nighttime-light inequality panel
use "$DATA/city_inequality.dta", clear
gen city_id = int(bc/100)
replace city_id = 11 if city_id == 1100
replace city_id = 12 if city_id == 1200
replace city_id = 31 if city_id == 3100
replace city_id = 50 if city_id == 5000
* county-level units directly under Hainan, Hubei and Xinjiang map to one code
* and have no robot data, so they are dropped
duplicates tag city_id year, gen(dup)
drop if dup > 0
keep city_id year gini theil atkinson

merge m:1 city_id using "$DATA/city_bartik_chn_rob.dta", keepusing(city_chn00_16_06) keep(match) nogen
merge m:1 city_id using "$DATA/city_bartik_eu_rob.dta",  keepusing(city_eu00_16_06)  keep(match) nogen
merge m:1 city_id using "$DATA/city_bartik_oth_rob.dta", keepusing(city_oth00_16_06) keep(match) nogen
merge m:1 city_id using "$DATA/city_share_all.dta", keepusing(share_c00) keep(master match) nogen
merge m:1 city_id using "$DATA/controls_99_03.dta", ///
    keepusing(gdpsh_ser_03 ave_wage_03 tot_pass_03 ict_emp_03 rd_ratio_03 for_share_03 water_int_03) ///
    keep(master match) nogen
merge m:1 city_id using `edu0', keep(master match) nogen
merge 1:1 city_id year using `cg', keep(master match) nogen
rename (city_chn00_16_06 city_eu00_16_06 city_oth00_16_06) (rob_chn rob_eu rob_oth)
drop if rob_chn == .

gen prov = cond(city_id < 100, city_id, int(city_id/100))
gen lwage03 = ln(ave_wage_03)
gen lpass03 = ln(tot_pass_03)
gen lict03  = ln(ict_emp_03)
bysort city_id: egen gini0  = mean(cond(inrange(year, 1992, 1994), gini, .))
bysort city_id: egen theil0 = mean(cond(inrange(year, 1992, 1994), theil, .))
bysort city_id: egen atk0   = mean(cond(inrange(year, 1992, 1994), atkinson, .))
* HighEdu: base-period college density above the median across cities
preserve
keep city_id edu0
duplicates drop
qui sum edu0, detail
local edumed = r(p50)
restore
gen highedu = edu0 > `edumed' if edu0 < .

label var rob_chn "Robot exposure, China (2006-2016, 2000 shares)"
label var rob_eu  "Robot exposure, 6 EU countries"
label var rob_oth "Robot exposure, 5 other EU countries"
compress
save "$OUT/panel_city_1992_2023.dta", replace

* city-level descriptive statistics
preserve
keep if year == 2010
sum rob_chn rob_eu rob_oth share_c00 edu0
corr rob_chn rob_eu rob_oth share_c00
restore
tabstat gini theil atkinson, by(year) stat(mean sd n) format(%6.3f)

*------------------------------------------------------------------------------
* 2. Year-by-year cross-sectional 2SLS (Figure 2)
*------------------------------------------------------------------------------
tempname h
postfile `h' year b se using "$OUT/yearly_cross_section.dta", replace
forvalues y = 1992/2023 {
    qui ivreg2 gini (rob_chn = rob_eu rob_oth) i.prov if year == `y', robust partial(i.prov)
    post `h' (`y') (_b[rob_chn]) (_se[rob_chn])
}
postclose `h'

preserve
use "$OUT/yearly_cross_section.dta", clear
gen lo = b - 1.96*se
gen hi = b + 1.96*se
list, clean noobs
export delimited using "$OUT/yearly_cross_section.csv", replace
drop if year == 2022                   // abnormal year in the light series
twoway (rarea lo hi year, color(gs13)) (connected b year, mcolor(black) lcolor(black) msize(small)), ///
    yline(-0.1444, lpattern(dash) lcolor(black)) yline(0, lcolor(gs8)) ///
    xline(2013.5, lpattern(shortdash) lcolor(gs6)) ///
    xlabel(1992(4)2024) ylabel(-0.3(0.05)0.05, format(%4.2f)) ///
    xtitle("Year") ytitle("Coefficient on robot exposure") ///
    legend(off) graphregion(color(white)) plotregion(lcolor(black))
graph export "$OUT/F2_yearly_cross_section.png", replace width(2000)
restore

*------------------------------------------------------------------------------
* 3. Event study with city and province-by-year fixed effects, 2008 = 0 (Figure 3)
*    2022 is excluded because the average Gini falls by about a third that year
*------------------------------------------------------------------------------
forvalues y = 1992/2023 {
    if `y' != 2008 & `y' != 2022 gen ev`y' = rob_chn * (year == `y')
}
tempname e
postfile `e' str8 outcome year b se using "$OUT/event_study.dta", replace
foreach v in gini theil atkinson {
    reghdfe `v' ev* if year != 2022, absorb(city_id prov#year) cluster(city_id)
    forvalues y = 1992/2023 {
        if `y' == 2008 post `e' ("`v'") (`y') (0) (0)
        else if `y' != 2022 post `e' ("`v'") (`y') (_b[ev`y']) (_se[ev`y'])
    }
}
postclose `e'
drop ev*

preserve
use "$OUT/event_study.dta", clear
gen lo = b - 1.96*se
gen hi = b + 1.96*se
export delimited using "$OUT/event_study.csv", replace
local i = 0
foreach v in gini theil atkinson {
    local ++i
    local ttl : word `i' of "Gini" "Theil" "Atkinson"
    twoway (rcap lo hi year if outcome == "`v'", lcolor(gs9)) ///
           (scatter b year if outcome == "`v'", mcolor(black) msize(small)), ///
        yline(0, lcolor(gs8)) xline(2008, lpattern(dash) lcolor(gs6)) ///
        xline(2013.5, lpattern(shortdash) lcolor(gs6)) ///
        xlabel(1992(4)2024, labsize(small)) ylabel(, labsize(small) format(%4.2f)) ///
        title("`ttl'", size(medium)) xtitle("") ytitle("Coefficient", size(small)) ///
        legend(off) graphregion(color(white)) plotregion(lcolor(black)) name(es_`v', replace)
}
graph combine es_gini es_theil es_atkinson, cols(1) xsize(6) ysize(8) graphregion(color(white))
graph export "$OUT/F3_event_study.png", replace width(1800)
restore

*------------------------------------------------------------------------------
* 4. Within-city 2SLS, equation (3) (Table 4)
*    Panel A, main window 2000-2013: 2010-2013 against 2000-2009
*    Panel B, placebo window 1992-2004: 2000-2004 against 1992-1999
*    Panel C, later window 2014-2023 excluding 2022: exposure x linear trend
*------------------------------------------------------------------------------
global B03 gdpsh_ser_03 lwage03 lpass03 lict03 rd_ratio_03 for_share_03 water_int_03

cap program drop did_iv
program define did_iv
    * arguments: outcome, initial-inequality variable, sample condition, period variable, model prefix
    args y y0 cond per tag
    cap drop rp rp_eu rp_oth cx_*
    gen rp     = rob_chn * `per'
    gen rp_eu  = rob_eu  * `per'
    gen rp_oth = rob_oth * `per'
    eststo `tag'a: ivreghdfe `y' (rp = rp_eu rp_oth) if `cond', ///
        absorb(city_id prov#year) cluster(city_id)
    foreach c of global B03 {
        gen cx_`c' = `c' * `per'
    }
    gen cx_init = `y0' * `per'
    eststo `tag'b: ivreghdfe `y' cx_* (rp = rp_eu rp_oth) if `cond', ///
        absorb(city_id prov#year) cluster(city_id)
end

gen post_dmsp = inrange(year, 2010, 2013)
gen post_pre  = inrange(year, 2000, 2004)
gen trend_v   = year - 2014

eststo clear
local i = 0
foreach v in gini theil atkinson {
    local ++i
    local v0 : word `i' of gini0 theil0 atk0
    did_iv `v' `v0' "inrange(year,2000,2013)"            post_dmsp A`i'
    did_iv `v' `v0' "year<=2004"                         post_pre  P`i'
    did_iv `v' `v0' "year>=2014 & year!=2022"            trend_v   V`i'
}
foreach p in A P V {
    esttab `p'1a `p'1b `p'2a `p'2b `p'3a `p'3b using "$OUT/T4_panel_`p'.csv", replace keep(rp) ///
        b(%7.4f) se(%7.4f) star(* 0.1 ** 0.05 *** 0.01) stats(N N_clust widstat jp) csv
    esttab `p'1a `p'1b `p'2a `p'2b `p'3a `p'3b, keep(rp) b(%7.4f) se(%7.4f) ///
        star(* 0.1 ** 0.05 *** 0.01) stats(N N_clust widstat jp) title("Panel `p'")
}

* Robustness (Section 4.4): OLS, province clustering, later window cut at 2019
cap drop rp
gen rp = rob_chn * post_dmsp
eststo R1: reghdfe gini rp if inrange(year,2000,2013), absorb(city_id prov#year) cluster(city_id)
eststo R2: reghdfe gini rp if inrange(year,2000,2013), absorb(city_id prov#year) cluster(prov)
replace rp = rob_chn * trend_v
eststo R3: reghdfe gini rp if inrange(year,2014,2019), absorb(city_id prov#year) cluster(city_id)
eststo R4: reghdfe gini rp if year>=2014 & year!=2022, absorb(city_id prov#year) cluster(prov)
esttab R1 R2 R3 R4 using "$OUT/robustness.csv", replace keep(rp) ///
    b(%7.4f) se(%7.4f) star(* 0.1 ** 0.05 *** 0.01) stats(N N_clust) csv
esttab R1 R2 R3 R4, keep(rp) b(%7.4f) se(%7.4f) star(* 0.1 ** 0.05 *** 0.01) stats(N N_clust)

*------------------------------------------------------------------------------
* 5. Higher education (Table 5)
*    (1) within-city OLS of the Gini on college density, 2009-2015
*    (2) robot exposure and college density: exposure x trend, 2009-2015
*    (3)-(4) main-window estimates split by base-period college density
*------------------------------------------------------------------------------
gen cg100 = coll_grad / 100            // college students per 10,000 residents, in hundreds
eststo clear
eststo E1: reghdfe gini cg100 if inrange(year,2009,2015), absorb(city_id year) cluster(city_id)
cap drop rp rp_eu rp_oth
gen rp     = rob_chn * (year - 2009)
gen rp_eu  = rob_eu  * (year - 2009)
gen rp_oth = rob_oth * (year - 2009)
eststo E2: ivreghdfe cg100 (rp = rp_eu rp_oth) if inrange(year,2009,2015), ///
    absorb(city_id prov#year) cluster(city_id)

drop rp rp_eu rp_oth
gen rp      = rob_chn * post_dmsp
gen rp_eu   = rob_eu  * post_dmsp
gen rp_oth  = rob_oth * post_dmsp
gen rph     = rp * highedu
gen rph_eu  = rp_eu * highedu
gen rph_oth = rp_oth * highedu
gen hpost   = highedu * post_dmsp
eststo E3: ivreghdfe gini hpost (rp rph = rp_eu rp_oth rph_eu rph_oth) if inrange(year,2000,2013), ///
    absorb(city_id prov#year) cluster(city_id)
lincom rp + rph
estadd scalar b_high  = r(estimate)
estadd scalar se_high = r(se)
estadd scalar p_high  = r(p)
eststo E4: ivreghdfe theil hpost (rp rph = rp_eu rp_oth rph_eu rph_oth) if inrange(year,2000,2013), ///
    absorb(city_id prov#year) cluster(city_id)
lincom rp + rph
estadd scalar b_high  = r(estimate)
estadd scalar se_high = r(se)
estadd scalar p_high  = r(p)
esttab E1 E2 E3 E4 using "$OUT/T5_education.csv", replace keep(cg100 rp rph) ///
    b(%7.4f) se(%7.4f) star(* 0.1 ** 0.05 *** 0.01) ///
    stats(N N_clust widstat b_high se_high p_high) csv
esttab E1 E2 E3 E4, keep(cg100 rp rph) b(%7.4f) se(%7.4f) star(* 0.1 ** 0.05 *** 0.01) ///
    stats(N N_clust widstat b_high se_high p_high)

* group sizes
preserve
keep city_id highedu edu0
duplicates drop
tab highedu, missing
sum edu0, detail
restore

log close
