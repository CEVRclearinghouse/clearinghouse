// This program is built to evaluate the long-term costs and benefits of AUD Treatment
// Anirban Basu & David Kim

set more off
capture mata:  mata drop AUDmodel()

version 13.0

mata:

real matrix AUDmodel(matrix basepop, /// Base Population
					 matrix remission, matrix futuretrt, matrix futureremission, matrix recover, matrix relapse, /// Treatment Probabilities
					 matrix aud_crime, matrix abst_crime, matrix aud_motor, matrix abst_motor, /// State-dependent Probabilities 
					 matrix aud_labor, matrix aud_labor_age, matrix abst_labor, matrix abst_labor_age, /// State-dependent Probabilities
					 matrix arrest, matrix conviction, matrix releaserate, matrix injail_futuretrt, /// Crime & Legal Probabilities  
					 matrix mortality, matrix mortality_d, matrix aud_mortality, matrix aud_mortality_d, matrix injail_mortality,  ///Mortality
					 matrix trtqol, matrix aud_qol, matrix abst_qol, matrix jail_qol, /// Quality of Life
					 matrix trtcosts, matrix timecosts, matrix futuretrtcosts, ///Treatment Costs
					 matrix aud_costs, matrix abst_costs, matrix eol_costs, /// Healthcare and End-of-Life care costs
					 matrix jail_costs, matrix justice_costs, matrix crime_costs, matrix crime_costs_qol, matrix motor_costs, /// Legal Costs 
					 matrix prodc, matrix prodc_age, matrix consmp, matrix consmp_age, /// Societal Costs
					 scalar horizon, scalar prnt) 		{	
					
/* Initializations */
real scalar  rho, ages, dm, dm1, trucatnum, beta
rho = 0.03			// discount rate

/*year */
costs = trtcosts :+ timecosts 
qalys = J(5,10000,0)
lys= J(5,10000,0)
aud_yr = J(5,10000,0)
emply_yr = J(5,10000,0)
jail_yr = J(5,10000,0)
dead = J(5,10000,0)
prod = J(5,10000,0)
consm = J(5,10000,0)
legal =J(5,10000,0)
crime_num=J(5,10000,0)
mva_num=J(5,10000,0)
curr_convicted = J(5,10000,0)
injail_aud = J(5,10000,0)
injail_abst = J(5,10000,0)
healthcosts = J(5,10000,0)

survival1 = J(horizon+1,10000,1)
survival2 = J(horizon+1,10000,1)
survival3 = J(horizon+1,10000,1)
survival4 = J(horizon+1,10000,1)
survival5 = J(horizon+1,10000,1)

//After initital AUD treatment
abst = remission
aud = 1:- remission

// Carry these pop through Markov model
for (i=1; i<= horizon; i++) {
	beta = (1/((1+rho)^(i-1)))
	/* release from injail & add to alive non-institutionalized population */
	aud = aud :+ injail_aud:*releaserate 
	abst = abst :+ injail_abst:*releaserate
	
	injail_aud = injail_aud:*(1:-releaserate)
	injail_abst = injail_abst:*(1:-releaserate)

	/* TRANSITION BETWEEN AUD & ABST */
	
	// in non-instituionalized population
	temp1 = aud:*futuretrt		// some get treatment
	temp2 = aud:*recover 		// some naturally recover
	temp3 = abst:*relapse		// some relapse

	costs = costs :+ (temp1:*futuretrtcosts):*beta

	healthcosts = healthcosts :+ (temp1:*futuretrtcosts):*beta
	aud = aud :+ temp3 :- (temp1:*futureremission) :- temp2
	abst = abst :+ (temp1:*futureremission) :+ temp2 :- temp3
	
	// within jail
	temp1 =  injail_aud:*injail_futuretrt	//some get treatment in future
	temp2 = injail_aud:*recover
	costs = costs :+ temp1:*futuretrtcosts:*beta
	
	healthcosts = healthcosts :+ (temp1:*futuretrtcosts):*beta
	injail_aud = injail_aud :- temp1:*futureremission :- temp2
	injail_abst = injail_abst :+ temp1:*futureremission :+ temp2

	/* QALYs & COSTS among non-institutionalized who die */
	temp = aud:*aud_mortality_d[i,.] :+ abst:*mortality_d[i,.]
	dead = dead :+ temp
	lys = lys :+ (temp:*0.5):*beta
	qalys = qalys :+ ((aud:*aud_mortality_d[i,.]:*aud_qol :+ abst:*mortality_d[i,.]:*abst_qol):*0.5):*beta
	costs  = costs :+ (eol_costs:*(temp)):*beta
	
	healthcosts = healthcosts :+ ((aud:*aud_mortality_d[i,.]:*aud_costs[i,.] :+ abst:*mortality_d[i,.]:*abst_costs[i,.]):*0.5 :+ eol_costs:*(temp)):*beta
	
	prod = prod :+ ((aud:*aud_mortality_d[i,.]:*aud_labor_age[i,.]:*prodc):+(abst:*mortality_d[i,.]:*abst_labor_age[i,.]:*prodc)):*0.5:*beta
	consm = consm :+ ((aud:*aud_mortality_d[i,.]:*consmp_age[i,.]):+(abst:*mortality_d[i,.]:*consmp_age[i,.])):*0.5:*beta

	costs = costs :- (aud:*aud_mortality_d[i,.]:*aud_labor_age[i,.]:*prod :+ abst:*mortality_d[i,.]:*abst_labor_age[i,.]:*prod):*0.5:*beta  ///
			:+ (aud:*aud_mortality_d[i,.]:*consmp_age[i,.] :+ abst:*mortality_d[i,.]:*consmp_age[i,.]):*0.5:*beta

	aud = 	aud :* (1 :- aud_mortality_d[i,.])
	abst = abst :* (1 :- mortality_d[i,.])

	// Health care & Productivity & consumption
	costs  = costs :+ ((aud:*aud_costs[i,.] :+ abst:*abst_costs[i,.])):*beta
	
	healthcosts = healthcosts :+ ((aud:*aud_costs[i,.] :+ abst:*abst_costs[i,.])):*beta
	
	prod = prod :+ (aud:*aud_labor_age[i,.]:*prodc :+ abst:*abst_labor_age[i,.]:*prodc):*beta
	consm = consm :+ (aud:*consmp_age[i,.] :+ abst:*consmp_age[i,.]):*beta
	costs = costs :- (aud:*aud_labor_age[i,.]:*prod :+ abst:*abst_labor_age[i,.]:*prod):*beta  :+ (aud:*consmp_age[i,.] :+ abst:*consmp_age[i,.]):*beta
	lys = lys :+ (aud :+ abst):*beta
	qalys = qalys :+ (aud:*aud_qol :+ abst:*abst_qol):*beta

	/* QALYs & COSTS among inmates who die*/ 
	dead = dead :+ (injail_aud:*(aud_mortality_d[i,.] :+ injail_mortality) :+ injail_abst:*(mortality_d[i,.] :+ injail_mortality))
	lys = lys :+ ((injail_aud:*(aud_mortality_d[i,.] :+ injail_mortality) :+ injail_abst:*(mortality_d[i,.] :+ injail_mortality)):*0.5):*beta
	qalys = qalys :+ ((injail_aud:*(aud_mortality_d[i,.] :+ injail_mortality):*(aud_qol:+jail_qol) :+ injail_abst:*(mortality_d[i,.] :+ injail_mortality):*(abst_qol:+jail_qol)):*0.5):*beta	

	healthcosts = healthcosts :+ (eol_costs:*(injail_aud:*(aud_mortality_d[i,.] :+ injail_mortality) :+ injail_abst:*(mortality_d[i,.] :+ injail_mortality))):*beta   ///

	legal = legal :+ (jail_costs:*(injail_aud:*(aud_mortality_d[i,.] :+ injail_mortality) :+ injail_abst:*(mortality_d[i,.] :+ injail_mortality))*0.5):*beta

	costs  = costs :+ (eol_costs:*(injail_aud:*(aud_mortality_d[i,.] :+ injail_mortality) :+ injail_abst:*(mortality_d[i,.] :+ injail_mortality))):*beta   ///
						:+ (jail_costs:*(injail_aud:*(aud_mortality_d[i,.] :+ injail_mortality) :+ injail_abst:*(mortality_d[i,.] :+ injail_mortality))*0.5):*beta
	
	injail_aud = injail_aud:*(1:-(aud_mortality_d[i,.] :+ injail_mortality))
	injail_abst = injail_abst:*(1:-(mortality_d[i,.] :+ injail_mortality))

	// QALYs & COSTS among those who remain incarcerated and do not die
	lys = lys :+ (injail_aud :+ injail_abst):*beta
	qalys = qalys :+ (injail_aud:*(aud_qol:+jail_qol) :+ injail_abst:*(abst_qol:+jail_qol)):*beta
	legal = legal :+ (jail_costs:*(injail_aud :+ injail_abst)):*beta
	costs = costs :+ ((injail_aud:*aud_costs[i,.] :+ injail_abst:*abst_costs[i,.]) :+ jail_costs:* (injail_aud :+ injail_abst)):*beta
	
	healthcosts = healthcosts :+ (injail_aud:*aud_costs[i,.] :+ injail_abst:*abst_costs[i,.]):*beta

	/* QALYs & COSTS & BEHAVIOR among non-instituionalized in current year who do not die */
	/* Commit crime & get arrested & incarcerated or not*/ 

	legal = legal :+ (aud:*(aud_crime:*(crime_costs:+crime_costs_qol):+ aud_motor:*motor_costs :+ aud_crime:*arrest:*justice_costs :+ aud_crime:*arrest:*conviction:*jail_costs:*0.5) ///
					  :+ abst:*(abst_crime:*(crime_costs:+crime_costs_qol):+ abst_motor:*motor_costs :+ abst_crime:*arrest:*justice_costs :+ abst_crime:*arrest:*conviction:*jail_costs:*0.5)):*beta
	costs = costs :+ (aud:*(aud_crime:*(crime_costs:+crime_costs_qol) :+ aud_motor:*motor_costs :+ aud_crime:*arrest:*justice_costs :+ aud_crime:*arrest:*conviction:*jail_costs:*0.5) ///
					  :+ abst:*(abst_crime:*(crime_costs:+crime_costs_qol) :+ abst_motor:*motor_costs :+ abst_crime:*arrest:*justice_costs :+ abst_crime:*arrest:*conviction:*jail_costs:*0.5)):*beta

	injail_aud = injail_aud :+ aud:*(aud_crime:*arrest:*conviction)
	injail_abst = injail_abst :+ abst:*(abst_crime:*arrest:*conviction)
	
	aud = aud:*(1 :- (aud_crime:*arrest:*conviction))
	abst = abst:*(1 :- (abst_crime:*arrest:*conviction))
	
	survival1[i+1,.] = 1:- dead[1,.]
	survival2[i+1,.] = 1:- dead[2,.]
	survival3[i+1,.] = 1:- dead[3,.]
	survival4[i+1,.] = 1:- dead[4,.]
	survival5[i+1,.] = 1:- dead[5,.]

/*Impact Inventory Estimate*/
	aud_yr = aud_yr:+aud // Years in AUD State
	emply_yr = emply_yr :+ ((aud:*aud_labor_age[i,.]):+(abst:*abst_labor_age[i,.])) // Years Employed
	jail_yr = jail_yr :+ (injail_aud :+ injail_abst) // Years in Jail
	crime_num = crime_num :+ ((aud:*aud_crime):+(abst:*abst_crime)) // # of Crimes
	mva_num = mva_num :+ ((aud:*aud_motor):+(abst:*abst_motor)) // # of Motor Vehicle Accident 
	}
	
healthcosts2 = costs :- (trtcosts:+ timecosts :+ legal :- prod :+ consm)

sumcost1 = trtcosts :+ healthcosts
sumcost2 = sumcost1 :+ timecosts
sumcost3 = sumcost2 :- prod :+ consm
sumcost4 = sumcost3 :+ legal


/*To Generate PSA Results*/
//Data
printf("\n")
printf("Overall Costs - Healthcare Sector Perspective \n")
sumcost1_tx = transposeonly(sumcost1)
sumcost1_tx[1..10,.]
printf("\n")
printf("Overall Costs - Societal Perspective \n")
sumcost4_tx = transposeonly(sumcost4)
sumcost4_tx[1..10,.]

printf("\n")
printf("Overall QALYs Gained \n")
qalys_tx = transposeonly(qalys)
qalys_tx[1..10,.]

//Net Monetary Benefit@$50,000/QALY 
nmb1_50K = (qalys:*50000) :- sumcost1
nmb4_50K = (qalys:*50000) :- sumcost4

//Net Monetary Benefit@$100,000/QALY 
nmb1_100K = (qalys:*100000) :- sumcost1
nmb4_100K = (qalys:*100000) :- sumcost4

if (prnt==1) {
printf("\n")
printf("Impact Inventory: Life Years, QALYS, Years in AUD, Years Employed \n")
for (i=1; i<= 5; i++) {
printf("%f", round(mean(lys[i,.]'),0.01)), printf(" ("), printf("%f", round(sqrt(variance(lys[i,.]')),0.01)), printf("),  "), ///
printf("%f", round(mean(qalys[i,.]'),0.01)), printf(" ("), printf("%f", round(sqrt(variance(qalys[i,.]')),0.01)), printf("),  "), ///
printf("%f", round(mean(aud_yr[i,.]'),0.01)), printf(" ("), printf("%f", round(sqrt(variance(aud_yr[i,.]')),0.01)), printf("),  "), ///
printf("%f", round(mean(emply_yr[i,.]'),0.01)), printf(" ("), printf("%f", round(sqrt(variance(emply_yr[i,.]')),0.01)), printf(")  \n")
	}

printf("\n")
printf("Impact Inventory: Years in Jail, # of Crimes, # of MVA \n")
for (i=1; i<= 5; i++) {
printf("%f", round(mean(jail_yr[i,.]'),0.01)), printf(" ("), printf("%f", round(sqrt(variance(jail_yr[i,.]')),0.01)), printf("),  "), /// 
printf("%f", round(mean(crime_num[i,.]'),0.01)), printf(" ("), printf("%f", round(sqrt(variance(crime_num[i,.]')),0.01)), printf("),  "), /// 
printf("%f", round(mean(mva_num[i,.]'),0.01)), printf(" ("), printf("%f", round(sqrt(variance(mva_num[i,.]')),0.01)), printf(") \n") 
	}

printf("\n")
printf("Impact Inventory: Treatment Costs, Time Costs, Healthcare Costs \n")
for (i=1; i<= 5; i++) {
printf("%f", round(mean(trtcosts[i,.]'),1)), printf(" ("), printf("%f", round(sqrt(variance(trtcosts[i,.]')),1)), printf("),  "), /// 
printf("%f", round(mean(timecosts[i,.]'),1)), printf(" ("), printf("%f", round(sqrt(variance(timecosts[i,.]')),1)), printf("),  "), /// 
printf("%f", round(mean(healthcosts[i,.]'),1)), printf(" ("), printf("%f", round(sqrt(variance(healthcosts[i,.]')),1)), printf(") \n") 
	}
	
printf("\n")
printf("Impact Inventory: Productivity Costs, Consumption Costs, Legal Costs \n")
for (i=1; i<= 5; i++) {
printf("%f", round(mean(prod[i,.]'),1)), printf(" ("), printf("%f", round(sqrt(variance(prod[i,.]')),1)), printf("),  "), /// 
printf("%f", round(mean(consm[i,.]'),1)), printf(" ("), printf("%f", round(sqrt(variance(consm[i,.]')),1)), printf("),  "), /// 
printf("%f", round(mean(legal[i,.]'),1)), printf(" ("), printf("%f", round(sqrt(variance(legal[i,.]')),1)), printf(")  \n") 
	}
	
printf("\n")
printf("Total Costs: Healthcare Sector and Societal \n")
for (i=1; i<= 5; i++) {
printf("%f", round(mean(sumcost1[i,.]'),1)), printf(" ("), printf("%f", round(sqrt(variance(sumcost1[i,.]')),1)), printf("),  "), ///
printf("%f", round(mean(sumcost4[i,.]'),1)), printf(" ("), printf("%f", round(sqrt(variance(sumcost4[i,.]')),1)), printf(")  \n")
	}

printf("\n")
printf("NMB@100K - Healthcare Sector and Societal \n")
for (i=1; i<= 5; i++) {
printf("%f", round(mean(nmb1_100K[i,.]'),1)), printf(" ("), printf("%f", round(sqrt(variance(nmb1_100K[i,.]')),1)), printf("),  "), ///
printf("%f", round(mean(nmb4_100K[i,.]'),1)), printf(" ("), printf("%f", round(sqrt(variance(nmb4_100K[i,.]')),1)), printf(") \n")
	}
	
printf("\n")
printf("Survival Rates over the Time Horizon\n")
mean(survival1')' , mean(survival2')' ,mean(survival3')' ,mean(survival4')' ,mean(survival5')'  
}
maxnmb = colmax(nmb1_50K)
voi1 = maxnmb :- nmb1_50K[5,.] 

maxnmb = colmax(nmb4_50K)
voi4 = maxnmb :- nmb4_50K[5,.] 

printf("\n")
printf("Per Subject Value of Infomation Estimates\n")
printf("Healthcare Sector and Societal\n")
if (prnt==1) {
printf("%f ", round(mean(voi1'),1)), printf("%f ", round(mean(voi4'),1))
}
voi = mean(voi1'), mean(voi4')

return(voi)
}

mata mosave AUDmodel(), dir(C:\Users\ddkim62\Desktop) replace
end 

******************************
**********INPUT DATA**********
******************************

global devsize = 10000
global horizon = 56	/* ages 45 to 100 */
mata: basepop = J(1,1,1) 	
///  annual populaton size by age  (46X1)
mata: horizon = J(1,1, $horizon)

/***************** BETA DISTRIBUTIONS
alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
beta <- alpha * (1 / mu - 1) 
************************************/

/*Treatment Probabilities*/
mata: rseed(123456789)
mata: remission = 	J(1,$devsize, 1):*rbeta(1,$devsize, 9.58, 61.9) \ /// Remission rates for CBI
					J(1,$devsize, 1):*rbeta(1,$devsize, 22.9, 123) \ /// Remission rates for MM
					J(1,$devsize, 1):*rbeta(1,$devsize, 22.3, 81.8) \ /// Remission rates for MM+Nalt
					J(1,$devsize, 1):*rbeta(1,$devsize, 18.3, 77.3) \ /// Remission rates for MM+ Acam
					J(1,$devsize, 1):*rbeta(1,$devsize, 15.8, 73.9)   //  Remission rates for MM+Nalt+Acap

mata: rseed(123456789)
mata: futuretrt = J(5,$devsize, 1):*rbeta(1,$devsize, 122.7, 2275) 

mata: rseed(123456789)
mata: futureremission  = J(5,$devsize, 1):*rbeta(1,$devsize, 31.4, 68.6) 	 	

mata: rseed(123456789)
mata: recover = J(5,$devsize, 1):*rbeta(1,$devsize, 66.0, 574.6)
 
mata: rseed(123456789)
mata: relapse =J(5,$devsize, 1):*rbeta(1,$devsize, 30, 970) 		

/*State-dependent probabilites*/

mata: rseed(123456789)							
mata: aud_crime = J(5,$devsize, 1):*rbeta(1,$devsize, 1.767, 8.233)	

mata: rseed(123456789)							
mata: abst_crime = J(5,$devsize, 1):*rbeta(1,$devsize, 3.246, 96.754)	

mata: rseed(123456789)							
mata: aud_motor = J(5,$devsize, 1):*rbeta(1,$devsize, 5.2, 94.8)	

mata: rseed(123456789)							
mata: abst_motor = J(5,$devsize, 1):*rbeta(1,$devsize, 12.6, 987.4)	

*Labor Force Particilation Rates
mata: rseed(123456789)							
mata: aud_labor = J(5,$devsize, 1):*rbeta(1,$devsize, 66.3, 33.7)	

*Age-specific LFP among AUD individuals
mata: aud_labor_age = 	rbeta(15,$devsize,78.6,21.4) \ /// age 45-59
						rbeta(5,$devsize,50.0,50.0) \ /// age 60-64
						rbeta(5,$devsize,27.4,72.6) \ /// age 65-69
						rbeta(5,$devsize,16.7,83.3) \ /// age 70-74
						rbeta(25,$devsize,6.5,93.5) \ /// age 75+
						J(1,$devsize,0) /*age > 100*/
mata: aud_labor_age = aud_labor_age[1..$horizon,.]
	
*LFP Among LRD individuals 
mata: rseed(123456789)							
mata: abst_labor = J(5,$devsize, 1):*rbeta(1,$devsize, 0.719, 0.281)

*Age-specific LFP among LRD individuals
mata: abst_labor_age = 	rbeta(15,$devsize,68.4,31.6) \ /// age 45-59
						rbeta(5,$devsize,57.1,42.9) \ /// age 60-64
						rbeta(5,$devsize,27.4,72.6) \ /// age 65-69
						rbeta(5,$devsize,16.7,83.3) \ /// age 70-74
						rbeta(25,$devsize,6.5,93.5) \ /// age 75+
						J(1,$devsize,0) /*age > 100*/
mata: abst_labor_age = abst_labor_age[1..$horizon,.]
						
/*Crime & Legal Probabilities*/
mata: rseed(123456789)							
mata: arrest = J(5,$devsize, 1):*rbeta(1,$devsize, 2.13, 7.87)	

mata: rseed(123456789)							
mata: conviction = J(5,$devsize, 1):*rbeta(1,$devsize, 3.38, 6.62)	

mata: rseed(123456789)
mata: releaserate = J(5,$devsize, 1):*rbeta(1,$devsize, 4.3, 5.7)

mata: rseed(123456789)
mata: injail_futuretrt = J(5,$devsize, 0)

/*Mortality*/
*Age-specific general mortality rates - deterministic*
mata: mortality_d = J(5,$devsize,0.0033) \ /// age 46-50
					J(5,$devsize,0.0051) \ /// age 51-55
					J(5,$devsize,0.0073) \ /// age 56-60
					J(5,$devsize,0.0105) \ /// age 61-65
					J(5,$devsize,0.0158) \ /// age 66-70
					J(5,$devsize,0.0240) \ /// age 71-75
					J(5,$devsize,0.0379) \ /// age 76-80
					J(5,$devsize,0.0618) \ /// age 81-85
					J(5,$devsize,0.1050) \ /// age 86-90
					J(5,$devsize,0.1709) \ /// age 91-95
					J(5,$devsize,0.2589) \ /// age 96-100
					J(1,$devsize,1.0000) /*age > 100*/
mata: mortality_d = mortality_d[1..$horizon,.]			
					
mata: rseed(123456789)					
mata: mortality = 	rbeta(5,$devsize,33,9967) \ /// age 46-50
					rbeta(5,$devsize,51,9949) \ /// age 51-55
					rbeta(5,$devsize,73,9927) \ /// age 56-60
					rbeta(5,$devsize,105,9895) \ /// age 61-65
					rbeta(5,$devsize,158,9842) \ /// age 66-70
					rbeta(5,$devsize,240,9760) \ /// age 71-75
					rbeta(5,$devsize,379,9621) \ /// age 76-80
					rbeta(5,$devsize,618,9382) \ /// age 81-85
					rbeta(5,$devsize,1050,8950) \ /// age 86-90
					rbeta(5,$devsize,1709,8291) \ /// age 91-95
					rbeta(5,$devsize,2589,7411) \ /// age 96-100
					J(1,$devsize,1.0000) /*age > 100*/
mata: mortality = mortality[1..$horizon,.]

*Age-specific relative-risk of mortality among AUD patients
mata: rseed(123456789)
mata: aud_mortality_RR = 	rnormal(5,$devsize,4.70,0.62) \ /// Age 46-50
							rnormal(10,$devsize,3.21,0.36) \ /// Age 51-60
							rnormal(40,$devsize,1.74,0.22) \ /// Age 60+
							J(1,$devsize,1.0000) /*Age 100+*/
mata: aud_mortality_RR = aud_mortality_RR[1..$horizon,.]

mata: aud_mortality_d = mortality_d:*aud_mortality_RR

mata: aud_mortality = mortality:*aud_mortality_RR

*Additive moretality rate in Jail - assuming to be zero
mata: rseed(123456789)
mata: injail_mortality = J(5,$devsize, 0) // additive mortality rate	

/*Quality of Life*/
mata: rseed(123456789)							
mata: aud_qol = J(5,$devsize, 1):*rbeta(1,$devsize, 37.7, 11.9)	

mata: rseed(123456789)							
mata: abst_qol = J(5,$devsize, 1):*rbeta(1,$devsize, 310.8, 72.9)

mata: rseed(123456789)
mata: jail_qol = J(5,$devsize, -1):*rbeta(1,$devsize, 325, 4010)

//# of Total Patients
mata: n_event = (157\309\309\303\305)

/*NAUSEA*/
mata: rseed(123456789)
//Nausea-related HRQoL Decrements 
mata: nausea_qol = J(5,10000,-1):*rbeta(1,10000,3.9,191.2)
mata: rseed(123456789)
mata: nausea_events = J(1,10000,0) \ /// # of Events for CBI Only
			 rbinomial(1,10000,309,0.21) \ /// For MM Only with Placebo
			 rbinomial(1,10000,309,0.34) \ /// For MM + Naltrexone
			 rbinomial(1,10000,303,0.24) \ /// For MM + Acamprosate
			 rbinomial(1,10000,305,0.42) //For MM + Nalt + Acam
			 
//Estimate the Expected Probability of Having Nausea 
mata: nausea_p = nausea_events :/ n_event
//Estimate the Expected HrQOL decrements due to Nausea
mata: mean_nausea_qol = nausea_p:*nausea_qol

/*VOMITING*/
//Vomiting-related HRQoL Decrements 
mata: rseed(123456789)
mata: vomit_qol = J(5,10000,-1):*rbeta(1,10000,3.9,191.2)
//Binomial Distributions Based on # of Events
mata: rseed(123456789)
mata: vomit_events = J(1,10000,0) \ /// # of Events for CBI Only
			 rbinomial(1,10000,309,0.09) \ /// For MM Only with Placebo
			 rbinomial(1,10000,309,0.15) \ /// For MM + Naltrexone
			 rbinomial(1,10000,303,0.09) \ /// For MM + Acamprosate
			 rbinomial(1,10000,305,0.18) //For MM + Nalt + Acam
			 
//Estimate the Expected Probability of Having Nausea 
mata: vomit_p = vomit_events :/ n_event
//Estimate the Expected HrQOL decrements due to Nausea
mata: mean_vomit_qol = vomit_p:*vomit_qol

/*DIARRHEA*/
//Vomiting-related HRQoL Decrements 
mata: rseed(123456789)
mata: diarrhea_qol = J(5,10000,-1):*rbeta(1,10000,3.9,191.2)
//Binomial Distributions Based on # of Events
mata: rseed(123456789)
mata: diarrhea_events = J(1,10000,0) \ /// # of Events for CBI Only
			 rbinomial(1,10000,309,0.35) \ /// For MM Only with Placebo
			 rbinomial(1,10000,309,0.31) \ /// For MM + Naltrexone
			 rbinomial(1,10000,303,0.65) \ /// For MM + Acamprosate
			 rbinomial(1,10000,305,0.56) //For MM + Nalt + Acam
			 
//Estimate the Expected Probability of Having Nausea 
mata: diarrhea_p = diarrhea_events :/ n_event
//Estimate the Expected HrQOL decrements due to Nausea
mata: mean_diarrhea_qol = diarrhea_p:*diarrhea_qol

/*Loss of Appetite*/
//Vomiting-related HRQoL Decrements 
mata: rseed(123456789)
mata: lossappetitie_qol = J(5,10000,-1):*rbeta(1,10000,15.0,235)
//Binomial Distributions Based on # of Events
mata: rseed(123456789)
mata: lossappetitie_events = J(1,10000,0) \ /// # of Events for CBI Only
			 rbinomial(1,10000,309,0.13) \ /// For MM Only with Placebo
			 rbinomial(1,10000,309,0.21) \ /// For MM + Naltrexone
			 rbinomial(1,10000,303,0.19) \ /// For MM + Acamprosate
			 rbinomial(1,10000,305,0.25) //For MM + Nalt + Acam
			 
//Estimate the Expected Probability of Having Nausea 
mata: lossappetitie_p = lossappetitie_events :/ n_event
//Estimate the Expected HrQOL decrements due to Nausea
mata: mean_lossappetitie_qol = lossappetitie_p:*lossappetitie_qol

//Overal Disutility of Adverse Events during Treatment
mata: trt_disutility = mean_nausea_qol :+ mean_vomit_qol :+ ///
		mean_diarrhea_qol :+ mean_lossappetitie_qol

mata: trtqol = aud_qol :+ trt_disutility	

/*Treatment Costs*/	

*Costs of Treatment
mata: rseed(123456789)					
mata: tr = 	rnormal(1,$devsize, 640, 18.1) \ ///
			rnormal(1,$devsize, 474, 7.52) \ ///
			rnormal(1,$devsize, 777, 19.5) \ ///
			rnormal(1,$devsize, 866, 22.4) \ ///
			rnormal(1,$devsize, 1162, 36.7) 
			
mata: trtcosts = J(5,$devsize, 1):*(tr:*(tr:>0))			

*Time costs related to treatment
mata: rseed(123456789)					
mata: tr = 	rnormal(1,$devsize, 1164, 102) \ ///
			rnormal(1,$devsize, 1014, 86.9) \ ///
			rnormal(1,$devsize, 1090, 111) \ ///
			rnormal(1,$devsize, 1128, 143) \ ///
			rnormal(1,$devsize, 1046, 101)
			
mata: timecosts = J(5,$devsize, 1):*(tr:*(tr:>0))			

*Future Treatment costs (average of all tree 
mata: rseed(123456789)	
mata: tr = rnormal(1,$devsize, 1872, 252)			
mata: futuretrtcosts = J(5,$devsize, 1):*(tr:*(tr:>0))	

/*Future Healthcare costs*/
*Among General population - based on MEPS 2003-2012
mata: rseed(123456789)					
mata: tr = 	rnormal(10,$devsize, 5039, 91) \ /// age 46-55
			rnormal(10,$devsize, 7780, 188) \ /// age 56-65
			rnormal(10,$devsize, 9814, 183) \ /// age 66-75
			rnormal(25,$devsize, 12037, 191) \ /// age 75+
			rnormal(1,$devsize, 12037, 191) /*Age 100+*/
			
mata: tr = tr[1..$horizon,.]
mata: abst_costs = J($horizon,$devsize, 1):*(tr:*(tr:>0))

*Among AUD patients - based on MEPS 2003-2012
mata: rseed(123456789)	
mata: tr = 	rnormal(20,$devsize, 12230, 1420) \ /// age 46-65
			rnormal(35,$devsize, 31196, 8944) \ /// age 65+
			rnormal(1,$devsize, 31196, 8944) /*Age 100+*/

mata: tr = tr[1..$horizon,.]		
mata: aud_costs = J($horizon,$devsize, 1):*(tr:*(tr:>0))						

*End-of-Life Care Expenditures
mata: rseed(123456789)					
mata: eol_costs = rnormal(1,$devsize, 37421, 4857)		
mata: eol_costs = J(5,$devsize, 1):*eol_costs	

/*Legal Costs*/
mata: rseed(123456789)	
mata: jail_costs = rnormal(1,$devsize, 34710, 9235)		
mata: jail_costs = J(5,$devsize, 1):*jail_costs

mata: rseed(123456789)					
mata: justice_costs = rnormal(1,$devsize, 1445, 270)		
mata: justice_costs = J(5,$devsize, 1):*justice_costs	

mata: rseed(123456789)					
mata: crime_costs = rnormal(1,$devsize, 1568, 314)		
mata: crime_costs = J(5,$devsize, 1):*crime_costs	

mata: rseed(123456789)					
mata: crime_costs_qol = rnormal(1,$devsize, 14245, 2849)		
mata: crime_costs_qol = J(5,$devsize, 1):*crime_costs_qol

mata: rseed(123456789)					
mata: motor_costs = rnormal(1,$devsize, 8581, 1716)		
mata: motor_costs = J(5,$devsize, 1):*motor_costs	

/*Societal Costs*/
*Average Annual Earnings
mata: rseed(123456789)	
*Adjusting for Fringe Benefits
mata: benefit_rate = 0.462
						
mata: prodc = J(5,$devsize, 1):*rnormal(1, $devsize, 50017, 500)
mata: prodc = (1+benefit_rate)*prodc
*Age-specific Average Annual Earnings
mata: rseed(123456789)
mata: p_c = rnormal(10,$devsize, 52749, 721) \ /// age 46-55
			rnormal(10,$devsize, 53329, 825) \ /// age 56-65
			rnormal(10,$devsize, 41317, 1568) \ /// age 66-75
			rnormal(25,$devsize, 37300, 3432) \ /// age 75+
			rnormal(1,$devsize, 37300, 3432) /*Age 100+*/
mata: p_c = p_c[1..$horizon,.]
			
mata: prodc_age = J($horizon,$devsize,1):*(p_c:*(p_c:>0))
mata: prodc_age = (1+benefit_rate):*prodc_age
mata: prodc_age[.,1..7]
	

*Average Consumption Costs (Non-healthcare Expenditure)
mata: rseed(123456789)							
mata: consmp = J(5,$devsize, 1):*rnormal(1, $devsize, 18988, 330)

*Age-specific Average Annual Earnings
mata: rseed(123456789)
mata: c_c = rnormal(10,$devsize, 21009, 683) \ /// age 46-55
			rnormal(10,$devsize, 24530, 707) \ /// age 56-65
			rnormal(10,$devsize, 21878, 901) \ /// age 66-75
			rnormal(25,$devsize, 18420, 914) \ /// age 75+
			rnormal(1,$devsize, 18420, 914) /*Age 100+*/
mata: c_c = c_c[1..$horizon,.]
mata: consmp_age = J($horizon,$devsize,1):*(c_c:*(c_c:>0))	
mata: consmp_age[.,1..7]

*************************************
**********END OF DATA INPUT**********
*************************************

mata: voi = AUDmodel(	basepop, 	///  annual populaton size by age  (1X1)
						/// Treatment Probabilities
						remission, 	/// remission probability by trt (4X10000)
						futuretrt, 	///  pr of seeking future trt (4X10000)
						futureremission, /// Avg. pr(remission from future trt) (4X10000)
						recover, 	///  pr(recovery from AUD) (4X10000)
						relapse, 	///  pr(relapse from abst) (4X10000)
						///State-dependent probabilities
						aud_crime,		/// Pr(crime |AUD) (4X10000)
						abst_crime,		/// Pr(crime |ABST) (4X10000)
						aud_motor,		/// Pr(motor acc |AUD) (4X10000)
						abst_motor,		/// Pr(motor acc |ABST) (4X10000)
						aud_labor,		/// Pr(labor market |AUD) (4X10000)
						aud_labor_age, 	///	Pr(labor market |AUD) by Age (horizonX10000)
						abst_labor,		/// Pr(labor market |ABST) (4X10000)
						abst_labor_age, /// Pr(labor market |ABST) by Age (horizonX10000)
						///Crime&Legal Probabilities
						arrest,				/// Pr(arest|crime) (4X10000)
						conviction,			/// Pr(conviction|arrest) (4x10000) 
						releaserate, 		///  release rate from jail (4X10000)
						injail_futuretrt, 	///  pr of seeking future trt in jail (4X10000)
						///Mortality
						mortality,			/// Pr(AUD mortality) by age (horizonX10000) 
						mortality_d, 		///	Pr(AUD mortality) by age (horizonX10000) - Deterministic
						aud_mortality, 		/// Pr(AUD mortality) by age (horizonX10000)
						aud_mortality_d, 	/// Pr(AUD mortality) by age (horizonX10000)
						injail_mortality, 	/// Added Pr(mortality[i,.] in jail) (4X10000)
						///Quality of Life
						trtqol,			/// QOL for trtstate (4X10000)
						aud_qol,		/// QOL for aud state (4X10000)
						abst_qol,		/// QOL for abst state (4X10000)
						jail_qol,       /// QOL decrements for being incarcerated (4X10000)
						///Treatment costs
						trtcosts,		/// Avg. costs of current trt (4X10000)
						timecosts,		/// Avg. costs of time seeking trt (4X10000)
						futuretrtcosts, /// Avg. costs of treat in future (4X10000)
						///Future Healthcare costs
						aud_costs,		/// Costs for aud state (horizonX10000)
						abst_costs,		/// Costs for abst state (horizonX10000)
						///End-of-Life Care Expenditures
						eol_costs,		/// Costs for eol state (4X10000)
						///Legal Costs
						jail_costs,		/// Costs for jail state (4X10000)
						justice_costs,	/// Costs for arrest (4X10000)
						crime_costs,	/// Costs for crime (tangilble costs) (4X10000)
						crime_costs_qol, /// Costs for crime (intangilble costs) (4X10000)
						motor_costs,	/// Costs for motor vehicle accidents (4X10000)		
						///Societal Costs
						prodc,		/// Avg. productivity (4X10000)
						prodc_age, 	/// Age-specific productivity (horizon*10000)
						consmp,		/// Avg. Consumption  (4X10000)
						consmp_age,	/// Age-specific consumption (horizon*10000)
						///
						horizon, ///
						1)
