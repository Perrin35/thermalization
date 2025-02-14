OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7186681) q[0];
sx q[0];
rz(-1.0942425) q[0];
sx q[0];
rz(-2.8835468) q[0];
rz(2.1482422) q[1];
sx q[1];
rz(1.8408096) q[1];
sx q[1];
rz(9.169133) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48707552) q[0];
sx q[0];
rz(-2.1571113) q[0];
sx q[0];
rz(2.4357027) q[0];
x q[1];
rz(-1.6338111) q[2];
sx q[2];
rz(-1.9587226) q[2];
sx q[2];
rz(0.20252075) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.62754831) q[1];
sx q[1];
rz(-1.6395634) q[1];
sx q[1];
rz(-1.7538944) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3575451) q[3];
sx q[3];
rz(-2.0590326) q[3];
sx q[3];
rz(2.8247716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51097441) q[2];
sx q[2];
rz(-1.7642085) q[2];
sx q[2];
rz(2.0261436) q[2];
rz(-0.014178064) q[3];
sx q[3];
rz(-1.7970128) q[3];
sx q[3];
rz(0.43629638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2071335) q[0];
sx q[0];
rz(-2.5196228) q[0];
sx q[0];
rz(1.2472664) q[0];
rz(2.6994052) q[1];
sx q[1];
rz(-1.7555883) q[1];
sx q[1];
rz(2.3242548) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8431664) q[0];
sx q[0];
rz(-2.0930556) q[0];
sx q[0];
rz(2.9199615) q[0];
rz(-pi) q[1];
rz(-1.1306612) q[2];
sx q[2];
rz(-1.8225742) q[2];
sx q[2];
rz(1.2680858) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4358959) q[1];
sx q[1];
rz(-1.4053939) q[1];
sx q[1];
rz(-2.4981605) q[1];
rz(-pi) q[2];
rz(-2.4659791) q[3];
sx q[3];
rz(-0.88425501) q[3];
sx q[3];
rz(-0.10304777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.77357972) q[2];
sx q[2];
rz(-2.7648338) q[2];
sx q[2];
rz(-0.22029857) q[2];
rz(-1.1822654) q[3];
sx q[3];
rz(-1.1743098) q[3];
sx q[3];
rz(-3.0784472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050506266) q[0];
sx q[0];
rz(-1.3171221) q[0];
sx q[0];
rz(-2.0029946) q[0];
rz(-2.7298722) q[1];
sx q[1];
rz(-0.76985923) q[1];
sx q[1];
rz(2.128111) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0419755) q[0];
sx q[0];
rz(-1.8238153) q[0];
sx q[0];
rz(-1.5861041) q[0];
rz(-1.8455681) q[2];
sx q[2];
rz(-0.45598511) q[2];
sx q[2];
rz(-2.6528461) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.24559969) q[1];
sx q[1];
rz(-0.87973824) q[1];
sx q[1];
rz(-0.57180239) q[1];
rz(-pi) q[2];
rz(-0.8441505) q[3];
sx q[3];
rz(-1.4733847) q[3];
sx q[3];
rz(1.0507492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.98075214) q[2];
sx q[2];
rz(-1.1557121) q[2];
sx q[2];
rz(-1.2551003) q[2];
rz(0.27215019) q[3];
sx q[3];
rz(-2.3641219) q[3];
sx q[3];
rz(3.0009771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.960152) q[0];
sx q[0];
rz(-0.93638268) q[0];
sx q[0];
rz(-0.61035672) q[0];
rz(2.4329674) q[1];
sx q[1];
rz(-1.0162063) q[1];
sx q[1];
rz(-1.5707387) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8886531) q[0];
sx q[0];
rz(-1.7084689) q[0];
sx q[0];
rz(2.987079) q[0];
rz(-2.8547457) q[2];
sx q[2];
rz(-1.45668) q[2];
sx q[2];
rz(0.6335887) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.849087) q[1];
sx q[1];
rz(-1.6189112) q[1];
sx q[1];
rz(-0.27166697) q[1];
x q[2];
rz(2.3409178) q[3];
sx q[3];
rz(-1.281637) q[3];
sx q[3];
rz(1.9137933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2107971) q[2];
sx q[2];
rz(-2.1288629) q[2];
sx q[2];
rz(-2.5861758) q[2];
rz(-2.4173315) q[3];
sx q[3];
rz(-1.9925502) q[3];
sx q[3];
rz(0.65653062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.637218) q[0];
sx q[0];
rz(-1.1906304) q[0];
sx q[0];
rz(0.36886886) q[0];
rz(0.24636191) q[1];
sx q[1];
rz(-1.3280222) q[1];
sx q[1];
rz(-1.7074283) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1808162) q[0];
sx q[0];
rz(-2.3276276) q[0];
sx q[0];
rz(-1.5778944) q[0];
rz(-pi) q[1];
rz(-0.35583115) q[2];
sx q[2];
rz(-0.47105481) q[2];
sx q[2];
rz(1.6182773) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7692657) q[1];
sx q[1];
rz(-2.7397836) q[1];
sx q[1];
rz(-3.1308082) q[1];
x q[2];
rz(1.3433387) q[3];
sx q[3];
rz(-1.47336) q[3];
sx q[3];
rz(-2.1506598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.65115923) q[2];
sx q[2];
rz(-1.8305402) q[2];
sx q[2];
rz(-1.0820214) q[2];
rz(-2.3467482) q[3];
sx q[3];
rz(-1.669603) q[3];
sx q[3];
rz(1.2580416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.865888) q[0];
sx q[0];
rz(-3.0401433) q[0];
sx q[0];
rz(-0.24359447) q[0];
rz(-0.98681915) q[1];
sx q[1];
rz(-2.3934264) q[1];
sx q[1];
rz(-2.2056244) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075386062) q[0];
sx q[0];
rz(-0.52156007) q[0];
sx q[0];
rz(0.88514502) q[0];
rz(-pi) q[1];
rz(-1.1963821) q[2];
sx q[2];
rz(-2.3585359) q[2];
sx q[2];
rz(-0.20314344) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82560357) q[1];
sx q[1];
rz(-0.88615075) q[1];
sx q[1];
rz(1.0427834) q[1];
rz(-pi) q[2];
rz(-1.5395079) q[3];
sx q[3];
rz(-0.60006053) q[3];
sx q[3];
rz(2.4158583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4794856) q[2];
sx q[2];
rz(-2.0027436) q[2];
sx q[2];
rz(-0.34995079) q[2];
rz(2.5838666) q[3];
sx q[3];
rz(-1.2099268) q[3];
sx q[3];
rz(0.85132712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6996985) q[0];
sx q[0];
rz(-0.63755578) q[0];
sx q[0];
rz(-0.30174524) q[0];
rz(2.8217577) q[1];
sx q[1];
rz(-1.6022857) q[1];
sx q[1];
rz(1.1357657) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2280884) q[0];
sx q[0];
rz(-2.4948409) q[0];
sx q[0];
rz(-0.089368377) q[0];
rz(-2.202353) q[2];
sx q[2];
rz(-1.9983091) q[2];
sx q[2];
rz(-1.7981426) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3150683) q[1];
sx q[1];
rz(-2.5365735) q[1];
sx q[1];
rz(1.7111163) q[1];
rz(-3.0301827) q[3];
sx q[3];
rz(-0.79395959) q[3];
sx q[3];
rz(0.68796989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4096421) q[2];
sx q[2];
rz(-0.66764098) q[2];
sx q[2];
rz(-0.14275924) q[2];
rz(1.7536633) q[3];
sx q[3];
rz(-1.1889435) q[3];
sx q[3];
rz(-2.4587542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1801572) q[0];
sx q[0];
rz(-0.016594369) q[0];
sx q[0];
rz(-0.55602443) q[0];
rz(2.9122638) q[1];
sx q[1];
rz(-1.8792968) q[1];
sx q[1];
rz(-1.75846) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62263238) q[0];
sx q[0];
rz(-0.83370249) q[0];
sx q[0];
rz(1.6181158) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6298619) q[2];
sx q[2];
rz(-2.0286273) q[2];
sx q[2];
rz(0.77564592) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.38249967) q[1];
sx q[1];
rz(-1.9803932) q[1];
sx q[1];
rz(-2.8252557) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9541627) q[3];
sx q[3];
rz(-1.1463506) q[3];
sx q[3];
rz(-0.92680537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.51656276) q[2];
sx q[2];
rz(-0.27250686) q[2];
sx q[2];
rz(1.2410835) q[2];
rz(0.062601335) q[3];
sx q[3];
rz(-1.4227941) q[3];
sx q[3];
rz(0.82714287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0144219) q[0];
sx q[0];
rz(-1.7272471) q[0];
sx q[0];
rz(0.14778368) q[0];
rz(0.5087018) q[1];
sx q[1];
rz(-1.3721507) q[1];
sx q[1];
rz(-0.37857372) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7258647) q[0];
sx q[0];
rz(-1.3368538) q[0];
sx q[0];
rz(-0.42629231) q[0];
rz(0.97855391) q[2];
sx q[2];
rz(-1.2662953) q[2];
sx q[2];
rz(-1.026012) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3562153) q[1];
sx q[1];
rz(-1.2357724) q[1];
sx q[1];
rz(-2.8122693) q[1];
rz(-0.76102961) q[3];
sx q[3];
rz(-1.6609523) q[3];
sx q[3];
rz(2.6554299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.96616894) q[2];
sx q[2];
rz(-0.95169008) q[2];
sx q[2];
rz(1.8625205) q[2];
rz(-1.7133948) q[3];
sx q[3];
rz(-2.1396075) q[3];
sx q[3];
rz(-2.9960347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3928423) q[0];
sx q[0];
rz(-2.5635283) q[0];
sx q[0];
rz(3.0774935) q[0];
rz(0.52866689) q[1];
sx q[1];
rz(-1.268498) q[1];
sx q[1];
rz(0.97506964) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3885766) q[0];
sx q[0];
rz(-2.574769) q[0];
sx q[0];
rz(-0.50811572) q[0];
x q[1];
rz(2.5717402) q[2];
sx q[2];
rz(-2.0790711) q[2];
sx q[2];
rz(0.13392553) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2960423) q[1];
sx q[1];
rz(-0.49440372) q[1];
sx q[1];
rz(0.47459666) q[1];
rz(-pi) q[2];
rz(-0.75210877) q[3];
sx q[3];
rz(-0.32882133) q[3];
sx q[3];
rz(0.11550918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9958682) q[2];
sx q[2];
rz(-2.4805562) q[2];
sx q[2];
rz(0.60689849) q[2];
rz(0.25035826) q[3];
sx q[3];
rz(-1.4162049) q[3];
sx q[3];
rz(-2.6355766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.09457) q[0];
sx q[0];
rz(-1.2170412) q[0];
sx q[0];
rz(1.8893597) q[0];
rz(-2.6429214) q[1];
sx q[1];
rz(-1.55232) q[1];
sx q[1];
rz(1.3943863) q[1];
rz(2.2398938) q[2];
sx q[2];
rz(-2.31291) q[2];
sx q[2];
rz(-0.55533021) q[2];
rz(-1.3648894) q[3];
sx q[3];
rz(-1.6905224) q[3];
sx q[3];
rz(-0.57991309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
