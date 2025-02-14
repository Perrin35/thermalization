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
rz(0.60740745) q[0];
sx q[0];
rz(-0.4439126) q[0];
sx q[0];
rz(-2.8191415) q[0];
rz(-1.9880265) q[1];
sx q[1];
rz(-1.1918951) q[1];
sx q[1];
rz(0.21790394) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26051021) q[0];
sx q[0];
rz(-0.54303193) q[0];
sx q[0];
rz(0.61123993) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.430021) q[2];
sx q[2];
rz(-2.8970538) q[2];
sx q[2];
rz(-1.0734347) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.17927781) q[1];
sx q[1];
rz(-1.4554108) q[1];
sx q[1];
rz(-0.25956599) q[1];
rz(2.3090906) q[3];
sx q[3];
rz(-2.7232086) q[3];
sx q[3];
rz(0.42319059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0107062) q[2];
sx q[2];
rz(-2.71038) q[2];
sx q[2];
rz(0.19403379) q[2];
rz(2.7583097) q[3];
sx q[3];
rz(-0.88519874) q[3];
sx q[3];
rz(-1.6398199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5782769) q[0];
sx q[0];
rz(-1.0551772) q[0];
sx q[0];
rz(2.03736) q[0];
rz(-0.68179321) q[1];
sx q[1];
rz(-1.7742523) q[1];
sx q[1];
rz(1.404748) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.800398) q[0];
sx q[0];
rz(-1.826108) q[0];
sx q[0];
rz(-1.4394151) q[0];
rz(-1.0114298) q[2];
sx q[2];
rz(-2.1678143) q[2];
sx q[2];
rz(1.9950713) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3225783) q[1];
sx q[1];
rz(-2.0657263) q[1];
sx q[1];
rz(1.9274345) q[1];
rz(-pi) q[2];
rz(-1.6854755) q[3];
sx q[3];
rz(-1.3175829) q[3];
sx q[3];
rz(0.51451433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3746987) q[2];
sx q[2];
rz(-0.99504343) q[2];
sx q[2];
rz(1.2332756) q[2];
rz(2.342566) q[3];
sx q[3];
rz(-0.34501758) q[3];
sx q[3];
rz(-0.12921216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98053539) q[0];
sx q[0];
rz(-2.0687456) q[0];
sx q[0];
rz(0.55915731) q[0];
rz(-1.1178389) q[1];
sx q[1];
rz(-1.585377) q[1];
sx q[1];
rz(0.1112172) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1532899) q[0];
sx q[0];
rz(-2.7493101) q[0];
sx q[0];
rz(0.28052237) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8410013) q[2];
sx q[2];
rz(-2.1779354) q[2];
sx q[2];
rz(2.6532113) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.049455482) q[1];
sx q[1];
rz(-2.4087445) q[1];
sx q[1];
rz(1.2867828) q[1];
rz(-pi) q[2];
rz(-2.7989332) q[3];
sx q[3];
rz(-2.7268134) q[3];
sx q[3];
rz(0.77769731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3126276) q[2];
sx q[2];
rz(-1.6964922) q[2];
sx q[2];
rz(1.1980537) q[2];
rz(-1.4231921) q[3];
sx q[3];
rz(-1.4462024) q[3];
sx q[3];
rz(-2.6660582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51197416) q[0];
sx q[0];
rz(-2.8422575) q[0];
sx q[0];
rz(-2.3053115) q[0];
rz(-1.5654927) q[1];
sx q[1];
rz(-1.0944159) q[1];
sx q[1];
rz(-2.2907168) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99573409) q[0];
sx q[0];
rz(-0.99829095) q[0];
sx q[0];
rz(-0.03972223) q[0];
rz(0.81340547) q[2];
sx q[2];
rz(-2.4508173) q[2];
sx q[2];
rz(-1.9936313) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5454514) q[1];
sx q[1];
rz(-0.79789466) q[1];
sx q[1];
rz(1.376838) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9662764) q[3];
sx q[3];
rz(-1.3787601) q[3];
sx q[3];
rz(-0.78498299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1026844) q[2];
sx q[2];
rz(-0.85997283) q[2];
sx q[2];
rz(0.50722185) q[2];
rz(-2.475259) q[3];
sx q[3];
rz(-0.75080502) q[3];
sx q[3];
rz(0.91608086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0136593) q[0];
sx q[0];
rz(-1.6935231) q[0];
sx q[0];
rz(-1.6130945) q[0];
rz(-1.5799892) q[1];
sx q[1];
rz(-2.0785619) q[1];
sx q[1];
rz(-2.9772421) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8678169) q[0];
sx q[0];
rz(-1.4901596) q[0];
sx q[0];
rz(3.0458661) q[0];
rz(-0.63173024) q[2];
sx q[2];
rz(-0.99020489) q[2];
sx q[2];
rz(-2.3931062) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.937035) q[1];
sx q[1];
rz(-1.2271082) q[1];
sx q[1];
rz(2.2930706) q[1];
rz(-pi) q[2];
rz(1.2772395) q[3];
sx q[3];
rz(-1.3428215) q[3];
sx q[3];
rz(-3.0346532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5483755) q[2];
sx q[2];
rz(-0.84725738) q[2];
sx q[2];
rz(0.66620052) q[2];
rz(0.14144746) q[3];
sx q[3];
rz(-1.5513159) q[3];
sx q[3];
rz(-2.1470054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2675466) q[0];
sx q[0];
rz(-0.62804896) q[0];
sx q[0];
rz(0.48126599) q[0];
rz(0.73356837) q[1];
sx q[1];
rz(-1.8736519) q[1];
sx q[1];
rz(3.0487294) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32321445) q[0];
sx q[0];
rz(-1.1284472) q[0];
sx q[0];
rz(-1.044892) q[0];
rz(-pi) q[1];
x q[1];
rz(0.038636908) q[2];
sx q[2];
rz(-1.5728356) q[2];
sx q[2];
rz(0.80261999) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.464092) q[1];
sx q[1];
rz(-0.94225303) q[1];
sx q[1];
rz(1.0586865) q[1];
x q[2];
rz(-2.1867083) q[3];
sx q[3];
rz(-0.43514565) q[3];
sx q[3];
rz(-2.8762115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7218472) q[2];
sx q[2];
rz(-1.6215308) q[2];
sx q[2];
rz(-1.9160371) q[2];
rz(3.0626152) q[3];
sx q[3];
rz(-2.0073399) q[3];
sx q[3];
rz(-1.2823766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6135898) q[0];
sx q[0];
rz(-0.71857518) q[0];
sx q[0];
rz(2.3336616) q[0];
rz(1.5241874) q[1];
sx q[1];
rz(-2.9831191) q[1];
sx q[1];
rz(0.70404109) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8600143) q[0];
sx q[0];
rz(-1.5905955) q[0];
sx q[0];
rz(-3.1327644) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.03116) q[2];
sx q[2];
rz(-3.1322479) q[2];
sx q[2];
rz(1.2197987) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0800533) q[1];
sx q[1];
rz(-2.1916917) q[1];
sx q[1];
rz(0.027959221) q[1];
x q[2];
rz(2.9765997) q[3];
sx q[3];
rz(-2.1793302) q[3];
sx q[3];
rz(-2.7138865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0713221) q[2];
sx q[2];
rz(-0.68838781) q[2];
sx q[2];
rz(0.88453156) q[2];
rz(-0.5136579) q[3];
sx q[3];
rz(-1.1762041) q[3];
sx q[3];
rz(-2.6562712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1270444) q[0];
sx q[0];
rz(-2.3764648) q[0];
sx q[0];
rz(0.86084086) q[0];
rz(-0.035482081) q[1];
sx q[1];
rz(-1.9396962) q[1];
sx q[1];
rz(1.4883581) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4299773) q[0];
sx q[0];
rz(-1.9693144) q[0];
sx q[0];
rz(0.43669711) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8143009) q[2];
sx q[2];
rz(-1.3483682) q[2];
sx q[2];
rz(3.0102863) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.70362857) q[1];
sx q[1];
rz(-1.9537821) q[1];
sx q[1];
rz(-0.18937892) q[1];
rz(-pi) q[2];
rz(-0.32522301) q[3];
sx q[3];
rz(-0.28984335) q[3];
sx q[3];
rz(-2.1397643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.52582467) q[2];
sx q[2];
rz(-1.1349698) q[2];
sx q[2];
rz(-1.5110678) q[2];
rz(-2.3903971) q[3];
sx q[3];
rz(-0.54371756) q[3];
sx q[3];
rz(-2.3403919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5649331) q[0];
sx q[0];
rz(-1.7486005) q[0];
sx q[0];
rz(-0.15889731) q[0];
rz(-1.502602) q[1];
sx q[1];
rz(-1.8105806) q[1];
sx q[1];
rz(1.9749036) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5983655) q[0];
sx q[0];
rz(-1.0002213) q[0];
sx q[0];
rz(-1.2735418) q[0];
x q[1];
rz(1.2449655) q[2];
sx q[2];
rz(-1.836317) q[2];
sx q[2];
rz(0.42761432) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1981956) q[1];
sx q[1];
rz(-1.5875729) q[1];
sx q[1];
rz(-1.6086964) q[1];
rz(2.1871879) q[3];
sx q[3];
rz(-3.1052178) q[3];
sx q[3];
rz(2.006435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4022973) q[2];
sx q[2];
rz(-2.4810956) q[2];
sx q[2];
rz(2.8759549) q[2];
rz(1.1916676) q[3];
sx q[3];
rz(-1.3375514) q[3];
sx q[3];
rz(0.55408293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3091076) q[0];
sx q[0];
rz(-2.5256248) q[0];
sx q[0];
rz(0.68514222) q[0];
rz(-0.55662545) q[1];
sx q[1];
rz(-1.4645422) q[1];
sx q[1];
rz(2.305078) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5852528) q[0];
sx q[0];
rz(-1.9808021) q[0];
sx q[0];
rz(0.37081952) q[0];
rz(2.0903953) q[2];
sx q[2];
rz(-1.8406036) q[2];
sx q[2];
rz(-2.4003911) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3565607) q[1];
sx q[1];
rz(-2.0036657) q[1];
sx q[1];
rz(0.69286107) q[1];
rz(-1.8440866) q[3];
sx q[3];
rz(-1.7215842) q[3];
sx q[3];
rz(2.085872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4461925) q[2];
sx q[2];
rz(-1.7170898) q[2];
sx q[2];
rz(2.2851473) q[2];
rz(2.942318) q[3];
sx q[3];
rz(-2.1592185) q[3];
sx q[3];
rz(2.2299531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0675426) q[0];
sx q[0];
rz(-1.0150801) q[0];
sx q[0];
rz(-1.6652921) q[0];
rz(0.53786565) q[1];
sx q[1];
rz(-1.270351) q[1];
sx q[1];
rz(-1.7958633) q[1];
rz(0.79277586) q[2];
sx q[2];
rz(-2.1332827) q[2];
sx q[2];
rz(2.573043) q[2];
rz(2.2945401) q[3];
sx q[3];
rz(-1.77917) q[3];
sx q[3];
rz(2.8875881) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
