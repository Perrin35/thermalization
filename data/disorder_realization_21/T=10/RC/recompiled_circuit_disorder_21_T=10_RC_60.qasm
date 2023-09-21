OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0857467) q[0];
sx q[0];
rz(-0.081781713) q[0];
sx q[0];
rz(-2.6401289) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(-2.8191541) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8034536) q[0];
sx q[0];
rz(-2.8087466) q[0];
sx q[0];
rz(-1.8344318) q[0];
rz(2.2828322) q[2];
sx q[2];
rz(-2.3220064) q[2];
sx q[2];
rz(2.8263) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4336366) q[1];
sx q[1];
rz(-2.4383713) q[1];
sx q[1];
rz(-2.1232848) q[1];
rz(-0.97500719) q[3];
sx q[3];
rz(-2.4132204) q[3];
sx q[3];
rz(1.4260074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.50513187) q[2];
sx q[2];
rz(-2.5487066) q[2];
sx q[2];
rz(0.55603975) q[2];
rz(-2.3089144) q[3];
sx q[3];
rz(-1.6502389) q[3];
sx q[3];
rz(0.94579831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6933724) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(0.15727501) q[0];
rz(2.8804624) q[1];
sx q[1];
rz(-1.7938679) q[1];
sx q[1];
rz(-3.0325586) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029023829) q[0];
sx q[0];
rz(-1.3431664) q[0];
sx q[0];
rz(2.9665222) q[0];
rz(-pi) q[1];
rz(-1.7686339) q[2];
sx q[2];
rz(-1.3126144) q[2];
sx q[2];
rz(-2.3193662) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2687159) q[1];
sx q[1];
rz(-1.1040338) q[1];
sx q[1];
rz(0.66653911) q[1];
rz(-pi) q[2];
rz(1.7440967) q[3];
sx q[3];
rz(-1.0565851) q[3];
sx q[3];
rz(1.7931995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9033501) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(-1.2634574) q[2];
rz(-0.3271099) q[3];
sx q[3];
rz(-1.5771022) q[3];
sx q[3];
rz(1.9272778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.064421244) q[0];
sx q[0];
rz(-3.0922958) q[0];
sx q[0];
rz(-1.3431312) q[0];
rz(-2.893977) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(0.48167357) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047086296) q[0];
sx q[0];
rz(-1.1990093) q[0];
sx q[0];
rz(-2.0949754) q[0];
rz(-pi) q[1];
rz(0.9573612) q[2];
sx q[2];
rz(-1.8463677) q[2];
sx q[2];
rz(-2.7797109) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18322769) q[1];
sx q[1];
rz(-1.0291568) q[1];
sx q[1];
rz(-0.72802131) q[1];
x q[2];
rz(2.5211469) q[3];
sx q[3];
rz(-2.1944322) q[3];
sx q[3];
rz(0.40943957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8032916) q[2];
sx q[2];
rz(-0.81739601) q[2];
sx q[2];
rz(0.49989191) q[2];
rz(2.5806184) q[3];
sx q[3];
rz(-1.2596954) q[3];
sx q[3];
rz(-1.4612173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9445779) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(2.5894077) q[0];
rz(-1.5532956) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(-1.2447371) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23403215) q[0];
sx q[0];
rz(-1.910277) q[0];
sx q[0];
rz(1.5725122) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3696026) q[2];
sx q[2];
rz(-1.6987213) q[2];
sx q[2];
rz(-2.163137) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9089531) q[1];
sx q[1];
rz(-0.5852355) q[1];
sx q[1];
rz(1.9613683) q[1];
rz(-pi) q[2];
rz(-2.7987715) q[3];
sx q[3];
rz(-0.59844136) q[3];
sx q[3];
rz(3.0330021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.84919471) q[2];
sx q[2];
rz(-1.8820102) q[2];
sx q[2];
rz(-1.9909031) q[2];
rz(1.6644647) q[3];
sx q[3];
rz(-1.632558) q[3];
sx q[3];
rz(0.48294827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(3.0634336) q[0];
sx q[0];
rz(-2.3796191) q[0];
sx q[0];
rz(0.081469014) q[0];
rz(-3.07913) q[1];
sx q[1];
rz(-1.1413347) q[1];
sx q[1];
rz(1.6385471) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3969288) q[0];
sx q[0];
rz(-1.443112) q[0];
sx q[0];
rz(-0.98608195) q[0];
rz(-pi) q[1];
rz(-0.91707768) q[2];
sx q[2];
rz(-0.13609016) q[2];
sx q[2];
rz(1.9249141) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4245783) q[1];
sx q[1];
rz(-1.6728405) q[1];
sx q[1];
rz(-1.8841519) q[1];
rz(-pi) q[2];
rz(-1.8893858) q[3];
sx q[3];
rz(-2.0584403) q[3];
sx q[3];
rz(-0.82061758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.5082671) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(1.2379237) q[2];
rz(-2.0189019) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(0.52156633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5313107) q[0];
sx q[0];
rz(-1.0045445) q[0];
sx q[0];
rz(2.8748728) q[0];
rz(-0.56089127) q[1];
sx q[1];
rz(-1.2979049) q[1];
sx q[1];
rz(0.7985324) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.610299) q[0];
sx q[0];
rz(-1.4302505) q[0];
sx q[0];
rz(-1.848624) q[0];
rz(-0.59136765) q[2];
sx q[2];
rz(-0.57833507) q[2];
sx q[2];
rz(2.8665198) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.861607) q[1];
sx q[1];
rz(-2.1415188) q[1];
sx q[1];
rz(-3.0793889) q[1];
rz(-0.60243209) q[3];
sx q[3];
rz(-0.89655399) q[3];
sx q[3];
rz(-2.7979421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16053998) q[2];
sx q[2];
rz(-1.8926228) q[2];
sx q[2];
rz(-2.7748761) q[2];
rz(-1.8803053) q[3];
sx q[3];
rz(-1.4533318) q[3];
sx q[3];
rz(-0.096207531) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325608) q[0];
sx q[0];
rz(-2.2213187) q[0];
sx q[0];
rz(2.5352056) q[0];
rz(-2.9442893) q[1];
sx q[1];
rz(-2.0154672) q[1];
sx q[1];
rz(-2.6775449) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34940091) q[0];
sx q[0];
rz(-0.95373017) q[0];
sx q[0];
rz(0.52857907) q[0];
rz(-pi) q[1];
rz(3.1214141) q[2];
sx q[2];
rz(-1.3300465) q[2];
sx q[2];
rz(-1.9764331) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0276427) q[1];
sx q[1];
rz(-2.0538035) q[1];
sx q[1];
rz(-0.34160683) q[1];
rz(-2.9052827) q[3];
sx q[3];
rz(-1.650562) q[3];
sx q[3];
rz(2.7304756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.93418926) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(2.8835473) q[2];
rz(-1.1856273) q[3];
sx q[3];
rz(-1.6262755) q[3];
sx q[3];
rz(-0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19514062) q[0];
sx q[0];
rz(-1.8608681) q[0];
sx q[0];
rz(-0.38129693) q[0];
rz(-3.0463468) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(-1.4415178) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1822752) q[0];
sx q[0];
rz(-1.5556941) q[0];
sx q[0];
rz(-3.0781151) q[0];
x q[1];
rz(-0.61261119) q[2];
sx q[2];
rz(-1.5393886) q[2];
sx q[2];
rz(2.999246) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.085233363) q[1];
sx q[1];
rz(-1.4409522) q[1];
sx q[1];
rz(1.4443881) q[1];
rz(1.5278682) q[3];
sx q[3];
rz(-1.9044442) q[3];
sx q[3];
rz(-1.4294525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9986481) q[2];
sx q[2];
rz(-2.7286077) q[2];
sx q[2];
rz(0.22658919) q[2];
rz(2.6930124) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(2.3118238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3806234) q[0];
sx q[0];
rz(-2.3801104) q[0];
sx q[0];
rz(1.7425849) q[0];
rz(-2.8245068) q[1];
sx q[1];
rz(-1.4750907) q[1];
sx q[1];
rz(2.1549966) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6331659) q[0];
sx q[0];
rz(-2.0619443) q[0];
sx q[0];
rz(1.3915865) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87395845) q[2];
sx q[2];
rz(-1.9947589) q[2];
sx q[2];
rz(1.6800113) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.36531891) q[1];
sx q[1];
rz(-0.96818189) q[1];
sx q[1];
rz(-1.7425294) q[1];
rz(-1.9429728) q[3];
sx q[3];
rz(-2.6937006) q[3];
sx q[3];
rz(-2.7152293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2150779) q[2];
sx q[2];
rz(-0.72700095) q[2];
sx q[2];
rz(-2.731936) q[2];
rz(-0.26327291) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(0.51945654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39919329) q[0];
sx q[0];
rz(-0.078646794) q[0];
sx q[0];
rz(1.4051399) q[0];
rz(0.82110226) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(-1.7260889) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9916519) q[0];
sx q[0];
rz(-1.2399925) q[0];
sx q[0];
rz(-1.8206235) q[0];
rz(-1.9344994) q[2];
sx q[2];
rz(-1.8038245) q[2];
sx q[2];
rz(1.9050913) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.43511697) q[1];
sx q[1];
rz(-2.8787328) q[1];
sx q[1];
rz(-0.95107066) q[1];
rz(1.2763001) q[3];
sx q[3];
rz(-2.6101972) q[3];
sx q[3];
rz(2.6607799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.71904174) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(-0.12410513) q[2];
rz(-2.1758046) q[3];
sx q[3];
rz(-1.6671168) q[3];
sx q[3];
rz(-2.4805099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60349764) q[0];
sx q[0];
rz(-2.8932543) q[0];
sx q[0];
rz(2.2809991) q[0];
rz(2.8339236) q[1];
sx q[1];
rz(-1.2528906) q[1];
sx q[1];
rz(1.2045592) q[1];
rz(1.3901426) q[2];
sx q[2];
rz(-1.3144819) q[2];
sx q[2];
rz(-1.1983295) q[2];
rz(0.23871213) q[3];
sx q[3];
rz(-2.5720027) q[3];
sx q[3];
rz(-3.0726074) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
