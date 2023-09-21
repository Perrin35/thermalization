OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(-3.0598109) q[0];
sx q[0];
rz(2.6401289) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(0.3224386) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8034536) q[0];
sx q[0];
rz(-0.33284602) q[0];
sx q[0];
rz(-1.8344318) q[0];
x q[1];
rz(-2.2519977) q[2];
sx q[2];
rz(-1.0729562) q[2];
sx q[2];
rz(-1.3537458) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7079561) q[1];
sx q[1];
rz(-0.70322137) q[1];
sx q[1];
rz(1.0183079) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6775042) q[3];
sx q[3];
rz(-0.9872735) q[3];
sx q[3];
rz(-2.1634963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6364608) q[2];
sx q[2];
rz(-2.5487066) q[2];
sx q[2];
rz(0.55603975) q[2];
rz(2.3089144) q[3];
sx q[3];
rz(-1.4913538) q[3];
sx q[3];
rz(0.94579831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44822025) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(-0.15727501) q[0];
rz(-0.26113025) q[1];
sx q[1];
rz(-1.3477247) q[1];
sx q[1];
rz(3.0325586) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63576525) q[0];
sx q[0];
rz(-0.286239) q[0];
sx q[0];
rz(2.2155227) q[0];
x q[1];
rz(2.5017545) q[2];
sx q[2];
rz(-2.8176762) q[2];
sx q[2];
rz(-1.6537635) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2687159) q[1];
sx q[1];
rz(-1.1040338) q[1];
sx q[1];
rz(-2.4750535) q[1];
rz(-2.845329) q[3];
sx q[3];
rz(-2.6014572) q[3];
sx q[3];
rz(-1.6903282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2382425) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(-1.2634574) q[2];
rz(2.8144828) q[3];
sx q[3];
rz(-1.5644904) q[3];
sx q[3];
rz(1.2143149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0771714) q[0];
sx q[0];
rz(-3.0922958) q[0];
sx q[0];
rz(1.3431312) q[0];
rz(0.24761565) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(-2.6599191) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7307229) q[0];
sx q[0];
rz(-1.0856837) q[0];
sx q[0];
rz(2.7184125) q[0];
rz(-pi) q[1];
rz(-2.8086497) q[2];
sx q[2];
rz(-0.98368401) q[2];
sx q[2];
rz(1.3981896) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.18322769) q[1];
sx q[1];
rz(-2.1124358) q[1];
sx q[1];
rz(-0.72802131) q[1];
x q[2];
rz(-0.89110156) q[3];
sx q[3];
rz(-0.84935969) q[3];
sx q[3];
rz(-2.6654055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8032916) q[2];
sx q[2];
rz(-2.3241966) q[2];
sx q[2];
rz(2.6417007) q[2];
rz(-2.5806184) q[3];
sx q[3];
rz(-1.2596954) q[3];
sx q[3];
rz(1.4612173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9445779) q[0];
sx q[0];
rz(-0.16600969) q[0];
sx q[0];
rz(-2.5894077) q[0];
rz(1.588297) q[1];
sx q[1];
rz(-2.242656) q[1];
sx q[1];
rz(1.2447371) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8053999) q[0];
sx q[0];
rz(-1.5724143) q[0];
sx q[0];
rz(-0.33948116) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99889465) q[2];
sx q[2];
rz(-0.23795393) q[2];
sx q[2];
rz(-1.1513125) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9089531) q[1];
sx q[1];
rz(-2.5563572) q[1];
sx q[1];
rz(-1.1802243) q[1];
x q[2];
rz(2.5707385) q[3];
sx q[3];
rz(-1.7613162) q[3];
sx q[3];
rz(1.1754456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2923979) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(1.1506895) q[2];
rz(-1.4771279) q[3];
sx q[3];
rz(-1.5090347) q[3];
sx q[3];
rz(2.6586444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-3.0634336) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(-3.0601236) q[0];
rz(3.07913) q[1];
sx q[1];
rz(-2.000258) q[1];
sx q[1];
rz(-1.5030456) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7446639) q[0];
sx q[0];
rz(-1.6984807) q[0];
sx q[0];
rz(-2.1555107) q[0];
rz(-pi) q[1];
rz(0.91707768) q[2];
sx q[2];
rz(-0.13609016) q[2];
sx q[2];
rz(-1.9249141) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0283708) q[1];
sx q[1];
rz(-1.2591259) q[1];
sx q[1];
rz(-3.034364) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50932192) q[3];
sx q[3];
rz(-1.29042) q[3];
sx q[3];
rz(2.544739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5082671) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(-1.2379237) q[2];
rz(2.0189019) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(2.6200263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6102819) q[0];
sx q[0];
rz(-1.0045445) q[0];
sx q[0];
rz(-2.8748728) q[0];
rz(-0.56089127) q[1];
sx q[1];
rz(-1.2979049) q[1];
sx q[1];
rz(-2.3430603) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.610299) q[0];
sx q[0];
rz(-1.4302505) q[0];
sx q[0];
rz(1.848624) q[0];
rz(-pi) q[1];
rz(0.59136765) q[2];
sx q[2];
rz(-0.57833507) q[2];
sx q[2];
rz(-2.8665198) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.861607) q[1];
sx q[1];
rz(-2.1415188) q[1];
sx q[1];
rz(3.0793889) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5391606) q[3];
sx q[3];
rz(-2.2450387) q[3];
sx q[3];
rz(-0.34365052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16053998) q[2];
sx q[2];
rz(-1.8926228) q[2];
sx q[2];
rz(-0.36671656) q[2];
rz(1.8803053) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(-0.096207531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7325608) q[0];
sx q[0];
rz(-0.92027396) q[0];
sx q[0];
rz(2.5352056) q[0];
rz(2.9442893) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(0.46404776) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34940091) q[0];
sx q[0];
rz(-2.1878625) q[0];
sx q[0];
rz(2.6130136) q[0];
rz(-0.020178528) q[2];
sx q[2];
rz(-1.3300465) q[2];
sx q[2];
rz(1.1651595) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6021767) q[1];
sx q[1];
rz(-0.5836986) q[1];
sx q[1];
rz(1.0023487) q[1];
rz(2.9052827) q[3];
sx q[3];
rz(-1.650562) q[3];
sx q[3];
rz(0.41111708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.93418926) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(0.25804538) q[2];
rz(-1.9559654) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(-0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19514062) q[0];
sx q[0];
rz(-1.8608681) q[0];
sx q[0];
rz(2.7602957) q[0];
rz(3.0463468) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(-1.7000748) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7521116) q[0];
sx q[0];
rz(-1.6342667) q[0];
sx q[0];
rz(1.5556637) q[0];
x q[1];
rz(-0.054585329) q[2];
sx q[2];
rz(-2.5282801) q[2];
sx q[2];
rz(1.7577946) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.86130202) q[1];
sx q[1];
rz(-0.18096563) q[1];
sx q[1];
rz(0.7678395) q[1];
rz(-pi) q[2];
rz(0.12318792) q[3];
sx q[3];
rz(-0.33629575) q[3];
sx q[3];
rz(1.5598701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9986481) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(2.9150035) q[2];
rz(-2.6930124) q[3];
sx q[3];
rz(-1.5357163) q[3];
sx q[3];
rz(-0.8297689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3806234) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(1.7425849) q[0];
rz(-2.8245068) q[1];
sx q[1];
rz(-1.4750907) q[1];
sx q[1];
rz(2.1549966) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87514585) q[0];
sx q[0];
rz(-0.52030021) q[0];
sx q[0];
rz(-0.32169028) q[0];
rz(-2.2676342) q[2];
sx q[2];
rz(-1.1468337) q[2];
sx q[2];
rz(1.4615814) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1074859) q[1];
sx q[1];
rz(-1.7120546) q[1];
sx q[1];
rz(-0.60955255) q[1];
rz(0.17296965) q[3];
sx q[3];
rz(-1.9860455) q[3];
sx q[3];
rz(-3.1239307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9265147) q[2];
sx q[2];
rz(-2.4145917) q[2];
sx q[2];
rz(-2.731936) q[2];
rz(2.8783197) q[3];
sx q[3];
rz(-1.8299088) q[3];
sx q[3];
rz(-0.51945654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7423994) q[0];
sx q[0];
rz(-0.078646794) q[0];
sx q[0];
rz(1.4051399) q[0];
rz(-0.82110226) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(1.7260889) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1499407) q[0];
sx q[0];
rz(-1.2399925) q[0];
sx q[0];
rz(1.8206235) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2070933) q[2];
sx q[2];
rz(-1.3377681) q[2];
sx q[2];
rz(-1.9050913) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6092097) q[1];
sx q[1];
rz(-1.4192974) q[1];
sx q[1];
rz(-1.3551559) q[1];
rz(-pi) q[2];
rz(-1.2763001) q[3];
sx q[3];
rz(-0.53139549) q[3];
sx q[3];
rz(2.6607799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71904174) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(-3.0174875) q[2];
rz(0.96578807) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(-0.66108274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60349764) q[0];
sx q[0];
rz(-2.8932543) q[0];
sx q[0];
rz(2.2809991) q[0];
rz(0.30766906) q[1];
sx q[1];
rz(-1.888702) q[1];
sx q[1];
rz(-1.9370334) q[1];
rz(-0.260367) q[2];
sx q[2];
rz(-1.7454864) q[2];
sx q[2];
rz(-2.7228552) q[2];
rz(-0.23871213) q[3];
sx q[3];
rz(-0.56959) q[3];
sx q[3];
rz(0.068985229) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
