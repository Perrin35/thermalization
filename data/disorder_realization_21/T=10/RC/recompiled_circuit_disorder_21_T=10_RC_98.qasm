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
rz(-0.50146377) q[0];
rz(1.4986829) q[1];
sx q[1];
rz(-2.745435) q[1];
sx q[1];
rz(-0.3224386) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8034536) q[0];
sx q[0];
rz(-0.33284602) q[0];
sx q[0];
rz(-1.3071609) q[0];
x q[1];
rz(-2.2519977) q[2];
sx q[2];
rz(-2.0686364) q[2];
sx q[2];
rz(-1.7878469) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7079561) q[1];
sx q[1];
rz(-2.4383713) q[1];
sx q[1];
rz(-1.0183079) q[1];
rz(-2.1665855) q[3];
sx q[3];
rz(-2.4132204) q[3];
sx q[3];
rz(-1.4260074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.50513187) q[2];
sx q[2];
rz(-2.5487066) q[2];
sx q[2];
rz(0.55603975) q[2];
rz(-0.83267823) q[3];
sx q[3];
rz(-1.4913538) q[3];
sx q[3];
rz(0.94579831) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44822025) q[0];
sx q[0];
rz(-1.6813261) q[0];
sx q[0];
rz(-2.9843176) q[0];
rz(2.8804624) q[1];
sx q[1];
rz(-1.7938679) q[1];
sx q[1];
rz(0.10903407) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.581668) q[0];
sx q[0];
rz(-1.4002869) q[0];
sx q[0];
rz(1.8018363) q[0];
rz(-pi) q[1];
rz(0.63983812) q[2];
sx q[2];
rz(-2.8176762) q[2];
sx q[2];
rz(1.6537635) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8728767) q[1];
sx q[1];
rz(-2.0375588) q[1];
sx q[1];
rz(0.66653911) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29626366) q[3];
sx q[3];
rz(-2.6014572) q[3];
sx q[3];
rz(1.6903282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2382425) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(1.2634574) q[2];
rz(-2.8144828) q[3];
sx q[3];
rz(-1.5771022) q[3];
sx q[3];
rz(1.2143149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0771714) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(-1.7984614) q[0];
rz(-2.893977) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(0.48167357) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4108698) q[0];
sx q[0];
rz(-2.0559089) q[0];
sx q[0];
rz(0.42318015) q[0];
rz(-2.0273676) q[2];
sx q[2];
rz(-0.66514981) q[2];
sx q[2];
rz(0.84012023) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2783918) q[1];
sx q[1];
rz(-2.26483) q[1];
sx q[1];
rz(-0.73514003) q[1];
rz(-pi) q[2];
rz(2.2504911) q[3];
sx q[3];
rz(-2.292233) q[3];
sx q[3];
rz(-0.47618714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3383011) q[2];
sx q[2];
rz(-0.81739601) q[2];
sx q[2];
rz(0.49989191) q[2];
rz(-0.56097427) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(-1.6803754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9445779) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(-2.5894077) q[0];
rz(-1.5532956) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(1.8968556) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9024076) q[0];
sx q[0];
rz(-2.8021078) q[0];
sx q[0];
rz(-3.1367338) q[0];
rz(-pi) q[1];
rz(-0.130529) q[2];
sx q[2];
rz(-1.7703238) q[2];
sx q[2];
rz(-0.56632698) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2326395) q[1];
sx q[1];
rz(-0.5852355) q[1];
sx q[1];
rz(1.1802243) q[1];
rz(-1.3454868) q[3];
sx q[3];
rz(-2.1300737) q[3];
sx q[3];
rz(-0.51636458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.84919471) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(-1.1506895) q[2];
rz(1.6644647) q[3];
sx q[3];
rz(-1.632558) q[3];
sx q[3];
rz(-2.6586444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078159049) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(-3.0601236) q[0];
rz(-0.062462656) q[1];
sx q[1];
rz(-1.1413347) q[1];
sx q[1];
rz(-1.6385471) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3969288) q[0];
sx q[0];
rz(-1.6984807) q[0];
sx q[0];
rz(2.1555107) q[0];
rz(-0.91707768) q[2];
sx q[2];
rz(-0.13609016) q[2];
sx q[2];
rz(-1.2166785) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0283708) q[1];
sx q[1];
rz(-1.2591259) q[1];
sx q[1];
rz(-3.034364) q[1];
x q[2];
rz(1.2522069) q[3];
sx q[3];
rz(-1.0831523) q[3];
sx q[3];
rz(0.82061758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.5082671) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(1.2379237) q[2];
rz(-1.1226908) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(2.6200263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5313107) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(-2.8748728) q[0];
rz(-0.56089127) q[1];
sx q[1];
rz(-1.2979049) q[1];
sx q[1];
rz(0.7985324) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7245367) q[0];
sx q[0];
rz(-0.31053156) q[0];
sx q[0];
rz(2.0470371) q[0];
x q[1];
rz(2.550225) q[2];
sx q[2];
rz(-0.57833507) q[2];
sx q[2];
rz(-0.27507281) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9763899) q[1];
sx q[1];
rz(-2.5678647) q[1];
sx q[1];
rz(-1.4742673) q[1];
rz(-0.800662) q[3];
sx q[3];
rz(-1.1122397) q[3];
sx q[3];
rz(2.3199905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.16053998) q[2];
sx q[2];
rz(-1.2489698) q[2];
sx q[2];
rz(2.7748761) q[2];
rz(1.2612873) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(0.096207531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-0.60638705) q[0];
rz(-2.9442893) q[1];
sx q[1];
rz(-2.0154672) q[1];
sx q[1];
rz(-2.6775449) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2460829) q[0];
sx q[0];
rz(-1.9946788) q[0];
sx q[0];
rz(-0.8830107) q[0];
rz(-pi) q[1];
rz(1.4887965) q[2];
sx q[2];
rz(-2.900015) q[2];
sx q[2];
rz(-1.8919924) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5210919) q[1];
sx q[1];
rz(-1.2695841) q[1];
sx q[1];
rz(-2.0786933) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9052827) q[3];
sx q[3];
rz(-1.650562) q[3];
sx q[3];
rz(-0.41111708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.93418926) q[2];
sx q[2];
rz(-1.0031909) q[2];
sx q[2];
rz(-0.25804538) q[2];
rz(-1.1856273) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(-3.1380838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.946452) q[0];
sx q[0];
rz(-1.8608681) q[0];
sx q[0];
rz(-2.7602957) q[0];
rz(0.095245846) q[1];
sx q[1];
rz(-0.97243273) q[1];
sx q[1];
rz(1.4415178) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7521116) q[0];
sx q[0];
rz(-1.6342667) q[0];
sx q[0];
rz(1.585929) q[0];
rz(2.5289815) q[2];
sx q[2];
rz(-1.5393886) q[2];
sx q[2];
rz(2.999246) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4691094) q[1];
sx q[1];
rz(-1.6961349) q[1];
sx q[1];
rz(-0.13087665) q[1];
x q[2];
rz(-0.12318792) q[3];
sx q[3];
rz(-0.33629575) q[3];
sx q[3];
rz(-1.5598701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1429446) q[2];
sx q[2];
rz(-2.7286077) q[2];
sx q[2];
rz(2.9150035) q[2];
rz(-0.44858027) q[3];
sx q[3];
rz(-1.5357163) q[3];
sx q[3];
rz(-2.3118238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7609693) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(-1.3990078) q[0];
rz(-0.31708583) q[1];
sx q[1];
rz(-1.6665019) q[1];
sx q[1];
rz(2.1549966) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87514585) q[0];
sx q[0];
rz(-2.6212924) q[0];
sx q[0];
rz(-2.8199024) q[0];
rz(-pi) q[1];
rz(0.53194745) q[2];
sx q[2];
rz(-0.94594687) q[2];
sx q[2];
rz(-2.9192386) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4793195) q[1];
sx q[1];
rz(-2.5179177) q[1];
sx q[1];
rz(-0.24346607) q[1];
x q[2];
rz(-1.9916233) q[3];
sx q[3];
rz(-1.4126561) q[3];
sx q[3];
rz(-1.6588253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2150779) q[2];
sx q[2];
rz(-2.4145917) q[2];
sx q[2];
rz(0.40965664) q[2];
rz(-2.8783197) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(-0.51945654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7423994) q[0];
sx q[0];
rz(-0.078646794) q[0];
sx q[0];
rz(-1.4051399) q[0];
rz(-2.3204904) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(-1.7260889) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50354276) q[0];
sx q[0];
rz(-1.8068131) q[0];
sx q[0];
rz(0.34061265) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.892898) q[2];
sx q[2];
rz(-1.2173614) q[2];
sx q[2];
rz(-0.24662185) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7064757) q[1];
sx q[1];
rz(-0.26285989) q[1];
sx q[1];
rz(0.95107066) q[1];
rz(-pi) q[2];
rz(-2.0831765) q[3];
sx q[3];
rz(-1.4231764) q[3];
sx q[3];
rz(-0.83422134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4225509) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(3.0174875) q[2];
rz(-0.96578807) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(0.66108274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60349764) q[0];
sx q[0];
rz(-2.8932543) q[0];
sx q[0];
rz(2.2809991) q[0];
rz(-0.30766906) q[1];
sx q[1];
rz(-1.2528906) q[1];
sx q[1];
rz(1.2045592) q[1];
rz(-0.60097354) q[2];
sx q[2];
rz(-2.8291611) q[2];
sx q[2];
rz(-0.57401382) q[2];
rz(-0.55660558) q[3];
sx q[3];
rz(-1.442933) q[3];
sx q[3];
rz(1.8419151) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
