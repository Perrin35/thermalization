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
rz(-3.0109316) q[0];
sx q[0];
rz(-1.8008512) q[0];
sx q[0];
rz(1.8774348) q[0];
rz(0.90509993) q[1];
sx q[1];
rz(2.0533419) q[1];
sx q[1];
rz(6.2702141) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94840996) q[0];
sx q[0];
rz(-1.4140861) q[0];
sx q[0];
rz(-2.66314) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6799203) q[2];
sx q[2];
rz(-1.4089917) q[2];
sx q[2];
rz(1.8617804) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4526635) q[1];
sx q[1];
rz(-0.48349342) q[1];
sx q[1];
rz(-0.74093282) q[1];
x q[2];
rz(1.3566586) q[3];
sx q[3];
rz(-0.99839568) q[3];
sx q[3];
rz(-2.1635087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26244792) q[2];
sx q[2];
rz(-1.6877354) q[2];
sx q[2];
rz(-3.0461779) q[2];
rz(-2.6573507) q[3];
sx q[3];
rz(-0.80317322) q[3];
sx q[3];
rz(-1.6835083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5937623) q[0];
sx q[0];
rz(-1.4365124) q[0];
sx q[0];
rz(-2.4875212) q[0];
rz(-1.7895128) q[1];
sx q[1];
rz(-0.65579954) q[1];
sx q[1];
rz(3.0316614) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39449542) q[0];
sx q[0];
rz(-1.5909068) q[0];
sx q[0];
rz(-3.1294401) q[0];
rz(0.97642297) q[2];
sx q[2];
rz(-2.4196673) q[2];
sx q[2];
rz(-0.41625574) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.47418865) q[1];
sx q[1];
rz(-1.1193573) q[1];
sx q[1];
rz(1.3787446) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5853048) q[3];
sx q[3];
rz(-2.6835052) q[3];
sx q[3];
rz(-2.213495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8249417) q[2];
sx q[2];
rz(-1.9330838) q[2];
sx q[2];
rz(1.4671154) q[2];
rz(-1.0351099) q[3];
sx q[3];
rz(-2.3605774) q[3];
sx q[3];
rz(-1.3181814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.3159981) q[0];
sx q[0];
rz(-0.23402973) q[0];
sx q[0];
rz(-2.3509534) q[0];
rz(2.3979893) q[1];
sx q[1];
rz(-1.5433106) q[1];
sx q[1];
rz(2.5620983) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42946896) q[0];
sx q[0];
rz(-2.3050024) q[0];
sx q[0];
rz(2.3267158) q[0];
rz(-pi) q[1];
rz(-2.8786009) q[2];
sx q[2];
rz(-1.8557669) q[2];
sx q[2];
rz(2.6927352) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5011141) q[1];
sx q[1];
rz(-1.5334629) q[1];
sx q[1];
rz(1.3126612) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73852818) q[3];
sx q[3];
rz(-2.4721842) q[3];
sx q[3];
rz(-2.6311292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6774595) q[2];
sx q[2];
rz(-0.83659283) q[2];
sx q[2];
rz(0.38132384) q[2];
rz(2.2903806) q[3];
sx q[3];
rz(-0.92654735) q[3];
sx q[3];
rz(-2.4066431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6094991) q[0];
sx q[0];
rz(-0.61436009) q[0];
sx q[0];
rz(-0.36439782) q[0];
rz(-0.20026194) q[1];
sx q[1];
rz(-1.7297144) q[1];
sx q[1];
rz(-0.099954896) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6969589) q[0];
sx q[0];
rz(-1.4958188) q[0];
sx q[0];
rz(3.1211389) q[0];
x q[1];
rz(1.7968788) q[2];
sx q[2];
rz(-0.43856171) q[2];
sx q[2];
rz(-1.7600702) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2435088) q[1];
sx q[1];
rz(-1.5931411) q[1];
sx q[1];
rz(3.0490655) q[1];
rz(-pi) q[2];
rz(-0.95619802) q[3];
sx q[3];
rz(-1.0805939) q[3];
sx q[3];
rz(2.4407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0692856) q[2];
sx q[2];
rz(-2.0733209) q[2];
sx q[2];
rz(3.0266673) q[2];
rz(1.140444) q[3];
sx q[3];
rz(-1.2830257) q[3];
sx q[3];
rz(0.97525245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2045778) q[0];
sx q[0];
rz(-2.1890722) q[0];
sx q[0];
rz(0.17247795) q[0];
rz(-1.0109488) q[1];
sx q[1];
rz(-0.5609678) q[1];
sx q[1];
rz(2.9867461) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.895899) q[0];
sx q[0];
rz(-1.2749199) q[0];
sx q[0];
rz(-0.22801836) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2482951) q[2];
sx q[2];
rz(-0.66993827) q[2];
sx q[2];
rz(-2.8084286) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.973792) q[1];
sx q[1];
rz(-1.1607045) q[1];
sx q[1];
rz(-1.4136259) q[1];
rz(-1.6825292) q[3];
sx q[3];
rz(-1.3936685) q[3];
sx q[3];
rz(-2.2551409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.054472063) q[2];
sx q[2];
rz(-0.27721578) q[2];
sx q[2];
rz(-0.37230125) q[2];
rz(-2.7436658) q[3];
sx q[3];
rz(-1.9979265) q[3];
sx q[3];
rz(3.1411689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.953124) q[0];
sx q[0];
rz(-0.57250452) q[0];
sx q[0];
rz(1.5931801) q[0];
rz(0.36987034) q[1];
sx q[1];
rz(-1.7431755) q[1];
sx q[1];
rz(0.63873783) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67154449) q[0];
sx q[0];
rz(-0.28451583) q[0];
sx q[0];
rz(-2.1485062) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9119092) q[2];
sx q[2];
rz(-0.55160597) q[2];
sx q[2];
rz(2.6831339) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7914541) q[1];
sx q[1];
rz(-1.119507) q[1];
sx q[1];
rz(2.7103488) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3494598) q[3];
sx q[3];
rz(-1.9895305) q[3];
sx q[3];
rz(1.8804388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0686191) q[2];
sx q[2];
rz(-2.0899453) q[2];
sx q[2];
rz(0.2529141) q[2];
rz(0.84826338) q[3];
sx q[3];
rz(-1.4784644) q[3];
sx q[3];
rz(-3.0566791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10097583) q[0];
sx q[0];
rz(-0.210013) q[0];
sx q[0];
rz(-0.14895359) q[0];
rz(1.8303998) q[1];
sx q[1];
rz(-1.7313749) q[1];
sx q[1];
rz(3.0874918) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7145007) q[0];
sx q[0];
rz(-2.6339871) q[0];
sx q[0];
rz(2.9604594) q[0];
rz(-0.04345036) q[2];
sx q[2];
rz(-0.80294007) q[2];
sx q[2];
rz(2.174365) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91840345) q[1];
sx q[1];
rz(-1.4610054) q[1];
sx q[1];
rz(-2.4415183) q[1];
rz(-pi) q[2];
rz(1.0955174) q[3];
sx q[3];
rz(-1.5500808) q[3];
sx q[3];
rz(-1.4521862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3397303) q[2];
sx q[2];
rz(-1.170155) q[2];
sx q[2];
rz(2.1745963) q[2];
rz(2.4407834) q[3];
sx q[3];
rz(-0.98027027) q[3];
sx q[3];
rz(2.7907659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12857777) q[0];
sx q[0];
rz(-1.9676493) q[0];
sx q[0];
rz(-1.5482192) q[0];
rz(2.7633527) q[1];
sx q[1];
rz(-2.389237) q[1];
sx q[1];
rz(0.38984782) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5704351) q[0];
sx q[0];
rz(-1.7661372) q[0];
sx q[0];
rz(-0.51049149) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8748449) q[2];
sx q[2];
rz(-1.5549193) q[2];
sx q[2];
rz(-0.35922633) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2193824) q[1];
sx q[1];
rz(-1.2212906) q[1];
sx q[1];
rz(0.52095682) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7933361) q[3];
sx q[3];
rz(-1.3218109) q[3];
sx q[3];
rz(0.96446645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0491911) q[2];
sx q[2];
rz(-1.1213877) q[2];
sx q[2];
rz(-1.8828877) q[2];
rz(-2.8973268) q[3];
sx q[3];
rz(-1.786307) q[3];
sx q[3];
rz(0.99610966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8213537) q[0];
sx q[0];
rz(-1.2293674) q[0];
sx q[0];
rz(1.7852596) q[0];
rz(-2.9755196) q[1];
sx q[1];
rz(-2.1898654) q[1];
sx q[1];
rz(0.43620268) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22101519) q[0];
sx q[0];
rz(-1.2599143) q[0];
sx q[0];
rz(-2.6665844) q[0];
rz(-pi) q[1];
rz(-2.2460552) q[2];
sx q[2];
rz(-2.6288599) q[2];
sx q[2];
rz(-0.89287478) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.055298565) q[1];
sx q[1];
rz(-2.3196967) q[1];
sx q[1];
rz(2.8248252) q[1];
x q[2];
rz(-2.7049191) q[3];
sx q[3];
rz(-1.5092351) q[3];
sx q[3];
rz(-1.7583282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0210375) q[2];
sx q[2];
rz(-1.1641116) q[2];
sx q[2];
rz(0.54350054) q[2];
rz(0.049526878) q[3];
sx q[3];
rz(-1.6560358) q[3];
sx q[3];
rz(-1.4878976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.6108342) q[0];
sx q[0];
rz(-0.13228358) q[0];
sx q[0];
rz(3.0623867) q[0];
rz(-0.40009701) q[1];
sx q[1];
rz(-0.85405093) q[1];
sx q[1];
rz(-1.0265464) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0745747) q[0];
sx q[0];
rz(-1.3368784) q[0];
sx q[0];
rz(-0.83252711) q[0];
x q[1];
rz(-2.7753749) q[2];
sx q[2];
rz(-1.4819659) q[2];
sx q[2];
rz(-0.41761145) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0594008) q[1];
sx q[1];
rz(-0.33723661) q[1];
sx q[1];
rz(2.853501) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9273716) q[3];
sx q[3];
rz(-0.35839265) q[3];
sx q[3];
rz(0.58979366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2058699) q[2];
sx q[2];
rz(-1.8368072) q[2];
sx q[2];
rz(1.4523466) q[2];
rz(2.3116889) q[3];
sx q[3];
rz(-0.87703505) q[3];
sx q[3];
rz(-0.95480603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6185388) q[0];
sx q[0];
rz(-1.9970311) q[0];
sx q[0];
rz(1.507623) q[0];
rz(1.2670831) q[1];
sx q[1];
rz(-2.0590084) q[1];
sx q[1];
rz(0.72437292) q[1];
rz(2.3054068) q[2];
sx q[2];
rz(-2.4124574) q[2];
sx q[2];
rz(-1.3395723) q[2];
rz(-1.2290365) q[3];
sx q[3];
rz(-2.1038341) q[3];
sx q[3];
rz(-2.9084222) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
