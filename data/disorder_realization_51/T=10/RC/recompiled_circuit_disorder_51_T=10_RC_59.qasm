OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5470619) q[0];
sx q[0];
rz(-2.0201595) q[0];
sx q[0];
rz(2.9601331) q[0];
rz(2.060086) q[1];
sx q[1];
rz(-0.67343155) q[1];
sx q[1];
rz(-1.0531309) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5232789) q[0];
sx q[0];
rz(-1.7821454) q[0];
sx q[0];
rz(1.1002512) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5187831) q[2];
sx q[2];
rz(-1.8138759) q[2];
sx q[2];
rz(2.958975) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6932106) q[1];
sx q[1];
rz(-2.3464077) q[1];
sx q[1];
rz(-0.27547142) q[1];
rz(-pi) q[2];
rz(2.9795322) q[3];
sx q[3];
rz(-1.9063623) q[3];
sx q[3];
rz(3.0367362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9782605) q[2];
sx q[2];
rz(-1.1085199) q[2];
sx q[2];
rz(1.367761) q[2];
rz(1.0129499) q[3];
sx q[3];
rz(-0.84665853) q[3];
sx q[3];
rz(-0.071454123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.0435836) q[0];
sx q[0];
rz(-1.4494891) q[0];
sx q[0];
rz(-0.59511551) q[0];
rz(-2.0960506) q[1];
sx q[1];
rz(-1.7273993) q[1];
sx q[1];
rz(-1.5140623) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94905108) q[0];
sx q[0];
rz(-0.66883123) q[0];
sx q[0];
rz(-2.2632954) q[0];
rz(-pi) q[1];
rz(-1.5397649) q[2];
sx q[2];
rz(-2.4040262) q[2];
sx q[2];
rz(-0.34027983) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6914312) q[1];
sx q[1];
rz(-1.3988004) q[1];
sx q[1];
rz(1.3732234) q[1];
rz(-pi) q[2];
rz(-2.0239003) q[3];
sx q[3];
rz(-0.11370224) q[3];
sx q[3];
rz(-0.5773471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65511584) q[2];
sx q[2];
rz(-2.1003508) q[2];
sx q[2];
rz(0.79616037) q[2];
rz(-2.1697309) q[3];
sx q[3];
rz(-2.4367417) q[3];
sx q[3];
rz(2.9698353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7650836) q[0];
sx q[0];
rz(-2.3628545) q[0];
sx q[0];
rz(0.064095108) q[0];
rz(-0.31072101) q[1];
sx q[1];
rz(-1.4704082) q[1];
sx q[1];
rz(1.6832738) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1761988) q[0];
sx q[0];
rz(-1.6178693) q[0];
sx q[0];
rz(1.4233627) q[0];
rz(-pi) q[1];
rz(-2.4086191) q[2];
sx q[2];
rz(-2.3419215) q[2];
sx q[2];
rz(0.016499585) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7160733) q[1];
sx q[1];
rz(-2.7287373) q[1];
sx q[1];
rz(-3.0070261) q[1];
rz(1.5879257) q[3];
sx q[3];
rz(-1.8855842) q[3];
sx q[3];
rz(-1.1312315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98214275) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(-2.4625051) q[2];
rz(1.7689765) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8452334) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(-2.0171719) q[0];
rz(-1.9212978) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(1.4845928) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9254018) q[0];
sx q[0];
rz(-1.5201709) q[0];
sx q[0];
rz(-0.030406818) q[0];
x q[1];
rz(1.053327) q[2];
sx q[2];
rz(-1.7382858) q[2];
sx q[2];
rz(-2.6094764) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.9323737) q[1];
sx q[1];
rz(-2.0062431) q[1];
sx q[1];
rz(-2.2030764) q[1];
x q[2];
rz(1.2403537) q[3];
sx q[3];
rz(-2.2664865) q[3];
sx q[3];
rz(2.2479923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1241887) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(1.0162639) q[2];
rz(1.4034363) q[3];
sx q[3];
rz(-1.4973463) q[3];
sx q[3];
rz(-2.0531634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6338585) q[0];
sx q[0];
rz(-1.2861179) q[0];
sx q[0];
rz(0.39988363) q[0];
rz(1.1625066) q[1];
sx q[1];
rz(-1.3299273) q[1];
sx q[1];
rz(-0.17366017) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6337316) q[0];
sx q[0];
rz(-1.5989688) q[0];
sx q[0];
rz(0.072576056) q[0];
rz(-pi) q[1];
rz(-1.3599672) q[2];
sx q[2];
rz(-1.1846917) q[2];
sx q[2];
rz(2.2112276) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.793321) q[1];
sx q[1];
rz(-1.4655359) q[1];
sx q[1];
rz(-0.028491032) q[1];
rz(0.33186121) q[3];
sx q[3];
rz(-1.3734986) q[3];
sx q[3];
rz(2.8062537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.48012039) q[2];
sx q[2];
rz(-1.1602594) q[2];
sx q[2];
rz(-2.373467) q[2];
rz(2.2875732) q[3];
sx q[3];
rz(-1.7212399) q[3];
sx q[3];
rz(2.1693726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.0356692) q[0];
sx q[0];
rz(-2.8978455) q[0];
sx q[0];
rz(-1.4703898) q[0];
rz(0.51180965) q[1];
sx q[1];
rz(-0.51135951) q[1];
sx q[1];
rz(1.9981729) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9914046) q[0];
sx q[0];
rz(-2.9692869) q[0];
sx q[0];
rz(-2.6323787) q[0];
x q[1];
rz(-1.2406232) q[2];
sx q[2];
rz(-2.1483443) q[2];
sx q[2];
rz(1.4616879) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9547564) q[1];
sx q[1];
rz(-1.3969159) q[1];
sx q[1];
rz(1.9435005) q[1];
x q[2];
rz(2.0538123) q[3];
sx q[3];
rz(-2.5452168) q[3];
sx q[3];
rz(3.0697825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64289552) q[2];
sx q[2];
rz(-0.84474793) q[2];
sx q[2];
rz(1.5303622) q[2];
rz(1.4536084) q[3];
sx q[3];
rz(-2.0724824) q[3];
sx q[3];
rz(2.8924275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11944184) q[0];
sx q[0];
rz(-2.3972153) q[0];
sx q[0];
rz(1.4402333) q[0];
rz(-0.72921905) q[1];
sx q[1];
rz(-1.1518642) q[1];
sx q[1];
rz(-1.1332606) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2204809) q[0];
sx q[0];
rz(-2.8163914) q[0];
sx q[0];
rz(-3.0082506) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1256345) q[2];
sx q[2];
rz(-0.71231132) q[2];
sx q[2];
rz(-0.51770891) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3164191) q[1];
sx q[1];
rz(-2.2146041) q[1];
sx q[1];
rz(2.7979147) q[1];
rz(-pi) q[2];
rz(0.53448581) q[3];
sx q[3];
rz(-0.30189461) q[3];
sx q[3];
rz(-1.2045977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7363654) q[2];
sx q[2];
rz(-2.0183125) q[2];
sx q[2];
rz(0.48842946) q[2];
rz(-1.8296261) q[3];
sx q[3];
rz(-0.043881504) q[3];
sx q[3];
rz(2.5002938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0768123) q[0];
sx q[0];
rz(-1.8490054) q[0];
sx q[0];
rz(0.07117614) q[0];
rz(-3.1094303) q[1];
sx q[1];
rz(-1.3379438) q[1];
sx q[1];
rz(1.2088998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9933388) q[0];
sx q[0];
rz(-0.78972602) q[0];
sx q[0];
rz(1.2475616) q[0];
x q[1];
rz(1.8078631) q[2];
sx q[2];
rz(-1.3915239) q[2];
sx q[2];
rz(-0.2722424) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4603685) q[1];
sx q[1];
rz(-1.7011233) q[1];
sx q[1];
rz(0.016101109) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6505359) q[3];
sx q[3];
rz(-2.4591082) q[3];
sx q[3];
rz(1.3852937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.20748392) q[2];
sx q[2];
rz(-2.948163) q[2];
sx q[2];
rz(1.5709546) q[2];
rz(-2.2682244) q[3];
sx q[3];
rz(-1.4040754) q[3];
sx q[3];
rz(-2.7895555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97380012) q[0];
sx q[0];
rz(-1.6163102) q[0];
sx q[0];
rz(2.8299676) q[0];
rz(2.3198126) q[1];
sx q[1];
rz(-0.59097925) q[1];
sx q[1];
rz(-1.6315546) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8845997) q[0];
sx q[0];
rz(-1.25367) q[0];
sx q[0];
rz(1.3604128) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19599234) q[2];
sx q[2];
rz(-1.2336944) q[2];
sx q[2];
rz(-2.5981284) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3703287) q[1];
sx q[1];
rz(-0.9653829) q[1];
sx q[1];
rz(-1.3243115) q[1];
rz(-pi) q[2];
rz(-0.66442499) q[3];
sx q[3];
rz(-2.8492152) q[3];
sx q[3];
rz(-1.9362621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8682378) q[2];
sx q[2];
rz(-2.6028825) q[2];
sx q[2];
rz(-0.81595016) q[2];
rz(-2.6319035) q[3];
sx q[3];
rz(-0.78521252) q[3];
sx q[3];
rz(1.9036487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2508535) q[0];
sx q[0];
rz(-1.4287404) q[0];
sx q[0];
rz(-2.9113286) q[0];
rz(0.62581217) q[1];
sx q[1];
rz(-0.9451378) q[1];
sx q[1];
rz(0.65840107) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2347682) q[0];
sx q[0];
rz(-2.0567729) q[0];
sx q[0];
rz(-2.4776239) q[0];
rz(-pi) q[1];
rz(2.7294331) q[2];
sx q[2];
rz(-2.4240652) q[2];
sx q[2];
rz(1.5199496) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.876993) q[1];
sx q[1];
rz(-2.0473285) q[1];
sx q[1];
rz(-2.9546253) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53604605) q[3];
sx q[3];
rz(-1.8744161) q[3];
sx q[3];
rz(2.924078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4663503) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(2.8038483) q[2];
rz(-2.1181469) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(1.9256928) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52453775) q[0];
sx q[0];
rz(-0.90159566) q[0];
sx q[0];
rz(3.0774075) q[0];
rz(1.2383923) q[1];
sx q[1];
rz(-1.7108142) q[1];
sx q[1];
rz(-1.6731813) q[1];
rz(2.1970489) q[2];
sx q[2];
rz(-1.7164451) q[2];
sx q[2];
rz(-3.0838983) q[2];
rz(1.6347377) q[3];
sx q[3];
rz(-2.1422374) q[3];
sx q[3];
rz(-0.16850785) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
