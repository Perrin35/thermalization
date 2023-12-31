OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0467779) q[0];
sx q[0];
rz(-1.0682286) q[0];
sx q[0];
rz(2.6775223) q[0];
rz(-1.1820473) q[1];
sx q[1];
rz(-3.074475) q[1];
sx q[1];
rz(-1.2844515) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099079236) q[0];
sx q[0];
rz(-2.6455542) q[0];
sx q[0];
rz(-0.92143671) q[0];
rz(-pi) q[1];
rz(3.0797144) q[2];
sx q[2];
rz(-2.6510694) q[2];
sx q[2];
rz(2.2762736) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.97548188) q[1];
sx q[1];
rz(-1.2088641) q[1];
sx q[1];
rz(-0.48717498) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9908882) q[3];
sx q[3];
rz(-1.5871443) q[3];
sx q[3];
rz(-0.63936641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7444732) q[2];
sx q[2];
rz(-0.93966547) q[2];
sx q[2];
rz(0.31952566) q[2];
rz(0.5784936) q[3];
sx q[3];
rz(-0.47839034) q[3];
sx q[3];
rz(2.4676676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4085061) q[0];
sx q[0];
rz(-1.5698313) q[0];
sx q[0];
rz(-0.52655667) q[0];
rz(0.56354848) q[1];
sx q[1];
rz(-1.5310409) q[1];
sx q[1];
rz(0.79663509) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63882534) q[0];
sx q[0];
rz(-0.728038) q[0];
sx q[0];
rz(1.0685705) q[0];
rz(1.37155) q[2];
sx q[2];
rz(-1.0052048) q[2];
sx q[2];
rz(2.9233962) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.14256515) q[1];
sx q[1];
rz(-1.5836645) q[1];
sx q[1];
rz(2.5343115) q[1];
rz(-pi) q[2];
rz(2.1809686) q[3];
sx q[3];
rz(-1.989813) q[3];
sx q[3];
rz(0.39715365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9600296) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(-0.78655085) q[2];
rz(-0.49318796) q[3];
sx q[3];
rz(-1.1816053) q[3];
sx q[3];
rz(0.44737679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2717473) q[0];
sx q[0];
rz(-2.2646876) q[0];
sx q[0];
rz(-1.7012117) q[0];
rz(2.4213743) q[1];
sx q[1];
rz(-0.59599382) q[1];
sx q[1];
rz(-0.70297855) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1709135) q[0];
sx q[0];
rz(-1.850607) q[0];
sx q[0];
rz(-0.12165102) q[0];
x q[1];
rz(-0.59962745) q[2];
sx q[2];
rz(-2.5587974) q[2];
sx q[2];
rz(-1.8011013) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2762201) q[1];
sx q[1];
rz(-0.11905383) q[1];
sx q[1];
rz(2.7369376) q[1];
rz(-pi) q[2];
rz(-0.17351563) q[3];
sx q[3];
rz(-0.95458889) q[3];
sx q[3];
rz(0.71615744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26677033) q[2];
sx q[2];
rz(-1.6569933) q[2];
sx q[2];
rz(-0.91153574) q[2];
rz(-2.2276145) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(-0.82733697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5660969) q[0];
sx q[0];
rz(-0.63593447) q[0];
sx q[0];
rz(-0.77600586) q[0];
rz(-1.2674747) q[1];
sx q[1];
rz(-1.1088561) q[1];
sx q[1];
rz(-2.5783096) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11382064) q[0];
sx q[0];
rz(-1.8437244) q[0];
sx q[0];
rz(0.91822894) q[0];
rz(-pi) q[1];
rz(1.9070542) q[2];
sx q[2];
rz(-0.90869892) q[2];
sx q[2];
rz(-1.3627571) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5717585) q[1];
sx q[1];
rz(-0.27184871) q[1];
sx q[1];
rz(3.1175201) q[1];
rz(-pi) q[2];
rz(2.8068845) q[3];
sx q[3];
rz(-2.3957806) q[3];
sx q[3];
rz(1.0583744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6528066) q[2];
sx q[2];
rz(-0.93549171) q[2];
sx q[2];
rz(2.9768067) q[2];
rz(2.9131043) q[3];
sx q[3];
rz(-2.861664) q[3];
sx q[3];
rz(-0.94648186) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8673458) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(2.2891323) q[0];
rz(0.35119855) q[1];
sx q[1];
rz(-2.5636667) q[1];
sx q[1];
rz(1.9794827) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64415414) q[0];
sx q[0];
rz(-2.4442406) q[0];
sx q[0];
rz(-2.5028412) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2010872) q[2];
sx q[2];
rz(-2.5224707) q[2];
sx q[2];
rz(-1.8096015) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3106766) q[1];
sx q[1];
rz(-1.984593) q[1];
sx q[1];
rz(2.3526741) q[1];
rz(2.4200053) q[3];
sx q[3];
rz(-1.5820832) q[3];
sx q[3];
rz(-1.23502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.44624415) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(1.2472786) q[2];
rz(-0.013109664) q[3];
sx q[3];
rz(-1.40991) q[3];
sx q[3];
rz(2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2728249) q[0];
sx q[0];
rz(-1.2478963) q[0];
sx q[0];
rz(-0.75063467) q[0];
rz(1.1095095) q[1];
sx q[1];
rz(-0.37056071) q[1];
sx q[1];
rz(1.3060588) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40697843) q[0];
sx q[0];
rz(-2.4771871) q[0];
sx q[0];
rz(-2.2218496) q[0];
rz(-0.38988955) q[2];
sx q[2];
rz(-2.6460558) q[2];
sx q[2];
rz(1.3930266) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8522779) q[1];
sx q[1];
rz(-1.5631952) q[1];
sx q[1];
rz(-0.72035933) q[1];
rz(-pi) q[2];
rz(2.552794) q[3];
sx q[3];
rz(-1.4466803) q[3];
sx q[3];
rz(-2.0037946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7541472) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(-2.5409017) q[2];
rz(1.0026275) q[3];
sx q[3];
rz(-0.53835821) q[3];
sx q[3];
rz(2.7887662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3951185) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(-1.0466928) q[0];
rz(-1.612161) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(0.41710645) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5865267) q[0];
sx q[0];
rz(-1.546372) q[0];
sx q[0];
rz(1.3506372) q[0];
rz(-pi) q[1];
rz(0.76191683) q[2];
sx q[2];
rz(-0.41315213) q[2];
sx q[2];
rz(-2.6921536) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4579791) q[1];
sx q[1];
rz(-0.90404592) q[1];
sx q[1];
rz(2.1702106) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8700637) q[3];
sx q[3];
rz(-2.2021658) q[3];
sx q[3];
rz(2.805998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.0059011857) q[2];
sx q[2];
rz(-2.6689172) q[2];
sx q[2];
rz(-1.3674412) q[2];
rz(-0.55316365) q[3];
sx q[3];
rz(-1.7382712) q[3];
sx q[3];
rz(2.6161391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.816514) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(-1.9518071) q[0];
rz(-1.7143543) q[1];
sx q[1];
rz(-0.27806565) q[1];
sx q[1];
rz(-0.11238012) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053514078) q[0];
sx q[0];
rz(-1.3980165) q[0];
sx q[0];
rz(1.0900351) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24765315) q[2];
sx q[2];
rz(-0.34341771) q[2];
sx q[2];
rz(-1.0719971) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8466859) q[1];
sx q[1];
rz(-1.2290188) q[1];
sx q[1];
rz(2.3218367) q[1];
rz(-0.12179575) q[3];
sx q[3];
rz(-1.8276916) q[3];
sx q[3];
rz(1.6459873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0671976) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(-1.3486264) q[2];
rz(-1.2049234) q[3];
sx q[3];
rz(-1.7347074) q[3];
sx q[3];
rz(-0.19781923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39712054) q[0];
sx q[0];
rz(-1.67698) q[0];
sx q[0];
rz(0.034974139) q[0];
rz(2.2947625) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(2.2299178) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.864894) q[0];
sx q[0];
rz(-1.3008586) q[0];
sx q[0];
rz(0.2322659) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5981204) q[2];
sx q[2];
rz(-1.6795571) q[2];
sx q[2];
rz(0.47765884) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8988077) q[1];
sx q[1];
rz(-2.6344732) q[1];
sx q[1];
rz(-2.7125263) q[1];
rz(-pi) q[2];
rz(-0.53578844) q[3];
sx q[3];
rz(-2.0013323) q[3];
sx q[3];
rz(-2.146194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.70790616) q[2];
sx q[2];
rz(-2.6277379) q[2];
sx q[2];
rz(-0.24924499) q[2];
rz(-2.3748659) q[3];
sx q[3];
rz(-1.3336351) q[3];
sx q[3];
rz(2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7837759) q[0];
sx q[0];
rz(-2.0932842) q[0];
sx q[0];
rz(-0.37208474) q[0];
rz(2.5601939) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(1.4153597) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85522643) q[0];
sx q[0];
rz(-0.33272538) q[0];
sx q[0];
rz(-0.65441982) q[0];
rz(2.084311) q[2];
sx q[2];
rz(-1.0143177) q[2];
sx q[2];
rz(1.2531812) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.49299875) q[1];
sx q[1];
rz(-3.0393638) q[1];
sx q[1];
rz(-2.3962254) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1595702) q[3];
sx q[3];
rz(-0.62871274) q[3];
sx q[3];
rz(-0.94129896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8156478) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(-0.20467219) q[2];
rz(1.4137911) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(-1.0958825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50080147) q[0];
sx q[0];
rz(-1.4688107) q[0];
sx q[0];
rz(-2.2039913) q[0];
rz(-1.5564556) q[1];
sx q[1];
rz(-1.7592447) q[1];
sx q[1];
rz(-1.6978282) q[1];
rz(1.8150868) q[2];
sx q[2];
rz(-0.46637022) q[2];
sx q[2];
rz(-2.7844219) q[2];
rz(2.1429569) q[3];
sx q[3];
rz(-1.5012267) q[3];
sx q[3];
rz(-2.5517626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
