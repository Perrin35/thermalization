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
rz(-1.0567226) q[0];
sx q[0];
rz(-1.280008) q[0];
sx q[0];
rz(2.4155389) q[0];
rz(-2.5441406) q[1];
sx q[1];
rz(-0.51550454) q[1];
sx q[1];
rz(2.1093624) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3772349) q[0];
sx q[0];
rz(-0.66417686) q[0];
sx q[0];
rz(2.9773877) q[0];
rz(-3.1025508) q[2];
sx q[2];
rz(-1.0532951) q[2];
sx q[2];
rz(-1.7797178) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1223678) q[1];
sx q[1];
rz(-1.4488954) q[1];
sx q[1];
rz(0.31369622) q[1];
rz(-pi) q[2];
rz(-1.0606403) q[3];
sx q[3];
rz(-0.83612305) q[3];
sx q[3];
rz(-1.5165129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7008179) q[2];
sx q[2];
rz(-0.60443193) q[2];
sx q[2];
rz(-3.0536998) q[2];
rz(1.6739738) q[3];
sx q[3];
rz(-1.847495) q[3];
sx q[3];
rz(2.1427593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1217594) q[0];
sx q[0];
rz(-1.8353945) q[0];
sx q[0];
rz(-1.4922967) q[0];
rz(2.3536033) q[1];
sx q[1];
rz(-1.183527) q[1];
sx q[1];
rz(0.96510395) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0168646) q[0];
sx q[0];
rz(-1.6627321) q[0];
sx q[0];
rz(-0.73081907) q[0];
rz(-2.9847174) q[2];
sx q[2];
rz(-0.74293908) q[2];
sx q[2];
rz(2.3791831) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9763042) q[1];
sx q[1];
rz(-1.7563662) q[1];
sx q[1];
rz(-2.0363505) q[1];
rz(-0.4697462) q[3];
sx q[3];
rz(-1.6516017) q[3];
sx q[3];
rz(1.4572136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1284156) q[2];
sx q[2];
rz(-0.33856496) q[2];
sx q[2];
rz(1.817912) q[2];
rz(-2.2549818) q[3];
sx q[3];
rz(-1.9804136) q[3];
sx q[3];
rz(-0.19865856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.897268) q[0];
sx q[0];
rz(-0.97049814) q[0];
sx q[0];
rz(0.6955198) q[0];
rz(2.3059402) q[1];
sx q[1];
rz(-1.4202159) q[1];
sx q[1];
rz(0.64170352) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4280689) q[0];
sx q[0];
rz(-1.4260573) q[0];
sx q[0];
rz(2.3109396) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19410816) q[2];
sx q[2];
rz(-0.72065565) q[2];
sx q[2];
rz(-1.0247165) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56395951) q[1];
sx q[1];
rz(-2.0490987) q[1];
sx q[1];
rz(-0.35525124) q[1];
rz(-pi) q[2];
rz(2.4536127) q[3];
sx q[3];
rz(-1.1065346) q[3];
sx q[3];
rz(-1.9339428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9782763) q[2];
sx q[2];
rz(-2.0127608) q[2];
sx q[2];
rz(-1.8070492) q[2];
rz(-1.7415907) q[3];
sx q[3];
rz(-1.6440697) q[3];
sx q[3];
rz(1.1388206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1432994) q[0];
sx q[0];
rz(-0.87109733) q[0];
sx q[0];
rz(-2.368108) q[0];
rz(-3.0553014) q[1];
sx q[1];
rz(-2.6373865) q[1];
sx q[1];
rz(0.48826826) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7542758) q[0];
sx q[0];
rz(-1.5026706) q[0];
sx q[0];
rz(-2.678837) q[0];
rz(0.62520844) q[2];
sx q[2];
rz(-0.80179384) q[2];
sx q[2];
rz(1.9156574) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.99566702) q[1];
sx q[1];
rz(-2.2856376) q[1];
sx q[1];
rz(0.30024528) q[1];
rz(-0.78924307) q[3];
sx q[3];
rz(-2.5960138) q[3];
sx q[3];
rz(-0.18462791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6733072) q[2];
sx q[2];
rz(-2.0279334) q[2];
sx q[2];
rz(-2.9735273) q[2];
rz(2.7189861) q[3];
sx q[3];
rz(-2.1674619) q[3];
sx q[3];
rz(-0.15338038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3733658) q[0];
sx q[0];
rz(-0.90665561) q[0];
sx q[0];
rz(-1.9510829) q[0];
rz(0.52781421) q[1];
sx q[1];
rz(-1.0820791) q[1];
sx q[1];
rz(-3.1324918) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6884228) q[0];
sx q[0];
rz(-2.1707782) q[0];
sx q[0];
rz(1.021941) q[0];
x q[1];
rz(-0.060623364) q[2];
sx q[2];
rz(-1.5098808) q[2];
sx q[2];
rz(-0.9275533) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2796265) q[1];
sx q[1];
rz(-1.1512966) q[1];
sx q[1];
rz(-0.90200938) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6713773) q[3];
sx q[3];
rz(-2.145916) q[3];
sx q[3];
rz(-2.0755538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3147543) q[2];
sx q[2];
rz(-0.49586168) q[2];
sx q[2];
rz(-2.4763988) q[2];
rz(-2.1151755) q[3];
sx q[3];
rz(-1.0746936) q[3];
sx q[3];
rz(-2.3608666) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6219532) q[0];
sx q[0];
rz(-0.53549796) q[0];
sx q[0];
rz(2.7435379) q[0];
rz(-1.6710501) q[1];
sx q[1];
rz(-2.2649951) q[1];
sx q[1];
rz(-2.3614531) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3624293) q[0];
sx q[0];
rz(-2.0423959) q[0];
sx q[0];
rz(-2.7410024) q[0];
rz(-0.69177912) q[2];
sx q[2];
rz(-2.8304812) q[2];
sx q[2];
rz(-1.2300613) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8079905) q[1];
sx q[1];
rz(-1.2871075) q[1];
sx q[1];
rz(2.4391443) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8858812) q[3];
sx q[3];
rz(-1.8198593) q[3];
sx q[3];
rz(2.0307756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.857343) q[2];
sx q[2];
rz(-2.0390859) q[2];
sx q[2];
rz(-0.16172376) q[2];
rz(3.0297847) q[3];
sx q[3];
rz(-0.72607741) q[3];
sx q[3];
rz(0.73981729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9991456) q[0];
sx q[0];
rz(-1.4302) q[0];
sx q[0];
rz(-2.4430742) q[0];
rz(-1.0922208) q[1];
sx q[1];
rz(-0.75463808) q[1];
sx q[1];
rz(-3.0879367) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2630477) q[0];
sx q[0];
rz(-1.8302655) q[0];
sx q[0];
rz(2.3058913) q[0];
rz(2.0094107) q[2];
sx q[2];
rz(-1.918186) q[2];
sx q[2];
rz(1.5127104) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3813579) q[1];
sx q[1];
rz(-1.5173727) q[1];
sx q[1];
rz(-2.670856) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92412432) q[3];
sx q[3];
rz(-1.69083) q[3];
sx q[3];
rz(1.2965352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1412389) q[2];
sx q[2];
rz(-2.443479) q[2];
sx q[2];
rz(2.9290283) q[2];
rz(-2.8138748) q[3];
sx q[3];
rz(-2.1543584) q[3];
sx q[3];
rz(1.3535961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89606365) q[0];
sx q[0];
rz(-2.0241757) q[0];
sx q[0];
rz(-2.8118706) q[0];
rz(-2.3700736) q[1];
sx q[1];
rz(-2.3604269) q[1];
sx q[1];
rz(2.3158997) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86733183) q[0];
sx q[0];
rz(-1.4837588) q[0];
sx q[0];
rz(2.4768193) q[0];
x q[1];
rz(0.97445455) q[2];
sx q[2];
rz(-1.505841) q[2];
sx q[2];
rz(0.37298733) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.030071478) q[1];
sx q[1];
rz(-2.3276405) q[1];
sx q[1];
rz(1.7205283) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8844374) q[3];
sx q[3];
rz(-1.7622593) q[3];
sx q[3];
rz(-1.2068656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9416435) q[2];
sx q[2];
rz(-1.4120833) q[2];
sx q[2];
rz(1.4818209) q[2];
rz(0.09305772) q[3];
sx q[3];
rz(-1.8424282) q[3];
sx q[3];
rz(0.53946462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.900443) q[0];
sx q[0];
rz(-0.41541442) q[0];
sx q[0];
rz(0.39733091) q[0];
rz(-0.98995248) q[1];
sx q[1];
rz(-1.7918469) q[1];
sx q[1];
rz(1.0362524) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91401228) q[0];
sx q[0];
rz(-1.5591994) q[0];
sx q[0];
rz(-0.40475233) q[0];
x q[1];
rz(-0.35276995) q[2];
sx q[2];
rz(-1.4401541) q[2];
sx q[2];
rz(2.5545211) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.73309988) q[1];
sx q[1];
rz(-1.4816135) q[1];
sx q[1];
rz(-0.21987889) q[1];
rz(-pi) q[2];
rz(0.37452368) q[3];
sx q[3];
rz(-2.0595042) q[3];
sx q[3];
rz(-2.599409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1194666) q[2];
sx q[2];
rz(-0.72212044) q[2];
sx q[2];
rz(1.3163346) q[2];
rz(-0.25257603) q[3];
sx q[3];
rz(-1.3467237) q[3];
sx q[3];
rz(-0.36309567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13755688) q[0];
sx q[0];
rz(-2.3513849) q[0];
sx q[0];
rz(-2.9055415) q[0];
rz(-0.099418489) q[1];
sx q[1];
rz(-0.75629083) q[1];
sx q[1];
rz(-1.258446) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1374219) q[0];
sx q[0];
rz(-1.4930471) q[0];
sx q[0];
rz(2.581152) q[0];
rz(0.66345556) q[2];
sx q[2];
rz(-1.1642312) q[2];
sx q[2];
rz(2.6321509) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7459877) q[1];
sx q[1];
rz(-1.4561936) q[1];
sx q[1];
rz(0.68273165) q[1];
rz(-1.1443196) q[3];
sx q[3];
rz(-1.1955373) q[3];
sx q[3];
rz(0.29510185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1666169) q[2];
sx q[2];
rz(-0.48305837) q[2];
sx q[2];
rz(2.100259) q[2];
rz(0.58639041) q[3];
sx q[3];
rz(-0.47724637) q[3];
sx q[3];
rz(0.81048107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7753684) q[0];
sx q[0];
rz(-1.0564221) q[0];
sx q[0];
rz(-1.139241) q[0];
rz(0.99209039) q[1];
sx q[1];
rz(-1.6012259) q[1];
sx q[1];
rz(1.59902) q[1];
rz(-2.5217549) q[2];
sx q[2];
rz(-0.38417338) q[2];
sx q[2];
rz(-1.2585121) q[2];
rz(-2.0076892) q[3];
sx q[3];
rz(-1.0224167) q[3];
sx q[3];
rz(-2.2800048) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
