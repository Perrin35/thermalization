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
rz(-0.7260538) q[0];
rz(0.59745204) q[1];
sx q[1];
rz(-2.6260881) q[1];
sx q[1];
rz(-2.1093624) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2054133) q[0];
sx q[0];
rz(-1.4698615) q[0];
sx q[0];
rz(-2.4839782) q[0];
rz(-1.6392567) q[2];
sx q[2];
rz(-0.51883829) q[2];
sx q[2];
rz(-1.8585132) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1223678) q[1];
sx q[1];
rz(-1.4488954) q[1];
sx q[1];
rz(-0.31369622) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0606403) q[3];
sx q[3];
rz(-2.3054696) q[3];
sx q[3];
rz(-1.6250798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.44077474) q[2];
sx q[2];
rz(-2.5371607) q[2];
sx q[2];
rz(3.0536998) q[2];
rz(1.4676189) q[3];
sx q[3];
rz(-1.2940977) q[3];
sx q[3];
rz(-0.99883336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0198332) q[0];
sx q[0];
rz(-1.8353945) q[0];
sx q[0];
rz(1.649296) q[0];
rz(0.78798931) q[1];
sx q[1];
rz(-1.183527) q[1];
sx q[1];
rz(-0.96510395) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0168646) q[0];
sx q[0];
rz(-1.6627321) q[0];
sx q[0];
rz(-2.4107736) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4282706) q[2];
sx q[2];
rz(-2.3025142) q[2];
sx q[2];
rz(-2.1676262) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1652884) q[1];
sx q[1];
rz(-1.3852264) q[1];
sx q[1];
rz(1.1052422) q[1];
x q[2];
rz(1.4802259) q[3];
sx q[3];
rz(-2.0388868) q[3];
sx q[3];
rz(0.15453574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1284156) q[2];
sx q[2];
rz(-2.8030277) q[2];
sx q[2];
rz(1.817912) q[2];
rz(2.2549818) q[3];
sx q[3];
rz(-1.9804136) q[3];
sx q[3];
rz(-2.9429341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.897268) q[0];
sx q[0];
rz(-2.1710945) q[0];
sx q[0];
rz(-2.4460728) q[0];
rz(0.8356525) q[1];
sx q[1];
rz(-1.7213768) q[1];
sx q[1];
rz(-2.4998891) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4280689) q[0];
sx q[0];
rz(-1.4260573) q[0];
sx q[0];
rz(-0.83065303) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19410816) q[2];
sx q[2];
rz(-2.420937) q[2];
sx q[2];
rz(-1.0247165) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56395951) q[1];
sx q[1];
rz(-2.0490987) q[1];
sx q[1];
rz(-0.35525124) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66777457) q[3];
sx q[3];
rz(-2.3333356) q[3];
sx q[3];
rz(-2.2796352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.16331638) q[2];
sx q[2];
rz(-2.0127608) q[2];
sx q[2];
rz(-1.8070492) q[2];
rz(-1.4000019) q[3];
sx q[3];
rz(-1.497523) q[3];
sx q[3];
rz(-2.002772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9982933) q[0];
sx q[0];
rz(-0.87109733) q[0];
sx q[0];
rz(2.368108) q[0];
rz(-0.086291226) q[1];
sx q[1];
rz(-2.6373865) q[1];
sx q[1];
rz(-0.48826826) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7542758) q[0];
sx q[0];
rz(-1.6389221) q[0];
sx q[0];
rz(-0.46275567) q[0];
x q[1];
rz(2.1147229) q[2];
sx q[2];
rz(-0.94878093) q[2];
sx q[2];
rz(-1.1116456) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.5544071) q[1];
sx q[1];
rz(-2.3766365) q[1];
sx q[1];
rz(-1.2423963) q[1];
x q[2];
rz(1.9776439) q[3];
sx q[3];
rz(-1.9449807) q[3];
sx q[3];
rz(2.0896623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.46828541) q[2];
sx q[2];
rz(-1.1136592) q[2];
sx q[2];
rz(-0.16806531) q[2];
rz(-0.42260653) q[3];
sx q[3];
rz(-0.97413078) q[3];
sx q[3];
rz(-2.9882123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3733658) q[0];
sx q[0];
rz(-2.234937) q[0];
sx q[0];
rz(1.1905097) q[0];
rz(-0.52781421) q[1];
sx q[1];
rz(-2.0595136) q[1];
sx q[1];
rz(-3.1324918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4531698) q[0];
sx q[0];
rz(-0.97081447) q[0];
sx q[0];
rz(-1.021941) q[0];
x q[1];
rz(-1.6318237) q[2];
sx q[2];
rz(-1.5102855) q[2];
sx q[2];
rz(2.4946545) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.121351) q[1];
sx q[1];
rz(-0.96881908) q[1];
sx q[1];
rz(-2.6247128) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1808257) q[3];
sx q[3];
rz(-2.4159263) q[3];
sx q[3];
rz(-0.31455597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8268383) q[2];
sx q[2];
rz(-0.49586168) q[2];
sx q[2];
rz(-2.4763988) q[2];
rz(-1.0264171) q[3];
sx q[3];
rz(-1.0746936) q[3];
sx q[3];
rz(-0.78072602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6219532) q[0];
sx q[0];
rz(-2.6060947) q[0];
sx q[0];
rz(-2.7435379) q[0];
rz(-1.4705426) q[1];
sx q[1];
rz(-2.2649951) q[1];
sx q[1];
rz(2.3614531) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3624293) q[0];
sx q[0];
rz(-2.0423959) q[0];
sx q[0];
rz(0.40059025) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4498135) q[2];
sx q[2];
rz(-2.8304812) q[2];
sx q[2];
rz(1.9115314) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.33360219) q[1];
sx q[1];
rz(-1.8544852) q[1];
sx q[1];
rz(-2.4391443) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8858812) q[3];
sx q[3];
rz(-1.3217333) q[3];
sx q[3];
rz(-2.0307756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2842497) q[2];
sx q[2];
rz(-1.1025068) q[2];
sx q[2];
rz(2.9798689) q[2];
rz(0.11180793) q[3];
sx q[3];
rz(-2.4155152) q[3];
sx q[3];
rz(0.73981729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9991456) q[0];
sx q[0];
rz(-1.7113926) q[0];
sx q[0];
rz(2.4430742) q[0];
rz(2.0493719) q[1];
sx q[1];
rz(-2.3869546) q[1];
sx q[1];
rz(3.0879367) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.221401) q[0];
sx q[0];
rz(-2.2760253) q[0];
sx q[0];
rz(2.7979275) q[0];
rz(-pi) q[1];
rz(-2.2767792) q[2];
sx q[2];
rz(-2.5891782) q[2];
sx q[2];
rz(0.68589003) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.084890993) q[1];
sx q[1];
rz(-0.47353256) q[1];
sx q[1];
rz(-3.0242306) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14999203) q[3];
sx q[3];
rz(-2.2120471) q[3];
sx q[3];
rz(2.7771726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1412389) q[2];
sx q[2];
rz(-0.69811368) q[2];
sx q[2];
rz(-0.21256438) q[2];
rz(-2.8138748) q[3];
sx q[3];
rz(-0.98723427) q[3];
sx q[3];
rz(1.7879965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89606365) q[0];
sx q[0];
rz(-1.117417) q[0];
sx q[0];
rz(-2.8118706) q[0];
rz(-2.3700736) q[1];
sx q[1];
rz(-0.78116575) q[1];
sx q[1];
rz(0.82569295) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3701009) q[0];
sx q[0];
rz(-2.2326062) q[0];
sx q[0];
rz(1.6812117) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0631376) q[2];
sx q[2];
rz(-0.9758853) q[2];
sx q[2];
rz(1.1537781) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6439425) q[1];
sx q[1];
rz(-1.4621328) q[1];
sx q[1];
rz(0.76246271) q[1];
rz(1.7686033) q[3];
sx q[3];
rz(-1.8231492) q[3];
sx q[3];
rz(-0.41393241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9416435) q[2];
sx q[2];
rz(-1.7295094) q[2];
sx q[2];
rz(-1.6597718) q[2];
rz(-0.09305772) q[3];
sx q[3];
rz(-1.8424282) q[3];
sx q[3];
rz(-0.53946462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.900443) q[0];
sx q[0];
rz(-2.7261782) q[0];
sx q[0];
rz(0.39733091) q[0];
rz(2.1516402) q[1];
sx q[1];
rz(-1.3497458) q[1];
sx q[1];
rz(2.1053402) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62972126) q[0];
sx q[0];
rz(-2.7366834) q[0];
sx q[0];
rz(3.1121503) q[0];
rz(-pi) q[1];
rz(2.7781899) q[2];
sx q[2];
rz(-0.37523233) q[2];
sx q[2];
rz(-2.4979532) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81779382) q[1];
sx q[1];
rz(-1.3518055) q[1];
sx q[1];
rz(-1.4794255) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96817044) q[3];
sx q[3];
rz(-0.60634469) q[3];
sx q[3];
rz(1.9023638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1194666) q[2];
sx q[2];
rz(-0.72212044) q[2];
sx q[2];
rz(1.3163346) q[2];
rz(0.25257603) q[3];
sx q[3];
rz(-1.3467237) q[3];
sx q[3];
rz(-2.778497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(3.0040358) q[0];
sx q[0];
rz(-2.3513849) q[0];
sx q[0];
rz(0.23605119) q[0];
rz(-3.0421742) q[1];
sx q[1];
rz(-2.3853018) q[1];
sx q[1];
rz(1.8831467) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61533046) q[0];
sx q[0];
rz(-1.0122504) q[0];
sx q[0];
rz(1.6625151) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0706606) q[2];
sx q[2];
rz(-2.1720083) q[2];
sx q[2];
rz(-2.380033) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7459877) q[1];
sx q[1];
rz(-1.6853991) q[1];
sx q[1];
rz(0.68273165) q[1];
x q[2];
rz(-1.997273) q[3];
sx q[3];
rz(-1.9460554) q[3];
sx q[3];
rz(-2.8464908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1666169) q[2];
sx q[2];
rz(-2.6585343) q[2];
sx q[2];
rz(-1.0413337) q[2];
rz(-2.5552022) q[3];
sx q[3];
rz(-2.6643463) q[3];
sx q[3];
rz(-0.81048107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36622421) q[0];
sx q[0];
rz(-1.0564221) q[0];
sx q[0];
rz(-1.139241) q[0];
rz(2.1495023) q[1];
sx q[1];
rz(-1.5403668) q[1];
sx q[1];
rz(-1.5425727) q[1];
rz(0.61983776) q[2];
sx q[2];
rz(-0.38417338) q[2];
sx q[2];
rz(-1.2585121) q[2];
rz(-2.5358148) q[3];
sx q[3];
rz(-0.68690261) q[3];
sx q[3];
rz(1.5920832) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
