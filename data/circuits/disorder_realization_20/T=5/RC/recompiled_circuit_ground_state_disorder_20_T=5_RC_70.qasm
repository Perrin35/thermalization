OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5372758) q[0];
sx q[0];
rz(-0.24157) q[0];
sx q[0];
rz(-2.8085652) q[0];
rz(2.0060519) q[1];
sx q[1];
rz(-0.82692868) q[1];
sx q[1];
rz(-0.64396042) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6493131) q[0];
sx q[0];
rz(-2.4631073) q[0];
sx q[0];
rz(-1.1301103) q[0];
rz(-pi) q[1];
rz(-2.8415235) q[2];
sx q[2];
rz(-1.4208394) q[2];
sx q[2];
rz(1.7008925) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.17929303) q[1];
sx q[1];
rz(-0.4518309) q[1];
sx q[1];
rz(2.5746114) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7936739) q[3];
sx q[3];
rz(-1.9609465) q[3];
sx q[3];
rz(-1.8371234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.48646271) q[2];
sx q[2];
rz(-2.6772406) q[2];
sx q[2];
rz(0.74074024) q[2];
rz(1.5898534) q[3];
sx q[3];
rz(-2.430075) q[3];
sx q[3];
rz(0.5927425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084259) q[0];
sx q[0];
rz(-2.0080703) q[0];
sx q[0];
rz(-0.11696996) q[0];
rz(-2.6372657) q[1];
sx q[1];
rz(-1.802899) q[1];
sx q[1];
rz(2.8574944) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11782538) q[0];
sx q[0];
rz(-1.6540915) q[0];
sx q[0];
rz(-1.0040332) q[0];
rz(-0.57055497) q[2];
sx q[2];
rz(-0.88316702) q[2];
sx q[2];
rz(2.2372383) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2843079) q[1];
sx q[1];
rz(-1.6928288) q[1];
sx q[1];
rz(3.1067501) q[1];
rz(-1.8485214) q[3];
sx q[3];
rz(-1.3284995) q[3];
sx q[3];
rz(1.5838069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8154907) q[2];
sx q[2];
rz(-2.2559866) q[2];
sx q[2];
rz(-2.3808114) q[2];
rz(0.86959362) q[3];
sx q[3];
rz(-1.875501) q[3];
sx q[3];
rz(-2.6884955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78075439) q[0];
sx q[0];
rz(-1.1203082) q[0];
sx q[0];
rz(-0.25217062) q[0];
rz(-1.699327) q[1];
sx q[1];
rz(-2.1863054) q[1];
sx q[1];
rz(2.7853277) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2748791) q[0];
sx q[0];
rz(-1.5646114) q[0];
sx q[0];
rz(1.5112108) q[0];
rz(-0.54889955) q[2];
sx q[2];
rz(-1.9268914) q[2];
sx q[2];
rz(2.9402972) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9911451) q[1];
sx q[1];
rz(-2.4887062) q[1];
sx q[1];
rz(-1.8064524) q[1];
rz(-pi) q[2];
rz(-3.1396418) q[3];
sx q[3];
rz(-1.8271128) q[3];
sx q[3];
rz(2.0511829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0226655) q[2];
sx q[2];
rz(-1.7513195) q[2];
sx q[2];
rz(-1.6290132) q[2];
rz(-0.0096983612) q[3];
sx q[3];
rz(-1.9111948) q[3];
sx q[3];
rz(-2.5201918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4951204) q[0];
sx q[0];
rz(-0.54238129) q[0];
sx q[0];
rz(-2.8908308) q[0];
rz(1.7950902) q[1];
sx q[1];
rz(-0.52918068) q[1];
sx q[1];
rz(-1.4432602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0899857) q[0];
sx q[0];
rz(-2.3653226) q[0];
sx q[0];
rz(-1.2597741) q[0];
rz(-pi) q[1];
rz(0.40374229) q[2];
sx q[2];
rz(-1.0761257) q[2];
sx q[2];
rz(-2.1658989) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6175198) q[1];
sx q[1];
rz(-2.1012602) q[1];
sx q[1];
rz(-0.73442119) q[1];
rz(-pi) q[2];
rz(0.26224995) q[3];
sx q[3];
rz(-1.3104386) q[3];
sx q[3];
rz(2.8182056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.5556339) q[2];
sx q[2];
rz(-2.0363753) q[2];
sx q[2];
rz(-0.24331681) q[2];
rz(2.5569052) q[3];
sx q[3];
rz(-0.57716113) q[3];
sx q[3];
rz(3.1134591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1481767) q[0];
sx q[0];
rz(-2.7758444) q[0];
sx q[0];
rz(-2.962033) q[0];
rz(1.9163632) q[1];
sx q[1];
rz(-1.4045249) q[1];
sx q[1];
rz(-2.6833351) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97044635) q[0];
sx q[0];
rz(-1.2752646) q[0];
sx q[0];
rz(1.4895205) q[0];
x q[1];
rz(-0.29693895) q[2];
sx q[2];
rz(-0.7504979) q[2];
sx q[2];
rz(-2.3194882) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1128487) q[1];
sx q[1];
rz(-2.9874871) q[1];
sx q[1];
rz(2.0779559) q[1];
rz(-pi) q[2];
rz(-0.86991258) q[3];
sx q[3];
rz(-1.2906089) q[3];
sx q[3];
rz(2.3126569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7563584) q[2];
sx q[2];
rz(-2.7172654) q[2];
sx q[2];
rz(-2.2198086) q[2];
rz(-1.8900169) q[3];
sx q[3];
rz(-1.7332964) q[3];
sx q[3];
rz(0.061669953) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.09403041) q[0];
sx q[0];
rz(-0.91218364) q[0];
sx q[0];
rz(-2.9929274) q[0];
rz(2.8246763) q[1];
sx q[1];
rz(-0.47873679) q[1];
sx q[1];
rz(2.5247578) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5182523) q[0];
sx q[0];
rz(-3.0419311) q[0];
sx q[0];
rz(2.2355493) q[0];
rz(-pi) q[1];
rz(-0.3601893) q[2];
sx q[2];
rz(-1.9860886) q[2];
sx q[2];
rz(2.5915938) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4839263) q[1];
sx q[1];
rz(-1.0063419) q[1];
sx q[1];
rz(0.29740833) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0471763) q[3];
sx q[3];
rz(-1.3304119) q[3];
sx q[3];
rz(0.8951069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2769015) q[2];
sx q[2];
rz(-1.7128877) q[2];
sx q[2];
rz(-3.1332664) q[2];
rz(0.62638038) q[3];
sx q[3];
rz(-2.4255987) q[3];
sx q[3];
rz(-3.1410419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.83547) q[0];
sx q[0];
rz(-1.1518814) q[0];
sx q[0];
rz(-0.48208153) q[0];
rz(0.029190633) q[1];
sx q[1];
rz(-1.8455285) q[1];
sx q[1];
rz(-0.76404244) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4493443) q[0];
sx q[0];
rz(-1.6381761) q[0];
sx q[0];
rz(-1.199556) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2305626) q[2];
sx q[2];
rz(-1.5507959) q[2];
sx q[2];
rz(-0.80772142) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.205207) q[1];
sx q[1];
rz(-2.1644785) q[1];
sx q[1];
rz(-0.20505814) q[1];
rz(-pi) q[2];
rz(-0.0079844012) q[3];
sx q[3];
rz(-1.1698876) q[3];
sx q[3];
rz(-2.8454091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46334106) q[2];
sx q[2];
rz(-1.8222787) q[2];
sx q[2];
rz(-2.6574262) q[2];
rz(-0.68228996) q[3];
sx q[3];
rz(-0.65410084) q[3];
sx q[3];
rz(0.13474034) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084455647) q[0];
sx q[0];
rz(-0.77780044) q[0];
sx q[0];
rz(-1.8498259) q[0];
rz(0.46503398) q[1];
sx q[1];
rz(-2.6227622) q[1];
sx q[1];
rz(-0.30716392) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7114176) q[0];
sx q[0];
rz(-1.9924376) q[0];
sx q[0];
rz(2.1297852) q[0];
rz(-pi) q[1];
rz(-2.6382006) q[2];
sx q[2];
rz(-2.1648228) q[2];
sx q[2];
rz(0.6436178) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0867566) q[1];
sx q[1];
rz(-2.1385178) q[1];
sx q[1];
rz(2.9556429) q[1];
x q[2];
rz(0.028178111) q[3];
sx q[3];
rz(-2.5163076) q[3];
sx q[3];
rz(0.5042838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9610567) q[2];
sx q[2];
rz(-2.7013216) q[2];
sx q[2];
rz(-0.72009909) q[2];
rz(-2.2911206) q[3];
sx q[3];
rz(-1.974778) q[3];
sx q[3];
rz(-1.7299962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20486031) q[0];
sx q[0];
rz(-2.9731049) q[0];
sx q[0];
rz(2.4334461) q[0];
rz(-1.5180961) q[1];
sx q[1];
rz(-0.93815362) q[1];
sx q[1];
rz(2.785397) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1461648) q[0];
sx q[0];
rz(-1.08279) q[0];
sx q[0];
rz(-1.0900709) q[0];
rz(-1.541186) q[2];
sx q[2];
rz(-0.95311728) q[2];
sx q[2];
rz(-2.7838865) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.5070008) q[1];
sx q[1];
rz(-2.2193877) q[1];
sx q[1];
rz(-1.6084987) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0354) q[3];
sx q[3];
rz(-1.937791) q[3];
sx q[3];
rz(-1.3657686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0492964) q[2];
sx q[2];
rz(-2.3342817) q[2];
sx q[2];
rz(2.2558007) q[2];
rz(-0.93585912) q[3];
sx q[3];
rz(-1.9610619) q[3];
sx q[3];
rz(-2.9203171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7037999) q[0];
sx q[0];
rz(-0.047310345) q[0];
sx q[0];
rz(1.4495151) q[0];
rz(2.119078) q[1];
sx q[1];
rz(-2.6578564) q[1];
sx q[1];
rz(1.6922916) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78010192) q[0];
sx q[0];
rz(-1.5211443) q[0];
sx q[0];
rz(1.9832558) q[0];
x q[1];
rz(-2.7141391) q[2];
sx q[2];
rz(-1.9505672) q[2];
sx q[2];
rz(-2.3878271) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3811581) q[1];
sx q[1];
rz(-1.0236003) q[1];
sx q[1];
rz(0.24914279) q[1];
rz(-1.8362884) q[3];
sx q[3];
rz(-1.42783) q[3];
sx q[3];
rz(0.5393103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8861822) q[2];
sx q[2];
rz(-1.9516727) q[2];
sx q[2];
rz(-0.6748684) q[2];
rz(2.1700962) q[3];
sx q[3];
rz(-2.4392305) q[3];
sx q[3];
rz(1.6500047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3424727) q[0];
sx q[0];
rz(-1.5958888) q[0];
sx q[0];
rz(2.2842443) q[0];
rz(0.2324066) q[1];
sx q[1];
rz(-2.1288165) q[1];
sx q[1];
rz(-1.7850599) q[1];
rz(-0.41396285) q[2];
sx q[2];
rz(-2.064075) q[2];
sx q[2];
rz(1.0309564) q[2];
rz(-1.483199) q[3];
sx q[3];
rz(-0.73053737) q[3];
sx q[3];
rz(-2.4127816) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
