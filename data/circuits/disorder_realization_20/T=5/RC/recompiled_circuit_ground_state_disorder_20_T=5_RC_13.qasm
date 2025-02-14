OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.60431689) q[0];
sx q[0];
rz(3.3831626) q[0];
sx q[0];
rz(9.0917505) q[0];
rz(2.0060519) q[1];
sx q[1];
rz(-0.82692868) q[1];
sx q[1];
rz(-0.64396042) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0892093) q[0];
sx q[0];
rz(-2.174447) q[0];
sx q[0];
rz(-2.8103845) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30006914) q[2];
sx q[2];
rz(-1.4208394) q[2];
sx q[2];
rz(1.7008925) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2702858) q[1];
sx q[1];
rz(-1.8075004) q[1];
sx q[1];
rz(-0.38856296) q[1];
x q[2];
rz(-2.7936739) q[3];
sx q[3];
rz(-1.9609465) q[3];
sx q[3];
rz(-1.8371234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.48646271) q[2];
sx q[2];
rz(-0.4643521) q[2];
sx q[2];
rz(-0.74074024) q[2];
rz(-1.5517392) q[3];
sx q[3];
rz(-2.430075) q[3];
sx q[3];
rz(0.5927425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43316677) q[0];
sx q[0];
rz(-2.0080703) q[0];
sx q[0];
rz(-3.0246227) q[0];
rz(-0.50432694) q[1];
sx q[1];
rz(-1.802899) q[1];
sx q[1];
rz(0.28409827) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.741521) q[0];
sx q[0];
rz(-2.1353545) q[0];
sx q[0];
rz(-0.098640504) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79757046) q[2];
sx q[2];
rz(-1.1402545) q[2];
sx q[2];
rz(2.0883462) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.85728474) q[1];
sx q[1];
rz(-1.4487639) q[1];
sx q[1];
rz(3.1067501) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8900364) q[3];
sx q[3];
rz(-1.3013892) q[3];
sx q[3];
rz(-3.0602853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8154907) q[2];
sx q[2];
rz(-0.88560605) q[2];
sx q[2];
rz(-0.76078129) q[2];
rz(-2.271999) q[3];
sx q[3];
rz(-1.875501) q[3];
sx q[3];
rz(-2.6884955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3608383) q[0];
sx q[0];
rz(-1.1203082) q[0];
sx q[0];
rz(-0.25217062) q[0];
rz(-1.4422656) q[1];
sx q[1];
rz(-0.95528722) q[1];
sx q[1];
rz(2.7853277) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39922455) q[0];
sx q[0];
rz(-0.059905298) q[0];
sx q[0];
rz(1.4673047) q[0];
rz(1.9819471) q[2];
sx q[2];
rz(-2.0817588) q[2];
sx q[2];
rz(-1.9821577) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.60881847) q[1];
sx q[1];
rz(-1.7131117) q[1];
sx q[1];
rz(-2.2102093) q[1];
rz(-0.0019508501) q[3];
sx q[3];
rz(-1.8271128) q[3];
sx q[3];
rz(1.0904097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0226655) q[2];
sx q[2];
rz(-1.7513195) q[2];
sx q[2];
rz(1.6290132) q[2];
rz(0.0096983612) q[3];
sx q[3];
rz(-1.9111948) q[3];
sx q[3];
rz(2.5201918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4951204) q[0];
sx q[0];
rz(-0.54238129) q[0];
sx q[0];
rz(0.25076184) q[0];
rz(1.3465025) q[1];
sx q[1];
rz(-2.612412) q[1];
sx q[1];
rz(-1.4432602) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0899857) q[0];
sx q[0];
rz(-0.7762701) q[0];
sx q[0];
rz(-1.8818186) q[0];
rz(0.94130959) q[2];
sx q[2];
rz(-2.513859) q[2];
sx q[2];
rz(2.8986487) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5574675) q[1];
sx q[1];
rz(-0.87601501) q[1];
sx q[1];
rz(2.4226339) q[1];
rz(-pi) q[2];
rz(2.3425927) q[3];
sx q[3];
rz(-0.36741396) q[3];
sx q[3];
rz(1.1297117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5859588) q[2];
sx q[2];
rz(-1.1052174) q[2];
sx q[2];
rz(-2.8982758) q[2];
rz(0.58468741) q[3];
sx q[3];
rz(-2.5644315) q[3];
sx q[3];
rz(3.1134591) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1481767) q[0];
sx q[0];
rz(-0.36574829) q[0];
sx q[0];
rz(-2.962033) q[0];
rz(1.9163632) q[1];
sx q[1];
rz(-1.7370677) q[1];
sx q[1];
rz(2.6833351) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97044635) q[0];
sx q[0];
rz(-1.2752646) q[0];
sx q[0];
rz(1.6520722) q[0];
rz(-pi) q[1];
x q[1];
rz(1.304428) q[2];
sx q[2];
rz(-0.86037105) q[2];
sx q[2];
rz(1.923234) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51650713) q[1];
sx q[1];
rz(-1.7053776) q[1];
sx q[1];
rz(-0.075304042) q[1];
rz(0.86991258) q[3];
sx q[3];
rz(-1.8509838) q[3];
sx q[3];
rz(2.3126569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3852343) q[2];
sx q[2];
rz(-0.42432722) q[2];
sx q[2];
rz(-2.2198086) q[2];
rz(1.2515757) q[3];
sx q[3];
rz(-1.4082963) q[3];
sx q[3];
rz(3.0799227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.09403041) q[0];
sx q[0];
rz(-0.91218364) q[0];
sx q[0];
rz(2.9929274) q[0];
rz(-2.8246763) q[1];
sx q[1];
rz(-2.6628559) q[1];
sx q[1];
rz(2.5247578) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7148833) q[0];
sx q[0];
rz(-1.6322109) q[0];
sx q[0];
rz(1.4922569) q[0];
x q[1];
rz(2.2451239) q[2];
sx q[2];
rz(-0.54276641) q[2];
sx q[2];
rz(1.3010058) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4839263) q[1];
sx q[1];
rz(-2.1352508) q[1];
sx q[1];
rz(-0.29740833) q[1];
rz(-2.0944164) q[3];
sx q[3];
rz(-1.3304119) q[3];
sx q[3];
rz(2.2464858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2769015) q[2];
sx q[2];
rz(-1.428705) q[2];
sx q[2];
rz(-3.1332664) q[2];
rz(-0.62638038) q[3];
sx q[3];
rz(-0.71599394) q[3];
sx q[3];
rz(0.00055073784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.83547) q[0];
sx q[0];
rz(-1.1518814) q[0];
sx q[0];
rz(2.6595111) q[0];
rz(3.112402) q[1];
sx q[1];
rz(-1.2960641) q[1];
sx q[1];
rz(-0.76404244) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4493443) q[0];
sx q[0];
rz(-1.5034165) q[0];
sx q[0];
rz(1.199556) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.025310658) q[2];
sx q[2];
rz(-0.9111852) q[2];
sx q[2];
rz(0.74756223) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5611539) q[1];
sx q[1];
rz(-0.62404666) q[1];
sx q[1];
rz(1.2777722) q[1];
rz(-pi) q[2];
rz(1.5519616) q[3];
sx q[3];
rz(-0.40098396) q[3];
sx q[3];
rz(2.8249521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.46334106) q[2];
sx q[2];
rz(-1.8222787) q[2];
sx q[2];
rz(2.6574262) q[2];
rz(0.68228996) q[3];
sx q[3];
rz(-0.65410084) q[3];
sx q[3];
rz(3.0068523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084455647) q[0];
sx q[0];
rz(-2.3637922) q[0];
sx q[0];
rz(1.2917668) q[0];
rz(0.46503398) q[1];
sx q[1];
rz(-2.6227622) q[1];
sx q[1];
rz(2.8344287) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7028593) q[0];
sx q[0];
rz(-0.68638681) q[0];
sx q[0];
rz(0.86875654) q[0];
rz(2.1910153) q[2];
sx q[2];
rz(-2.3831316) q[2];
sx q[2];
rz(1.72067) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7502082) q[1];
sx q[1];
rz(-2.5473875) q[1];
sx q[1];
rz(1.2886402) q[1];
rz(-pi) q[2];
rz(-1.5504595) q[3];
sx q[3];
rz(-0.94579783) q[3];
sx q[3];
rz(0.53903264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9610567) q[2];
sx q[2];
rz(-2.7013216) q[2];
sx q[2];
rz(2.4214936) q[2];
rz(-2.2911206) q[3];
sx q[3];
rz(-1.974778) q[3];
sx q[3];
rz(-1.7299962) q[3];
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
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20486031) q[0];
sx q[0];
rz(-2.9731049) q[0];
sx q[0];
rz(2.4334461) q[0];
rz(1.5180961) q[1];
sx q[1];
rz(-2.203439) q[1];
sx q[1];
rz(2.785397) q[1];
rz(-pi/2) q[2];
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
rz(-pi) q[1];
rz(-1.6004066) q[2];
sx q[2];
rz(-2.1884754) q[2];
sx q[2];
rz(-2.7838865) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.5070008) q[1];
sx q[1];
rz(-2.2193877) q[1];
sx q[1];
rz(-1.6084987) q[1];
x q[2];
rz(-0.40608866) q[3];
sx q[3];
rz(-1.1392987) q[3];
sx q[3];
rz(0.0270947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0922962) q[2];
sx q[2];
rz(-0.80731097) q[2];
sx q[2];
rz(-2.2558007) q[2];
rz(-0.93585912) q[3];
sx q[3];
rz(-1.9610619) q[3];
sx q[3];
rz(0.22127557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43779272) q[0];
sx q[0];
rz(-0.047310345) q[0];
sx q[0];
rz(1.4495151) q[0];
rz(-1.0225147) q[1];
sx q[1];
rz(-0.48373628) q[1];
sx q[1];
rz(1.449301) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67774862) q[0];
sx q[0];
rz(-0.41526702) q[0];
sx q[0];
rz(1.4474611) q[0];
x q[1];
rz(0.42745356) q[2];
sx q[2];
rz(-1.9505672) q[2];
sx q[2];
rz(0.7537656) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8359747) q[1];
sx q[1];
rz(-2.545649) q[1];
sx q[1];
rz(1.1862331) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0690156) q[3];
sx q[3];
rz(-0.3007362) q[3];
sx q[3];
rz(1.5141443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.25541043) q[2];
sx q[2];
rz(-1.9516727) q[2];
sx q[2];
rz(2.4667242) q[2];
rz(-2.1700962) q[3];
sx q[3];
rz(-0.70236218) q[3];
sx q[3];
rz(-1.491588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7991199) q[0];
sx q[0];
rz(-1.5958888) q[0];
sx q[0];
rz(2.2842443) q[0];
rz(2.9091861) q[1];
sx q[1];
rz(-1.0127761) q[1];
sx q[1];
rz(1.3565328) q[1];
rz(-0.41396285) q[2];
sx q[2];
rz(-2.064075) q[2];
sx q[2];
rz(1.0309564) q[2];
rz(-1.6583937) q[3];
sx q[3];
rz(-2.4110553) q[3];
sx q[3];
rz(0.72881107) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
