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
rz(-1.1355407) q[1];
sx q[1];
rz(3.9685213) q[1];
sx q[1];
rz(10.068738) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49227958) q[0];
sx q[0];
rz(-0.67848533) q[0];
sx q[0];
rz(1.1301103) q[0];
rz(-pi) q[1];
rz(0.47253387) q[2];
sx q[2];
rz(-2.807155) q[2];
sx q[2];
rz(2.5616733) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.8713069) q[1];
sx q[1];
rz(-1.3340923) q[1];
sx q[1];
rz(2.7530297) q[1];
rz(-pi) q[2];
rz(-1.1584362) q[3];
sx q[3];
rz(-1.8915911) q[3];
sx q[3];
rz(0.12925805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.48646271) q[2];
sx q[2];
rz(-2.6772406) q[2];
sx q[2];
rz(-0.74074024) q[2];
rz(-1.5898534) q[3];
sx q[3];
rz(-0.71151763) q[3];
sx q[3];
rz(0.5927425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
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
rz(-2.8574944) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.741521) q[0];
sx q[0];
rz(-2.1353545) q[0];
sx q[0];
rz(-0.098640504) q[0];
x q[1];
rz(-0.98911907) q[2];
sx q[2];
rz(-0.8627514) q[2];
sx q[2];
rz(3.0281554) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.71775466) q[1];
sx q[1];
rz(-1.6053797) q[1];
sx q[1];
rz(-1.4486905) q[1];
x q[2];
rz(2.8900364) q[3];
sx q[3];
rz(-1.3013892) q[3];
sx q[3];
rz(-0.081307383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8154907) q[2];
sx q[2];
rz(-2.2559866) q[2];
sx q[2];
rz(0.76078129) q[2];
rz(-0.86959362) q[3];
sx q[3];
rz(-1.2660916) q[3];
sx q[3];
rz(-2.6884955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78075439) q[0];
sx q[0];
rz(-2.0212845) q[0];
sx q[0];
rz(0.25217062) q[0];
rz(1.699327) q[1];
sx q[1];
rz(-2.1863054) q[1];
sx q[1];
rz(-2.7853277) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2748791) q[0];
sx q[0];
rz(-1.5646114) q[0];
sx q[0];
rz(-1.6303819) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1596456) q[2];
sx q[2];
rz(-2.0817588) q[2];
sx q[2];
rz(1.9821577) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9911451) q[1];
sx q[1];
rz(-0.65288645) q[1];
sx q[1];
rz(-1.8064524) q[1];
x q[2];
rz(-3.1396418) q[3];
sx q[3];
rz(-1.3144799) q[3];
sx q[3];
rz(1.0904097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0226655) q[2];
sx q[2];
rz(-1.7513195) q[2];
sx q[2];
rz(1.6290132) q[2];
rz(-3.1318943) q[3];
sx q[3];
rz(-1.9111948) q[3];
sx q[3];
rz(2.5201918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4951204) q[0];
sx q[0];
rz(-2.5992114) q[0];
sx q[0];
rz(-2.8908308) q[0];
rz(-1.3465025) q[1];
sx q[1];
rz(-2.612412) q[1];
sx q[1];
rz(1.4432602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051606962) q[0];
sx q[0];
rz(-0.7762701) q[0];
sx q[0];
rz(1.8818186) q[0];
x q[1];
rz(0.40374229) q[2];
sx q[2];
rz(-1.0761257) q[2];
sx q[2];
rz(-2.1658989) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5574675) q[1];
sx q[1];
rz(-0.87601501) q[1];
sx q[1];
rz(-2.4226339) q[1];
rz(-pi) q[2];
rz(0.79899995) q[3];
sx q[3];
rz(-0.36741396) q[3];
sx q[3];
rz(2.011881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5859588) q[2];
sx q[2];
rz(-1.1052174) q[2];
sx q[2];
rz(-2.8982758) q[2];
rz(2.5569052) q[3];
sx q[3];
rz(-2.5644315) q[3];
sx q[3];
rz(-3.1134591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1481767) q[0];
sx q[0];
rz(-2.7758444) q[0];
sx q[0];
rz(2.962033) q[0];
rz(1.2252294) q[1];
sx q[1];
rz(-1.7370677) q[1];
sx q[1];
rz(-2.6833351) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4438547) q[0];
sx q[0];
rz(-0.30618822) q[0];
sx q[0];
rz(2.8809887) q[0];
rz(-pi) q[1];
rz(0.29693895) q[2];
sx q[2];
rz(-0.7504979) q[2];
sx q[2];
rz(2.3194882) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0441664) q[1];
sx q[1];
rz(-1.6454182) q[1];
sx q[1];
rz(-1.4358372) q[1];
rz(0.36009501) q[3];
sx q[3];
rz(-2.2392139) q[3];
sx q[3];
rz(-2.6289712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7563584) q[2];
sx q[2];
rz(-2.7172654) q[2];
sx q[2];
rz(2.2198086) q[2];
rz(1.2515757) q[3];
sx q[3];
rz(-1.7332964) q[3];
sx q[3];
rz(0.061669953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.09403041) q[0];
sx q[0];
rz(-0.91218364) q[0];
sx q[0];
rz(0.14866522) q[0];
rz(2.8246763) q[1];
sx q[1];
rz(-2.6628559) q[1];
sx q[1];
rz(-2.5247578) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5182523) q[0];
sx q[0];
rz(-3.0419311) q[0];
sx q[0];
rz(0.90604337) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2451239) q[2];
sx q[2];
rz(-2.5988262) q[2];
sx q[2];
rz(1.8405869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4839263) q[1];
sx q[1];
rz(-2.1352508) q[1];
sx q[1];
rz(-0.29740833) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27583337) q[3];
sx q[3];
rz(-2.0778928) q[3];
sx q[3];
rz(-0.81229336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2769015) q[2];
sx q[2];
rz(-1.7128877) q[2];
sx q[2];
rz(-3.1332664) q[2];
rz(0.62638038) q[3];
sx q[3];
rz(-0.71599394) q[3];
sx q[3];
rz(-0.00055073784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.83547) q[0];
sx q[0];
rz(-1.9897113) q[0];
sx q[0];
rz(-0.48208153) q[0];
rz(-0.029190633) q[1];
sx q[1];
rz(-1.2960641) q[1];
sx q[1];
rz(2.3775502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0498063) q[0];
sx q[0];
rz(-2.7645664) q[0];
sx q[0];
rz(-1.3868807) q[0];
x q[1];
rz(-3.116282) q[2];
sx q[2];
rz(-0.9111852) q[2];
sx q[2];
rz(-0.74756223) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5611539) q[1];
sx q[1];
rz(-2.517546) q[1];
sx q[1];
rz(1.2777722) q[1];
rz(1.5519616) q[3];
sx q[3];
rz(-0.40098396) q[3];
sx q[3];
rz(2.8249521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6782516) q[2];
sx q[2];
rz(-1.8222787) q[2];
sx q[2];
rz(0.48416644) q[2];
rz(-2.4593027) q[3];
sx q[3];
rz(-2.4874918) q[3];
sx q[3];
rz(-3.0068523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(3.057137) q[0];
sx q[0];
rz(-2.3637922) q[0];
sx q[0];
rz(-1.8498259) q[0];
rz(2.6765587) q[1];
sx q[1];
rz(-2.6227622) q[1];
sx q[1];
rz(0.30716392) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7114176) q[0];
sx q[0];
rz(-1.9924376) q[0];
sx q[0];
rz(1.0118075) q[0];
x q[1];
rz(-0.91395821) q[2];
sx q[2];
rz(-1.1595396) q[2];
sx q[2];
rz(-0.62818254) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0548361) q[1];
sx q[1];
rz(-2.1385178) q[1];
sx q[1];
rz(-0.18594976) q[1];
rz(1.5504595) q[3];
sx q[3];
rz(-0.94579783) q[3];
sx q[3];
rz(-0.53903264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9610567) q[2];
sx q[2];
rz(-0.44027105) q[2];
sx q[2];
rz(0.72009909) q[2];
rz(2.2911206) q[3];
sx q[3];
rz(-1.1668147) q[3];
sx q[3];
rz(1.4115964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9367323) q[0];
sx q[0];
rz(-0.16848773) q[0];
sx q[0];
rz(-0.70814651) q[0];
rz(1.6234966) q[1];
sx q[1];
rz(-2.203439) q[1];
sx q[1];
rz(-2.785397) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66445456) q[0];
sx q[0];
rz(-1.9915446) q[0];
sx q[0];
rz(-0.53945213) q[0];
x q[1];
rz(-1.6004066) q[2];
sx q[2];
rz(-2.1884754) q[2];
sx q[2];
rz(-2.7838865) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5722288) q[1];
sx q[1];
rz(-2.4920643) q[1];
sx q[1];
rz(-3.0919051) q[1];
rz(-pi) q[2];
rz(1.1061927) q[3];
sx q[3];
rz(-1.2038017) q[3];
sx q[3];
rz(-1.7758241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0492964) q[2];
sx q[2];
rz(-0.80731097) q[2];
sx q[2];
rz(0.88579196) q[2];
rz(-2.2057335) q[3];
sx q[3];
rz(-1.9610619) q[3];
sx q[3];
rz(2.9203171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7037999) q[0];
sx q[0];
rz(-0.047310345) q[0];
sx q[0];
rz(-1.6920775) q[0];
rz(-1.0225147) q[1];
sx q[1];
rz(-0.48373628) q[1];
sx q[1];
rz(-1.6922916) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3614907) q[0];
sx q[0];
rz(-1.6204483) q[0];
sx q[0];
rz(-1.1583369) q[0];
rz(0.42745356) q[2];
sx q[2];
rz(-1.1910254) q[2];
sx q[2];
rz(-0.7537656) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.30561799) q[1];
sx q[1];
rz(-0.59594369) q[1];
sx q[1];
rz(-1.1862331) q[1];
rz(-pi) q[2];
rz(-2.0725771) q[3];
sx q[3];
rz(-0.3007362) q[3];
sx q[3];
rz(1.6274483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8861822) q[2];
sx q[2];
rz(-1.9516727) q[2];
sx q[2];
rz(0.6748684) q[2];
rz(-0.97149649) q[3];
sx q[3];
rz(-0.70236218) q[3];
sx q[3];
rz(-1.6500047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3424727) q[0];
sx q[0];
rz(-1.5457038) q[0];
sx q[0];
rz(-0.85734838) q[0];
rz(-0.2324066) q[1];
sx q[1];
rz(-1.0127761) q[1];
sx q[1];
rz(1.3565328) q[1];
rz(-1.0398374) q[2];
sx q[2];
rz(-1.9329484) q[2];
sx q[2];
rz(-0.74495391) q[2];
rz(0.84216778) q[3];
sx q[3];
rz(-1.6292059) q[3];
sx q[3];
rz(-0.90730351) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
