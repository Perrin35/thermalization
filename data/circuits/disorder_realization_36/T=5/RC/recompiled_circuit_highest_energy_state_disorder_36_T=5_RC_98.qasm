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
rz(2.2344196) q[0];
sx q[0];
rz(-0.89651674) q[0];
sx q[0];
rz(-0.040891115) q[0];
rz(-4.4019051) q[1];
sx q[1];
rz(2.6949096) q[1];
sx q[1];
rz(6.3781368) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0329566) q[0];
sx q[0];
rz(-2.6755736) q[0];
sx q[0];
rz(0.52435438) q[0];
rz(-2.7822438) q[2];
sx q[2];
rz(-2.5541325) q[2];
sx q[2];
rz(2.0094144) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.59254348) q[1];
sx q[1];
rz(-2.2993339) q[1];
sx q[1];
rz(1.0647573) q[1];
x q[2];
rz(2.0435752) q[3];
sx q[3];
rz(-1.9736787) q[3];
sx q[3];
rz(1.1927746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.2824668) q[2];
sx q[2];
rz(-2.5390415) q[2];
sx q[2];
rz(-1.1653384) q[2];
rz(0.91853842) q[3];
sx q[3];
rz(-0.10401741) q[3];
sx q[3];
rz(1.2134086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.7391881) q[0];
sx q[0];
rz(-0.63888752) q[0];
sx q[0];
rz(-0.30493394) q[0];
rz(1.0324427) q[1];
sx q[1];
rz(-2.86125) q[1];
sx q[1];
rz(0.84411821) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88972487) q[0];
sx q[0];
rz(-1.8093523) q[0];
sx q[0];
rz(0.15695928) q[0];
x q[1];
rz(1.779489) q[2];
sx q[2];
rz(-1.450453) q[2];
sx q[2];
rz(-2.142765) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5836778) q[1];
sx q[1];
rz(-0.47835438) q[1];
sx q[1];
rz(-0.5443404) q[1];
rz(2.5145766) q[3];
sx q[3];
rz(-1.5582784) q[3];
sx q[3];
rz(2.1345958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.9742763) q[2];
sx q[2];
rz(-2.4582489) q[2];
sx q[2];
rz(-2.5565476) q[2];
rz(-0.85917464) q[3];
sx q[3];
rz(-1.1480568) q[3];
sx q[3];
rz(1.1332716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1187196) q[0];
sx q[0];
rz(-0.82094231) q[0];
sx q[0];
rz(3.0430479) q[0];
rz(0.56602829) q[1];
sx q[1];
rz(-2.2955344) q[1];
sx q[1];
rz(0.50618323) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7166185) q[0];
sx q[0];
rz(-2.0831265) q[0];
sx q[0];
rz(1.1879735) q[0];
rz(-2.7012791) q[2];
sx q[2];
rz(-0.8699421) q[2];
sx q[2];
rz(0.93849692) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.652388) q[1];
sx q[1];
rz(-1.4474157) q[1];
sx q[1];
rz(1.3493933) q[1];
x q[2];
rz(-0.39301707) q[3];
sx q[3];
rz(-2.1073494) q[3];
sx q[3];
rz(2.4483333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11518662) q[2];
sx q[2];
rz(-1.2135442) q[2];
sx q[2];
rz(3.0565267) q[2];
rz(-2.3885942) q[3];
sx q[3];
rz(-2.4716061) q[3];
sx q[3];
rz(1.6116713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6104777) q[0];
sx q[0];
rz(-1.5014638) q[0];
sx q[0];
rz(-2.6016972) q[0];
rz(-3.1119697) q[1];
sx q[1];
rz(-2.3952775) q[1];
sx q[1];
rz(1.354904) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52068096) q[0];
sx q[0];
rz(-1.3518466) q[0];
sx q[0];
rz(-0.82100533) q[0];
rz(-pi) q[1];
rz(-0.37432713) q[2];
sx q[2];
rz(-0.74639635) q[2];
sx q[2];
rz(-2.057892) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3554509) q[1];
sx q[1];
rz(-2.1340573) q[1];
sx q[1];
rz(1.2652629) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79014312) q[3];
sx q[3];
rz(-1.6340268) q[3];
sx q[3];
rz(0.15925285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4386091) q[2];
sx q[2];
rz(-2.1712124) q[2];
sx q[2];
rz(-0.89540974) q[2];
rz(1.4488975) q[3];
sx q[3];
rz(-1.9796895) q[3];
sx q[3];
rz(-2.7321775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2153636) q[0];
sx q[0];
rz(-2.138593) q[0];
sx q[0];
rz(-0.0005501752) q[0];
rz(0.5254566) q[1];
sx q[1];
rz(-2.2976687) q[1];
sx q[1];
rz(2.1702683) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4068953) q[0];
sx q[0];
rz(-2.3427123) q[0];
sx q[0];
rz(-2.6968234) q[0];
rz(-pi) q[1];
rz(-1.9840365) q[2];
sx q[2];
rz(-1.8964502) q[2];
sx q[2];
rz(2.535274) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3974053) q[1];
sx q[1];
rz(-1.0013072) q[1];
sx q[1];
rz(-0.90267148) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9266403) q[3];
sx q[3];
rz(-1.264241) q[3];
sx q[3];
rz(-2.7988899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1725258) q[2];
sx q[2];
rz(-2.064164) q[2];
sx q[2];
rz(-1.7871008) q[2];
rz(-0.6012249) q[3];
sx q[3];
rz(-2.0982592) q[3];
sx q[3];
rz(-1.575298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4300267) q[0];
sx q[0];
rz(-0.30421782) q[0];
sx q[0];
rz(0.29997224) q[0];
rz(-0.9777588) q[1];
sx q[1];
rz(-2.561196) q[1];
sx q[1];
rz(-2.207644) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47563206) q[0];
sx q[0];
rz(-0.36132672) q[0];
sx q[0];
rz(2.7368109) q[0];
rz(-2.0522473) q[2];
sx q[2];
rz(-1.0100968) q[2];
sx q[2];
rz(-2.581832) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0505414) q[1];
sx q[1];
rz(-0.81650298) q[1];
sx q[1];
rz(1.8920808) q[1];
x q[2];
rz(1.2190231) q[3];
sx q[3];
rz(-1.9493503) q[3];
sx q[3];
rz(-1.5768676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0705491) q[2];
sx q[2];
rz(-2.1328378) q[2];
sx q[2];
rz(-1.1318995) q[2];
rz(-1.2567358) q[3];
sx q[3];
rz(-2.3102424) q[3];
sx q[3];
rz(-1.7251714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32135949) q[0];
sx q[0];
rz(-1.2437404) q[0];
sx q[0];
rz(2.5309122) q[0];
rz(0.75675476) q[1];
sx q[1];
rz(-0.38833955) q[1];
sx q[1];
rz(-1.4974219) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0038549) q[0];
sx q[0];
rz(-0.12659368) q[0];
sx q[0];
rz(2.9403482) q[0];
rz(-0.87074222) q[2];
sx q[2];
rz(-2.6079834) q[2];
sx q[2];
rz(-2.6117532) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.25631902) q[1];
sx q[1];
rz(-2.4838964) q[1];
sx q[1];
rz(-2.9672876) q[1];
rz(-1.6711216) q[3];
sx q[3];
rz(-0.85160321) q[3];
sx q[3];
rz(2.6503785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18621592) q[2];
sx q[2];
rz(-0.66408855) q[2];
sx q[2];
rz(0.60626283) q[2];
rz(0.22054211) q[3];
sx q[3];
rz(-1.8044148) q[3];
sx q[3];
rz(0.36673275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99899387) q[0];
sx q[0];
rz(-1.4936916) q[0];
sx q[0];
rz(-0.9077453) q[0];
rz(1.6084464) q[1];
sx q[1];
rz(-2.1870859) q[1];
sx q[1];
rz(3.122701) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.205446) q[0];
sx q[0];
rz(-0.48580446) q[0];
sx q[0];
rz(-2.8914408) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18757815) q[2];
sx q[2];
rz(-2.8370428) q[2];
sx q[2];
rz(1.2414602) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0844473) q[1];
sx q[1];
rz(-1.6020163) q[1];
sx q[1];
rz(-2.679326) q[1];
rz(-pi) q[2];
rz(2.5653061) q[3];
sx q[3];
rz(-2.2243847) q[3];
sx q[3];
rz(1.5022851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3733526) q[2];
sx q[2];
rz(-1.5529996) q[2];
sx q[2];
rz(0.24660435) q[2];
rz(1.2982347) q[3];
sx q[3];
rz(-1.6876829) q[3];
sx q[3];
rz(-1.8728144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3835555) q[0];
sx q[0];
rz(-2.0691431) q[0];
sx q[0];
rz(-2.2776336) q[0];
rz(1.2147238) q[1];
sx q[1];
rz(-1.1178958) q[1];
sx q[1];
rz(-0.58089677) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6462517) q[0];
sx q[0];
rz(-0.87336191) q[0];
sx q[0];
rz(2.4861885) q[0];
rz(-pi) q[1];
rz(1.5283818) q[2];
sx q[2];
rz(-0.73684947) q[2];
sx q[2];
rz(0.15351099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.249596) q[1];
sx q[1];
rz(-0.70887762) q[1];
sx q[1];
rz(-2.6061771) q[1];
rz(0.77928752) q[3];
sx q[3];
rz(-1.3715991) q[3];
sx q[3];
rz(-1.3785386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8860127) q[2];
sx q[2];
rz(-0.46147999) q[2];
sx q[2];
rz(0.18730051) q[2];
rz(-2.1646132) q[3];
sx q[3];
rz(-1.5005485) q[3];
sx q[3];
rz(1.8627082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63000694) q[0];
sx q[0];
rz(-0.66660175) q[0];
sx q[0];
rz(1.8774207) q[0];
rz(-2.1695747) q[1];
sx q[1];
rz(-0.7518026) q[1];
sx q[1];
rz(1.972563) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3384333) q[0];
sx q[0];
rz(-1.9040603) q[0];
sx q[0];
rz(-0.14305556) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3832267) q[2];
sx q[2];
rz(-2.6347199) q[2];
sx q[2];
rz(-0.015794347) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0824521) q[1];
sx q[1];
rz(-2.2603717) q[1];
sx q[1];
rz(-0.13549681) q[1];
rz(2.8037352) q[3];
sx q[3];
rz(-0.33379972) q[3];
sx q[3];
rz(2.5797082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6675889) q[2];
sx q[2];
rz(-1.7550125) q[2];
sx q[2];
rz(-0.37707314) q[2];
rz(0.22370473) q[3];
sx q[3];
rz(-2.5535899) q[3];
sx q[3];
rz(-1.912775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8857464) q[0];
sx q[0];
rz(-1.3930014) q[0];
sx q[0];
rz(-0.2572671) q[0];
rz(-2.5480351) q[1];
sx q[1];
rz(-2.3496353) q[1];
sx q[1];
rz(-1.4812462) q[1];
rz(-0.72102265) q[2];
sx q[2];
rz(-0.99698721) q[2];
sx q[2];
rz(0.1884603) q[2];
rz(2.138924) q[3];
sx q[3];
rz(-2.6780861) q[3];
sx q[3];
rz(-1.1939315) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
