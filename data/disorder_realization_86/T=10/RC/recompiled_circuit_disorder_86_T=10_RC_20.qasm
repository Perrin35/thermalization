OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7744301) q[0];
sx q[0];
rz(-0.91355938) q[0];
sx q[0];
rz(1.4120742) q[0];
rz(-2.9867759) q[1];
sx q[1];
rz(-0.59564367) q[1];
sx q[1];
rz(-1.6593978) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55573758) q[0];
sx q[0];
rz(-1.1475539) q[0];
sx q[0];
rz(1.7413571) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6719195) q[2];
sx q[2];
rz(-0.28684068) q[2];
sx q[2];
rz(-1.4925721) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2130148) q[1];
sx q[1];
rz(-1.9933812) q[1];
sx q[1];
rz(-1.0282474) q[1];
rz(0.10098884) q[3];
sx q[3];
rz(-1.0102934) q[3];
sx q[3];
rz(1.3274173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.98510629) q[2];
sx q[2];
rz(-2.6323695) q[2];
sx q[2];
rz(0.86581725) q[2];
rz(2.1872897) q[3];
sx q[3];
rz(-1.538397) q[3];
sx q[3];
rz(-1.8538063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99825478) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(0.026219333) q[0];
rz(1.5401309) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(-2.1781133) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022097691) q[0];
sx q[0];
rz(-0.95675981) q[0];
sx q[0];
rz(3.1387781) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8112667) q[2];
sx q[2];
rz(-2.1147554) q[2];
sx q[2];
rz(-2.3842173) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3669489) q[1];
sx q[1];
rz(-1.0732713) q[1];
sx q[1];
rz(0.36555396) q[1];
x q[2];
rz(-0.24641896) q[3];
sx q[3];
rz(-1.911947) q[3];
sx q[3];
rz(-2.6951172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6271237) q[2];
sx q[2];
rz(-1.1273948) q[2];
sx q[2];
rz(-0.13452402) q[2];
rz(0.7450122) q[3];
sx q[3];
rz(-0.22694215) q[3];
sx q[3];
rz(-2.1988595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9298252) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(2.3441558) q[0];
rz(-1.047661) q[1];
sx q[1];
rz(-0.14973775) q[1];
sx q[1];
rz(2.581596) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6277498) q[0];
sx q[0];
rz(-1.3686413) q[0];
sx q[0];
rz(-1.1802799) q[0];
rz(-0.30324869) q[2];
sx q[2];
rz(-1.5911284) q[2];
sx q[2];
rz(0.081239935) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8078976) q[1];
sx q[1];
rz(-1.9296608) q[1];
sx q[1];
rz(-1.3586504) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0873763) q[3];
sx q[3];
rz(-1.9994352) q[3];
sx q[3];
rz(-0.40070686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3893163) q[2];
sx q[2];
rz(-1.1976778) q[2];
sx q[2];
rz(0.17253549) q[2];
rz(-0.98207384) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(-1.0579695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(2.003222) q[0];
sx q[0];
rz(-2.0439742) q[0];
sx q[0];
rz(2.8570783) q[0];
rz(0.31670397) q[1];
sx q[1];
rz(-2.7088294) q[1];
sx q[1];
rz(-1.2987312) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5772229) q[0];
sx q[0];
rz(-1.3457314) q[0];
sx q[0];
rz(-2.0805012) q[0];
rz(-pi) q[1];
rz(-0.80438517) q[2];
sx q[2];
rz(-1.569869) q[2];
sx q[2];
rz(1.5915807) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.87151566) q[1];
sx q[1];
rz(-1.75711) q[1];
sx q[1];
rz(-0.99888505) q[1];
x q[2];
rz(3.0292547) q[3];
sx q[3];
rz(-1.2445407) q[3];
sx q[3];
rz(2.5374967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6056885) q[2];
sx q[2];
rz(-0.19583344) q[2];
sx q[2];
rz(-2.7569125) q[2];
rz(-0.7540594) q[3];
sx q[3];
rz(-2.0850756) q[3];
sx q[3];
rz(-1.4543021) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5383179) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(1.3866562) q[0];
rz(0.23100135) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(2.8447661) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41243991) q[0];
sx q[0];
rz(-1.531633) q[0];
sx q[0];
rz(-0.57106437) q[0];
x q[1];
rz(-0.16634059) q[2];
sx q[2];
rz(-0.74976774) q[2];
sx q[2];
rz(-1.2321842) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0592812) q[1];
sx q[1];
rz(-2.0356405) q[1];
sx q[1];
rz(2.6962198) q[1];
x q[2];
rz(0.71367587) q[3];
sx q[3];
rz(-1.5034961) q[3];
sx q[3];
rz(-2.0590559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7632873) q[2];
sx q[2];
rz(-1.3102691) q[2];
sx q[2];
rz(-0.39247593) q[2];
rz(1.9893507) q[3];
sx q[3];
rz(-2.4270054) q[3];
sx q[3];
rz(0.31744441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-1.0157938) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(0.75138599) q[0];
rz(1.8136576) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(2.5352535) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4088926) q[0];
sx q[0];
rz(-0.048763976) q[0];
sx q[0];
rz(-0.34838895) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7251284) q[2];
sx q[2];
rz(-2.205924) q[2];
sx q[2];
rz(-3.0481899) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.50748435) q[1];
sx q[1];
rz(-2.1180696) q[1];
sx q[1];
rz(-2.9299111) q[1];
rz(-pi) q[2];
rz(1.8984406) q[3];
sx q[3];
rz(-0.22938211) q[3];
sx q[3];
rz(2.4250507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.63885826) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(1.139337) q[2];
rz(1.4849439) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(-3.0373354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5320324) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(0.50810057) q[0];
rz(1.5787026) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(-0.79024822) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7471874) q[0];
sx q[0];
rz(-0.98441511) q[0];
sx q[0];
rz(-1.9275083) q[0];
x q[1];
rz(1.5590645) q[2];
sx q[2];
rz(-1.3571697) q[2];
sx q[2];
rz(0.27660433) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3527457) q[1];
sx q[1];
rz(-2.6585796) q[1];
sx q[1];
rz(2.5301945) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4942057) q[3];
sx q[3];
rz(-1.8772519) q[3];
sx q[3];
rz(-2.2341773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.053085176) q[2];
sx q[2];
rz(-2.6997456) q[2];
sx q[2];
rz(1.4132168) q[2];
rz(-1.6648071) q[3];
sx q[3];
rz(-2.1063185) q[3];
sx q[3];
rz(-0.061554519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7168032) q[0];
sx q[0];
rz(-3.1112818) q[0];
sx q[0];
rz(1.0472263) q[0];
rz(0.60910243) q[1];
sx q[1];
rz(-1.4139688) q[1];
sx q[1];
rz(1.75288) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7097276) q[0];
sx q[0];
rz(-1.0114397) q[0];
sx q[0];
rz(0.035168408) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65865626) q[2];
sx q[2];
rz(-2.5052862) q[2];
sx q[2];
rz(-0.051740019) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3125004) q[1];
sx q[1];
rz(-1.478985) q[1];
sx q[1];
rz(-1.5593668) q[1];
rz(-pi) q[2];
rz(0.29626131) q[3];
sx q[3];
rz(-2.1564266) q[3];
sx q[3];
rz(-1.1405917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9528815) q[2];
sx q[2];
rz(-2.7313576) q[2];
sx q[2];
rz(2.2593373) q[2];
rz(-1.4011718) q[3];
sx q[3];
rz(-1.975235) q[3];
sx q[3];
rz(1.9410979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-2.8700478) q[0];
sx q[0];
rz(-2.7247868) q[0];
sx q[0];
rz(-1.4260938) q[0];
rz(-3.0601314) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(0.55823278) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4687913) q[0];
sx q[0];
rz(-1.9580012) q[0];
sx q[0];
rz(2.0331435) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4633281) q[2];
sx q[2];
rz(-1.6245914) q[2];
sx q[2];
rz(-1.4449643) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.25591125) q[1];
sx q[1];
rz(-0.7464039) q[1];
sx q[1];
rz(-1.6649151) q[1];
rz(-pi) q[2];
rz(-0.48645143) q[3];
sx q[3];
rz(-1.4064186) q[3];
sx q[3];
rz(-2.1742976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2925064) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(-0.212184) q[2];
rz(2.9296181) q[3];
sx q[3];
rz(-0.68325716) q[3];
sx q[3];
rz(1.9395444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6367209) q[0];
sx q[0];
rz(-2.3175406) q[0];
sx q[0];
rz(1.5378392) q[0];
rz(2.3161855) q[1];
sx q[1];
rz(-2.4688265) q[1];
sx q[1];
rz(2.6182981) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16738811) q[0];
sx q[0];
rz(-0.61199576) q[0];
sx q[0];
rz(-0.1083072) q[0];
x q[1];
rz(-0.043838219) q[2];
sx q[2];
rz(-2.1423116) q[2];
sx q[2];
rz(0.29585719) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1112422) q[1];
sx q[1];
rz(-2.6551464) q[1];
sx q[1];
rz(-2.4187947) q[1];
rz(2.8966122) q[3];
sx q[3];
rz(-1.3462726) q[3];
sx q[3];
rz(-2.6905439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3045197) q[2];
sx q[2];
rz(-2.4261116) q[2];
sx q[2];
rz(-0.26930299) q[2];
rz(0.4942016) q[3];
sx q[3];
rz(-0.84635693) q[3];
sx q[3];
rz(2.0555029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3257278) q[0];
sx q[0];
rz(-1.5300735) q[0];
sx q[0];
rz(-1.6515401) q[0];
rz(-1.4670463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(-1.6124484) q[2];
sx q[2];
rz(-2.8786082) q[2];
sx q[2];
rz(-1.0515121) q[2];
rz(0.23444093) q[3];
sx q[3];
rz(-0.81080484) q[3];
sx q[3];
rz(-0.57639359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
