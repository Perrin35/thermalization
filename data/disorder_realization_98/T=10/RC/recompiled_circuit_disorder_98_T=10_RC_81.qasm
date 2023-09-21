OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7712819) q[0];
sx q[0];
rz(-0.9960649) q[0];
sx q[0];
rz(2.2709742) q[0];
rz(-7.3047819) q[1];
sx q[1];
rz(2.8586913) q[1];
sx q[1];
rz(18.999264) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5538841) q[0];
sx q[0];
rz(-2.2061078) q[0];
sx q[0];
rz(0.61650886) q[0];
x q[1];
rz(-2.1120464) q[2];
sx q[2];
rz(-2.0143348) q[2];
sx q[2];
rz(-2.8461547) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7357199) q[1];
sx q[1];
rz(-0.18438965) q[1];
sx q[1];
rz(-1.7806446) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2655067) q[3];
sx q[3];
rz(-2.5698834) q[3];
sx q[3];
rz(-1.3523462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8831138) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(-2.4374938) q[2];
rz(-2.1885833) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927521) q[0];
sx q[0];
rz(-1.623818) q[0];
sx q[0];
rz(-0.63252226) q[0];
rz(2.6951492) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(0.65223637) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3632293) q[0];
sx q[0];
rz(-3.1104381) q[0];
sx q[0];
rz(2.7345783) q[0];
rz(1.863443) q[2];
sx q[2];
rz(-1.7990944) q[2];
sx q[2];
rz(1.6739068) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3995041) q[1];
sx q[1];
rz(-1.3474476) q[1];
sx q[1];
rz(0.29466596) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5370876) q[3];
sx q[3];
rz(-1.0009871) q[3];
sx q[3];
rz(-1.2165716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5941045) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(1.1616421) q[2];
rz(-1.9836327) q[3];
sx q[3];
rz(-2.0722814) q[3];
sx q[3];
rz(1.6903711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0911672) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(2.2913349) q[0];
rz(-2.6440874) q[1];
sx q[1];
rz(-1.9583227) q[1];
sx q[1];
rz(-1.3495548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.772086) q[0];
sx q[0];
rz(-1.8096576) q[0];
sx q[0];
rz(0.68349616) q[0];
rz(-pi) q[1];
rz(-2.0850052) q[2];
sx q[2];
rz(-0.075767013) q[2];
sx q[2];
rz(-2.5839992) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6408491) q[1];
sx q[1];
rz(-1.8060246) q[1];
sx q[1];
rz(0.51847036) q[1];
rz(-1.1447103) q[3];
sx q[3];
rz(-1.5618556) q[3];
sx q[3];
rz(0.63250354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0321908) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(0.90399495) q[2];
rz(0.30113014) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(-1.3999456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6501453) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(1.4105463) q[0];
rz(-2.5097805) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(0.036380336) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7061493) q[0];
sx q[0];
rz(-0.94067803) q[0];
sx q[0];
rz(0.26421996) q[0];
rz(-pi) q[1];
rz(0.47469791) q[2];
sx q[2];
rz(-0.63650741) q[2];
sx q[2];
rz(-1.1324901) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1069146) q[1];
sx q[1];
rz(-0.87768302) q[1];
sx q[1];
rz(0.17130674) q[1];
rz(-pi) q[2];
rz(-2.7500238) q[3];
sx q[3];
rz(-0.25300004) q[3];
sx q[3];
rz(2.5391425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.92695421) q[2];
sx q[2];
rz(-1.4321233) q[2];
sx q[2];
rz(1.1882163) q[2];
rz(-0.67048091) q[3];
sx q[3];
rz(-1.9159578) q[3];
sx q[3];
rz(2.5454583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7017512) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(2.9751076) q[0];
rz(-0.75603756) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(-2.9072445) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4011824) q[0];
sx q[0];
rz(-0.64490841) q[0];
sx q[0];
rz(0.58437225) q[0];
rz(-pi) q[1];
rz(2.064803) q[2];
sx q[2];
rz(-1.6550118) q[2];
sx q[2];
rz(0.46352026) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8824132) q[1];
sx q[1];
rz(-0.83173527) q[1];
sx q[1];
rz(1.3684567) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0626416) q[3];
sx q[3];
rz(-1.8936994) q[3];
sx q[3];
rz(-0.80432804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4328737) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(0.22053545) q[2];
rz(-2.7045414) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(0.76550686) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7261312) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(-0.81714001) q[0];
rz(-0.56610402) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(-1.9979427) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9691539) q[0];
sx q[0];
rz(-1.6108496) q[0];
sx q[0];
rz(-2.7705631) q[0];
rz(2.8191889) q[2];
sx q[2];
rz(-0.69903261) q[2];
sx q[2];
rz(0.93271819) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7856577) q[1];
sx q[1];
rz(-0.4520843) q[1];
sx q[1];
rz(-1.362734) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46623047) q[3];
sx q[3];
rz(-2.5737408) q[3];
sx q[3];
rz(0.0044435244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3879261) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(2.6521519) q[2];
rz(0.22805452) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(-0.42603809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.56931) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(-1.2868767) q[0];
rz(-2.4781748) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(1.2333262) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7147303) q[0];
sx q[0];
rz(-0.38422248) q[0];
sx q[0];
rz(-1.5412488) q[0];
x q[1];
rz(0.8823422) q[2];
sx q[2];
rz(-2.2922278) q[2];
sx q[2];
rz(1.7989858) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.308097) q[1];
sx q[1];
rz(-1.6819681) q[1];
sx q[1];
rz(1.5566467) q[1];
rz(-0.60621467) q[3];
sx q[3];
rz(-2.5698235) q[3];
sx q[3];
rz(1.9619463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.104091) q[2];
sx q[2];
rz(-2.9064894) q[2];
sx q[2];
rz(1.0160149) q[2];
rz(0.070090381) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(-1.0664553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0163517) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(1.6960779) q[0];
rz(-2.9267172) q[1];
sx q[1];
rz(-2.3863249) q[1];
sx q[1];
rz(1.258237) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1930775) q[0];
sx q[0];
rz(-2.6876039) q[0];
sx q[0];
rz(1.2598739) q[0];
rz(-2.2082981) q[2];
sx q[2];
rz(-2.1142695) q[2];
sx q[2];
rz(2.9535434) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.82842365) q[1];
sx q[1];
rz(-1.0272044) q[1];
sx q[1];
rz(-1.6924752) q[1];
x q[2];
rz(2.1636837) q[3];
sx q[3];
rz(-2.3330354) q[3];
sx q[3];
rz(-2.9582994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68226472) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(-0.56274596) q[2];
rz(-0.19966666) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(-1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0257618) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(-0.2510221) q[0];
rz(-2.714278) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(3.1138611) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32556191) q[0];
sx q[0];
rz(-1.0849909) q[0];
sx q[0];
rz(-0.82657878) q[0];
rz(-pi) q[1];
rz(0.8717732) q[2];
sx q[2];
rz(-1.7773526) q[2];
sx q[2];
rz(-0.86824647) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9412649) q[1];
sx q[1];
rz(-2.1695734) q[1];
sx q[1];
rz(1.4501249) q[1];
x q[2];
rz(1.2003044) q[3];
sx q[3];
rz(-3.070188) q[3];
sx q[3];
rz(-2.2927473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.015908265) q[2];
sx q[2];
rz(-0.97970825) q[2];
sx q[2];
rz(-1.1431747) q[2];
rz(-2.9987191) q[3];
sx q[3];
rz(-1.5121127) q[3];
sx q[3];
rz(-0.85723248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3906355) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(2.6877158) q[0];
rz(-0.67165309) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(2.8840816) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47782183) q[0];
sx q[0];
rz(-2.0174694) q[0];
sx q[0];
rz(1.9985984) q[0];
rz(-pi) q[1];
rz(-2.0673413) q[2];
sx q[2];
rz(-1.0814582) q[2];
sx q[2];
rz(-0.41000965) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32138667) q[1];
sx q[1];
rz(-2.4452129) q[1];
sx q[1];
rz(-0.022298261) q[1];
x q[2];
rz(-2.8994843) q[3];
sx q[3];
rz(-2.0383516) q[3];
sx q[3];
rz(0.20413354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3283078) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(-0.60662398) q[2];
rz(0.47484067) q[3];
sx q[3];
rz(-2.1947221) q[3];
sx q[3];
rz(2.301208) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3474779) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(2.9150302) q[1];
sx q[1];
rz(-1.7270052) q[1];
sx q[1];
rz(-2.5546767) q[1];
rz(-0.46802855) q[2];
sx q[2];
rz(-1.9149018) q[2];
sx q[2];
rz(0.55185774) q[2];
rz(-0.26279454) q[3];
sx q[3];
rz(-1.1017208) q[3];
sx q[3];
rz(2.8337939) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];