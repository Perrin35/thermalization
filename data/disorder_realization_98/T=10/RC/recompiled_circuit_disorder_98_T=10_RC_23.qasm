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
rz(-0.87061849) q[0];
rz(-7.3047819) q[1];
sx q[1];
rz(2.8586913) q[1];
sx q[1];
rz(18.999264) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5877085) q[0];
sx q[0];
rz(-2.2061078) q[0];
sx q[0];
rz(0.61650886) q[0];
rz(-2.6354191) q[2];
sx q[2];
rz(-2.0548327) q[2];
sx q[2];
rz(1.6137705) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40587273) q[1];
sx q[1];
rz(-2.957203) q[1];
sx q[1];
rz(1.3609481) q[1];
rz(-pi) q[2];
rz(-1.0204131) q[3];
sx q[3];
rz(-1.4074416) q[3];
sx q[3];
rz(2.6640716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.25847882) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(0.70409888) q[2];
rz(-0.95300931) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(-1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927521) q[0];
sx q[0];
rz(-1.623818) q[0];
sx q[0];
rz(-0.63252226) q[0];
rz(0.44644341) q[1];
sx q[1];
rz(-1.4182785) q[1];
sx q[1];
rz(-2.4893563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5271877) q[0];
sx q[0];
rz(-1.5584649) q[0];
sx q[0];
rz(-3.1129818) q[0];
x q[1];
rz(-2.2488238) q[2];
sx q[2];
rz(-2.7724578) q[2];
sx q[2];
rz(0.54112753) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7420885) q[1];
sx q[1];
rz(-1.794145) q[1];
sx q[1];
rz(-0.29466596) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84516256) q[3];
sx q[3];
rz(-2.3361428) q[3];
sx q[3];
rz(1.0172539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5941045) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(-1.1616421) q[2];
rz(1.15796) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(1.4512216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0911672) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(-0.85025775) q[0];
rz(-0.49750528) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(1.7920378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.772086) q[0];
sx q[0];
rz(-1.8096576) q[0];
sx q[0];
rz(-2.4580965) q[0];
rz(-pi) q[1];
x q[1];
rz(0.037319855) q[2];
sx q[2];
rz(-1.63675) q[2];
sx q[2];
rz(-1.0730336) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.68223665) q[1];
sx q[1];
rz(-2.5767234) q[1];
sx q[1];
rz(0.45046803) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9968824) q[3];
sx q[3];
rz(-1.5797371) q[3];
sx q[3];
rz(-0.63250354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1094018) q[2];
sx q[2];
rz(-1.9058062) q[2];
sx q[2];
rz(0.90399495) q[2];
rz(-0.30113014) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(1.3999456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6501453) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(-1.4105463) q[0];
rz(-2.5097805) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(-3.1052123) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43544337) q[0];
sx q[0];
rz(-2.2009146) q[0];
sx q[0];
rz(0.26421996) q[0];
rz(-2.6668947) q[2];
sx q[2];
rz(-2.5050852) q[2];
sx q[2];
rz(-2.0091025) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.42602793) q[1];
sx q[1];
rz(-1.4392816) q[1];
sx q[1];
rz(2.2711666) q[1];
rz(-pi) q[2];
rz(1.4724457) q[3];
sx q[3];
rz(-1.3373168) q[3];
sx q[3];
rz(2.1360306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2146384) q[2];
sx q[2];
rz(-1.4321233) q[2];
sx q[2];
rz(1.9533763) q[2];
rz(-0.67048091) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(0.59613434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7017512) q[0];
sx q[0];
rz(-1.259946) q[0];
sx q[0];
rz(-0.16648509) q[0];
rz(0.75603756) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(-0.23434815) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4011824) q[0];
sx q[0];
rz(-2.4966842) q[0];
sx q[0];
rz(-2.5572204) q[0];
rz(-0.095586153) q[2];
sx q[2];
rz(-1.0786973) q[2];
sx q[2];
rz(-2.0795859) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8824132) q[1];
sx q[1];
rz(-2.3098574) q[1];
sx q[1];
rz(1.773136) q[1];
rz(-1.8946394) q[3];
sx q[3];
rz(-1.6456592) q[3];
sx q[3];
rz(2.4002241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.70871893) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(2.9210572) q[2];
rz(2.7045414) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(2.3760858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41546145) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(2.3244526) q[0];
rz(2.5754886) q[1];
sx q[1];
rz(-1.7931466) q[1];
sx q[1];
rz(-1.1436499) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2957942) q[0];
sx q[0];
rz(-2.7685071) q[0];
sx q[0];
rz(-0.11008115) q[0];
x q[1];
rz(2.8191889) q[2];
sx q[2];
rz(-2.44256) q[2];
sx q[2];
rz(-0.93271819) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3559349) q[1];
sx q[1];
rz(-0.4520843) q[1];
sx q[1];
rz(-1.7788586) q[1];
rz(1.2915217) q[3];
sx q[3];
rz(-1.0696628) q[3];
sx q[3];
rz(2.6078893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3879261) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(-2.6521519) q[2];
rz(0.22805452) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(-0.42603809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.56931) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(1.8547159) q[0];
rz(-2.4781748) q[1];
sx q[1];
rz(-1.5692915) q[1];
sx q[1];
rz(1.9082665) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7147303) q[0];
sx q[0];
rz(-0.38422248) q[0];
sx q[0];
rz(1.5412488) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2592505) q[2];
sx q[2];
rz(-0.8493648) q[2];
sx q[2];
rz(-1.7989858) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8334956) q[1];
sx q[1];
rz(-1.6819681) q[1];
sx q[1];
rz(1.5566467) q[1];
rz(-pi) q[2];
rz(0.60621467) q[3];
sx q[3];
rz(-2.5698235) q[3];
sx q[3];
rz(1.1796463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.104091) q[2];
sx q[2];
rz(-2.9064894) q[2];
sx q[2];
rz(-2.1255778) q[2];
rz(3.0715023) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(-2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12524097) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(-1.6960779) q[0];
rz(-0.21487543) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(1.258237) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94851516) q[0];
sx q[0];
rz(-0.45398871) q[0];
sx q[0];
rz(1.2598739) q[0];
rz(-0.6446722) q[2];
sx q[2];
rz(-2.1053227) q[2];
sx q[2];
rz(-1.0169741) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5453323) q[1];
sx q[1];
rz(-2.5858871) q[1];
sx q[1];
rz(0.19821367) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6120841) q[3];
sx q[3];
rz(-2.2141075) q[3];
sx q[3];
rz(-2.1852126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4593279) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(0.56274596) q[2];
rz(2.941926) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(2.0195885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11583081) q[0];
sx q[0];
rz(-2.0985726) q[0];
sx q[0];
rz(2.8905706) q[0];
rz(-0.42731467) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(0.02773157) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32556191) q[0];
sx q[0];
rz(-2.0566018) q[0];
sx q[0];
rz(2.3150139) q[0];
rz(-pi) q[1];
rz(-1.8856144) q[2];
sx q[2];
rz(-2.417649) q[2];
sx q[2];
rz(0.9418504) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30222826) q[1];
sx q[1];
rz(-1.4711958) q[1];
sx q[1];
rz(-2.539413) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9412882) q[3];
sx q[3];
rz(-0.071404608) q[3];
sx q[3];
rz(-2.2927473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.015908265) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(1.1431747) q[2];
rz(0.14287359) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(-2.2843602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3906355) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(-0.45387682) q[0];
rz(0.67165309) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(0.25751105) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6637708) q[0];
sx q[0];
rz(-1.1241233) q[0];
sx q[0];
rz(1.9985984) q[0];
x q[1];
rz(-2.0673413) q[2];
sx q[2];
rz(-1.0814582) q[2];
sx q[2];
rz(-0.41000965) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35044893) q[1];
sx q[1];
rz(-0.87462438) q[1];
sx q[1];
rz(1.552156) q[1];
rz(-pi) q[2];
rz(1.1274687) q[3];
sx q[3];
rz(-2.6192198) q[3];
sx q[3];
rz(-2.4362107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3283078) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(0.60662398) q[2];
rz(-2.666752) q[3];
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
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7941147) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(0.22656245) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(-0.46802855) q[2];
sx q[2];
rz(-1.9149018) q[2];
sx q[2];
rz(0.55185774) q[2];
rz(0.26279454) q[3];
sx q[3];
rz(-2.0398718) q[3];
sx q[3];
rz(-0.3077988) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
