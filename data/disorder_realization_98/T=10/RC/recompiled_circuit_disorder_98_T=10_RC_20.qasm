OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37031072) q[0];
sx q[0];
rz(4.1376576) q[0];
sx q[0];
rz(7.1538038) q[0];
rz(-1.0215966) q[1];
sx q[1];
rz(-0.28290132) q[1];
sx q[1];
rz(2.9918848) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7604231) q[0];
sx q[0];
rz(-1.086735) q[0];
sx q[0];
rz(2.3056187) q[0];
rz(-0.82586536) q[2];
sx q[2];
rz(-2.4561433) q[2];
sx q[2];
rz(-2.4862188) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.95853165) q[1];
sx q[1];
rz(-1.5325938) q[1];
sx q[1];
rz(1.7512291) q[1];
rz(-pi) q[2];
rz(1.2655067) q[3];
sx q[3];
rz(-0.57170924) q[3];
sx q[3];
rz(-1.3523462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8831138) q[2];
sx q[2];
rz(-3.0443865) q[2];
sx q[2];
rz(2.4374938) q[2];
rz(-2.1885833) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(1.4037508) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927521) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(0.63252226) q[0];
rz(-2.6951492) q[1];
sx q[1];
rz(-1.4182785) q[1];
sx q[1];
rz(-2.4893563) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1855542) q[0];
sx q[0];
rz(-1.5421876) q[0];
sx q[0];
rz(-1.5584598) q[0];
rz(1.863443) q[2];
sx q[2];
rz(-1.3424982) q[2];
sx q[2];
rz(1.4676859) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0374239) q[1];
sx q[1];
rz(-1.8579322) q[1];
sx q[1];
rz(-1.8038521) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2964301) q[3];
sx q[3];
rz(-0.80544986) q[3];
sx q[3];
rz(1.0172539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5941045) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(1.9799505) q[2];
rz(-1.9836327) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(-1.6903711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0911672) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(2.2913349) q[0];
rz(-0.49750528) q[1];
sx q[1];
rz(-1.9583227) q[1];
sx q[1];
rz(-1.7920378) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.772086) q[0];
sx q[0];
rz(-1.8096576) q[0];
sx q[0];
rz(2.4580965) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5047968) q[2];
sx q[2];
rz(-1.5335576) q[2];
sx q[2];
rz(-2.6413692) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68223665) q[1];
sx q[1];
rz(-2.5767234) q[1];
sx q[1];
rz(2.6911246) q[1];
rz(-pi) q[2];
rz(0.0098185929) q[3];
sx q[3];
rz(-1.9968642) q[3];
sx q[3];
rz(2.1992418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0321908) q[2];
sx q[2];
rz(-1.9058062) q[2];
sx q[2];
rz(-2.2375977) q[2];
rz(-2.8404625) q[3];
sx q[3];
rz(-1.3826933) q[3];
sx q[3];
rz(-1.7416471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6501453) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(1.4105463) q[0];
rz(0.63181216) q[1];
sx q[1];
rz(-1.3316863) q[1];
sx q[1];
rz(-0.036380336) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43544337) q[0];
sx q[0];
rz(-0.94067803) q[0];
sx q[0];
rz(2.8773727) q[0];
rz(-2.5600299) q[2];
sx q[2];
rz(-1.8459324) q[2];
sx q[2];
rz(3.0951701) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42602793) q[1];
sx q[1];
rz(-1.702311) q[1];
sx q[1];
rz(2.2711666) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4724457) q[3];
sx q[3];
rz(-1.3373168) q[3];
sx q[3];
rz(1.0055621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.92695421) q[2];
sx q[2];
rz(-1.4321233) q[2];
sx q[2];
rz(-1.1882163) q[2];
rz(-2.4711117) q[3];
sx q[3];
rz(-1.9159578) q[3];
sx q[3];
rz(-2.5454583) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7017512) q[0];
sx q[0];
rz(-1.259946) q[0];
sx q[0];
rz(-2.9751076) q[0];
rz(2.3855551) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(-2.9072445) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4011824) q[0];
sx q[0];
rz(-2.4966842) q[0];
sx q[0];
rz(-2.5572204) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7469823) q[2];
sx q[2];
rz(-0.50054769) q[2];
sx q[2];
rz(-1.2622152) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.25917945) q[1];
sx q[1];
rz(-0.83173527) q[1];
sx q[1];
rz(-1.3684567) q[1];
rz(-pi) q[2];
rz(1.2469532) q[3];
sx q[3];
rz(-1.6456592) q[3];
sx q[3];
rz(2.4002241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70871893) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(0.22053545) q[2];
rz(-0.43705127) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(2.3760858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41546145) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(-2.3244526) q[0];
rz(2.5754886) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(-1.9979427) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1724388) q[0];
sx q[0];
rz(-1.6108496) q[0];
sx q[0];
rz(-0.37102951) q[0];
x q[1];
rz(1.3104865) q[2];
sx q[2];
rz(-2.227265) q[2];
sx q[2];
rz(1.3442163) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.027187849) q[1];
sx q[1];
rz(-1.6611551) q[1];
sx q[1];
rz(-1.1272217) q[1];
rz(-pi) q[2];
rz(-2.6236344) q[3];
sx q[3];
rz(-1.8149788) q[3];
sx q[3];
rz(1.9675919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3879261) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(-0.4894408) q[2];
rz(-0.22805452) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(-2.7155546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.56931) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(-1.8547159) q[0];
rz(2.4781748) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(-1.2333262) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1713283) q[0];
sx q[0];
rz(-1.5818705) q[0];
sx q[0];
rz(1.1867255) q[0];
x q[1];
rz(-0.8823422) q[2];
sx q[2];
rz(-0.8493648) q[2];
sx q[2];
rz(1.7989858) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8804633) q[1];
sx q[1];
rz(-1.5848586) q[1];
sx q[1];
rz(-0.1111828) q[1];
rz(-pi) q[2];
rz(-0.60621467) q[3];
sx q[3];
rz(-2.5698235) q[3];
sx q[3];
rz(-1.1796463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.104091) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(-1.0160149) q[2];
rz(0.070090381) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(-1.0664553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12524097) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(1.6960779) q[0];
rz(2.9267172) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(1.258237) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94851516) q[0];
sx q[0];
rz(-2.6876039) q[0];
sx q[0];
rz(1.2598739) q[0];
x q[1];
rz(-0.6446722) q[2];
sx q[2];
rz(-1.03627) q[2];
sx q[2];
rz(-2.1246186) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5453323) q[1];
sx q[1];
rz(-2.5858871) q[1];
sx q[1];
rz(0.19821367) q[1];
rz(-2.6120841) q[3];
sx q[3];
rz(-2.2141075) q[3];
sx q[3];
rz(2.1852126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4593279) q[2];
sx q[2];
rz(-0.1846281) q[2];
sx q[2];
rz(0.56274596) q[2];
rz(-2.941926) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(-2.0195885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11583081) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(-2.8905706) q[0];
rz(2.714278) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(-3.1138611) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6513034) q[0];
sx q[0];
rz(-0.9285183) q[0];
sx q[0];
rz(2.519033) q[0];
rz(-pi) q[1];
rz(0.26720033) q[2];
sx q[2];
rz(-2.2520817) q[2];
sx q[2];
rz(-2.6097678) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2003277) q[1];
sx q[1];
rz(-0.97201921) q[1];
sx q[1];
rz(1.4501249) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1157007) q[3];
sx q[3];
rz(-1.6373487) q[3];
sx q[3];
rz(-0.47749146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1256844) q[2];
sx q[2];
rz(-0.97970825) q[2];
sx q[2];
rz(1.1431747) q[2];
rz(-0.14287359) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(-0.85723248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7509572) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(-2.6877158) q[0];
rz(-0.67165309) q[1];
sx q[1];
rz(-1.4524873) q[1];
sx q[1];
rz(0.25751105) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2430902) q[0];
sx q[0];
rz(-1.9542964) q[0];
sx q[0];
rz(-0.48454185) q[0];
rz(-pi) q[1];
rz(-2.0673413) q[2];
sx q[2];
rz(-1.0814582) q[2];
sx q[2];
rz(2.731583) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.820206) q[1];
sx q[1];
rz(-0.69637978) q[1];
sx q[1];
rz(3.1192944) q[1];
rz(-pi) q[2];
x q[2];
rz(2.014124) q[3];
sx q[3];
rz(-2.6192198) q[3];
sx q[3];
rz(-0.70538196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8132849) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(2.5349687) q[2];
rz(2.666752) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7941147) q[0];
sx q[0];
rz(-1.62962) q[0];
sx q[0];
rz(-0.99123065) q[0];
rz(-0.22656245) q[1];
sx q[1];
rz(-1.7270052) q[1];
sx q[1];
rz(-2.5546767) q[1];
rz(1.1889585) q[2];
sx q[2];
rz(-1.1321862) q[2];
sx q[2];
rz(1.95375) q[2];
rz(2.0444617) q[3];
sx q[3];
rz(-0.53285014) q[3];
sx q[3];
rz(0.22900029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
