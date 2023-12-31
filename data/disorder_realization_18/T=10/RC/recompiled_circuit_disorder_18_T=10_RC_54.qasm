OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3192531) q[0];
sx q[0];
rz(-2.1847794) q[0];
sx q[0];
rz(1.7083038) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(-1.911093) q[1];
sx q[1];
rz(1.9967611) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24554907) q[0];
sx q[0];
rz(-2.890812) q[0];
sx q[0];
rz(-1.8401237) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8561864) q[2];
sx q[2];
rz(-1.1661582) q[2];
sx q[2];
rz(-2.0602496) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.45935985) q[1];
sx q[1];
rz(-1.1182922) q[1];
sx q[1];
rz(-0.020999055) q[1];
rz(-0.47031109) q[3];
sx q[3];
rz(-0.59578005) q[3];
sx q[3];
rz(2.6252928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.11674374) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(-1.2343181) q[2];
rz(-1.0553137) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9556483) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(0.82988513) q[0];
rz(-2.8886967) q[1];
sx q[1];
rz(-2.3650832) q[1];
sx q[1];
rz(-2.2944962) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9053099) q[0];
sx q[0];
rz(-2.5295527) q[0];
sx q[0];
rz(0.73487868) q[0];
rz(-pi) q[1];
rz(-0.29909889) q[2];
sx q[2];
rz(-2.3397589) q[2];
sx q[2];
rz(-0.15660827) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3165247) q[1];
sx q[1];
rz(-1.0267338) q[1];
sx q[1];
rz(-1.6506509) q[1];
rz(-pi) q[2];
rz(2.0081677) q[3];
sx q[3];
rz(-1.2542033) q[3];
sx q[3];
rz(-2.7377627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.117924) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(2.4222597) q[2];
rz(-1.8524648) q[3];
sx q[3];
rz(-2.9057861) q[3];
sx q[3];
rz(1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8711202) q[0];
sx q[0];
rz(-1.8914762) q[0];
sx q[0];
rz(2.5701994) q[0];
rz(2.1748623) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(-2.6142696) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8427211) q[0];
sx q[0];
rz(-1.5639595) q[0];
sx q[0];
rz(-0.30928916) q[0];
rz(-pi) q[1];
rz(0.40076077) q[2];
sx q[2];
rz(-0.71360613) q[2];
sx q[2];
rz(2.1378627) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1151162) q[1];
sx q[1];
rz(-2.3079702) q[1];
sx q[1];
rz(-1.2181746) q[1];
rz(2.7819488) q[3];
sx q[3];
rz(-0.76562866) q[3];
sx q[3];
rz(-1.5776724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8335235) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(2.4148338) q[2];
rz(-2.1440992) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.5184901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961287) q[0];
sx q[0];
rz(-1.6587057) q[0];
sx q[0];
rz(-1.0035275) q[0];
rz(3.1009122) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(2.3153268) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4429312) q[0];
sx q[0];
rz(-2.1974265) q[0];
sx q[0];
rz(2.4311964) q[0];
rz(-pi) q[1];
rz(0.063886558) q[2];
sx q[2];
rz(-0.70110828) q[2];
sx q[2];
rz(1.0039767) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1288209) q[1];
sx q[1];
rz(-1.2342617) q[1];
sx q[1];
rz(-3.0613042) q[1];
rz(-1.9598947) q[3];
sx q[3];
rz(-2.1660921) q[3];
sx q[3];
rz(2.1967595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5740009) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(2.2272002) q[2];
rz(3.0363723) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(0.58925327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.88702622) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(-2.0078833) q[0];
rz(-0.0056313593) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(-0.61757913) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78695801) q[0];
sx q[0];
rz(-1.5411396) q[0];
sx q[0];
rz(-1.6164854) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5768443) q[2];
sx q[2];
rz(-0.72482938) q[2];
sx q[2];
rz(-0.85541475) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6981814) q[1];
sx q[1];
rz(-0.65084208) q[1];
sx q[1];
rz(-1.0990259) q[1];
x q[2];
rz(-0.17554749) q[3];
sx q[3];
rz(-1.5581589) q[3];
sx q[3];
rz(2.4059084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.79409838) q[2];
sx q[2];
rz(-1.8811052) q[2];
sx q[2];
rz(-2.9079672) q[2];
rz(-2.1485093) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(2.5036507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96930209) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(3.0517975) q[0];
rz(-0.88109294) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(0.18009137) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8343617) q[0];
sx q[0];
rz(-2.0198856) q[0];
sx q[0];
rz(2.2063072) q[0];
rz(-pi) q[1];
rz(2.1744556) q[2];
sx q[2];
rz(-1.4146311) q[2];
sx q[2];
rz(2.9784163) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1485968) q[1];
sx q[1];
rz(-2.8088048) q[1];
sx q[1];
rz(0.51698835) q[1];
rz(0.018499231) q[3];
sx q[3];
rz(-2.0634994) q[3];
sx q[3];
rz(0.25445709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.2580516) q[2];
sx q[2];
rz(-1.5774612) q[2];
sx q[2];
rz(1.1759261) q[2];
rz(-0.073143395) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(0.49889645) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.183855) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(1.4181597) q[0];
rz(-1.9765967) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(0.040239008) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38567625) q[0];
sx q[0];
rz(-1.3323116) q[0];
sx q[0];
rz(-0.083906108) q[0];
rz(-pi) q[1];
x q[1];
rz(0.046182403) q[2];
sx q[2];
rz(-1.1650411) q[2];
sx q[2];
rz(1.6398167) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8602627) q[1];
sx q[1];
rz(-1.7107692) q[1];
sx q[1];
rz(-0.94957385) q[1];
rz(-pi) q[2];
rz(1.2651029) q[3];
sx q[3];
rz(-1.8743268) q[3];
sx q[3];
rz(-2.2895253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0030901) q[2];
sx q[2];
rz(-2.3562727) q[2];
sx q[2];
rz(-0.17383943) q[2];
rz(1.7447757) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(0.85913908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.023507) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(1.7370976) q[0];
rz(1.2069758) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(1.5054024) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0907785) q[0];
sx q[0];
rz(-0.53508004) q[0];
sx q[0];
rz(-0.65061609) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82377388) q[2];
sx q[2];
rz(-2.0588377) q[2];
sx q[2];
rz(1.6549695) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9951524) q[1];
sx q[1];
rz(-2.4447828) q[1];
sx q[1];
rz(3.0909096) q[1];
rz(-0.32780111) q[3];
sx q[3];
rz(-0.9947239) q[3];
sx q[3];
rz(-2.2097261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3796842) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(-2.9910679) q[2];
rz(1.5395509) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6983011) q[0];
sx q[0];
rz(-1.1346096) q[0];
sx q[0];
rz(0.64569965) q[0];
rz(2.6422016) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(1.8766778) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0255233) q[0];
sx q[0];
rz(-0.86206305) q[0];
sx q[0];
rz(0.52225964) q[0];
rz(0.55391295) q[2];
sx q[2];
rz(-1.950112) q[2];
sx q[2];
rz(-1.5253138) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9727271) q[1];
sx q[1];
rz(-1.4251514) q[1];
sx q[1];
rz(-1.5345854) q[1];
rz(-pi) q[2];
rz(0.42322741) q[3];
sx q[3];
rz(-1.755135) q[3];
sx q[3];
rz(-1.1545899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.634793) q[2];
sx q[2];
rz(-2.1719077) q[2];
sx q[2];
rz(-0.52465049) q[2];
rz(0.55650416) q[3];
sx q[3];
rz(-1.5126712) q[3];
sx q[3];
rz(0.35486737) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6181347) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(-2.8961704) q[0];
rz(-0.55631176) q[1];
sx q[1];
rz(-1.5309155) q[1];
sx q[1];
rz(1.4642749) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7510371) q[0];
sx q[0];
rz(-1.3532191) q[0];
sx q[0];
rz(1.4575973) q[0];
x q[1];
rz(2.9173031) q[2];
sx q[2];
rz(-0.99594342) q[2];
sx q[2];
rz(1.526051) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97800868) q[1];
sx q[1];
rz(-2.2469236) q[1];
sx q[1];
rz(-0.34983695) q[1];
rz(3.0562923) q[3];
sx q[3];
rz(-1.7461667) q[3];
sx q[3];
rz(-0.264197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8993373) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(-2.396092) q[2];
rz(1.4238822) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2765008) q[0];
sx q[0];
rz(-1.6456974) q[0];
sx q[0];
rz(0.67281848) q[0];
rz(-1.8854234) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(-0.71184288) q[2];
sx q[2];
rz(-1.1887475) q[2];
sx q[2];
rz(2.9954994) q[2];
rz(-2.5916354) q[3];
sx q[3];
rz(-2.4933542) q[3];
sx q[3];
rz(-1.8174432) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
