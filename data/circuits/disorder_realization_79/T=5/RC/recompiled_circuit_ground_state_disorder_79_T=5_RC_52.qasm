OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3747342) q[0];
sx q[0];
rz(-2.2280966) q[0];
sx q[0];
rz(-1.7537533) q[0];
rz(-2.9517382) q[1];
sx q[1];
rz(-1.318765) q[1];
sx q[1];
rz(2.8554822) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0510584) q[0];
sx q[0];
rz(-1.5947475) q[0];
sx q[0];
rz(-2.6811203) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6880116) q[2];
sx q[2];
rz(-2.5273537) q[2];
sx q[2];
rz(-0.87488824) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.58264153) q[1];
sx q[1];
rz(-1.0233592) q[1];
sx q[1];
rz(-1.8758586) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0864618) q[3];
sx q[3];
rz(-1.1305439) q[3];
sx q[3];
rz(0.80894473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0949969) q[2];
sx q[2];
rz(-1.6511714) q[2];
sx q[2];
rz(-2.4563834) q[2];
rz(-1.4677216) q[3];
sx q[3];
rz(-2.0592212) q[3];
sx q[3];
rz(-0.26836029) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0005223) q[0];
sx q[0];
rz(-1.5474316) q[0];
sx q[0];
rz(-2.1970314) q[0];
rz(3.0682849) q[1];
sx q[1];
rz(-0.83275515) q[1];
sx q[1];
rz(-1.8836969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99745482) q[0];
sx q[0];
rz(-2.1228482) q[0];
sx q[0];
rz(2.6004418) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5340822) q[2];
sx q[2];
rz(-1.3943496) q[2];
sx q[2];
rz(0.29613972) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.98698178) q[1];
sx q[1];
rz(-1.4263337) q[1];
sx q[1];
rz(2.7491436) q[1];
x q[2];
rz(-3.1044534) q[3];
sx q[3];
rz(-1.5243153) q[3];
sx q[3];
rz(-2.6543736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32626095) q[2];
sx q[2];
rz(-1.3911284) q[2];
sx q[2];
rz(-2.1168671) q[2];
rz(2.4552086) q[3];
sx q[3];
rz(-2.4095583) q[3];
sx q[3];
rz(0.91275233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3150385) q[0];
sx q[0];
rz(-1.2737561) q[0];
sx q[0];
rz(2.3231373) q[0];
rz(1.9326899) q[1];
sx q[1];
rz(-1.5070288) q[1];
sx q[1];
rz(0.85982972) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8743961) q[0];
sx q[0];
rz(-2.1566297) q[0];
sx q[0];
rz(0.39975014) q[0];
x q[1];
rz(-2.4218338) q[2];
sx q[2];
rz(-1.9406332) q[2];
sx q[2];
rz(-1.0074266) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2261994) q[1];
sx q[1];
rz(-0.92297727) q[1];
sx q[1];
rz(-2.1658705) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1016261) q[3];
sx q[3];
rz(-0.2519603) q[3];
sx q[3];
rz(0.89706883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.39765057) q[2];
sx q[2];
rz(-0.33471477) q[2];
sx q[2];
rz(-3.0261377) q[2];
rz(0.3913106) q[3];
sx q[3];
rz(-1.0262998) q[3];
sx q[3];
rz(1.1439884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7855969) q[0];
sx q[0];
rz(-1.1859256) q[0];
sx q[0];
rz(-2.9010229) q[0];
rz(-1.3661512) q[1];
sx q[1];
rz(-1.4931449) q[1];
sx q[1];
rz(1.2496525) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8927596) q[0];
sx q[0];
rz(-0.9380149) q[0];
sx q[0];
rz(0.31676333) q[0];
rz(0.95486887) q[2];
sx q[2];
rz(-1.8661341) q[2];
sx q[2];
rz(-1.5952641) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3981957) q[1];
sx q[1];
rz(-0.92602387) q[1];
sx q[1];
rz(-1.866524) q[1];
rz(-pi) q[2];
rz(-2.4576538) q[3];
sx q[3];
rz(-1.9762282) q[3];
sx q[3];
rz(-2.8059174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8852641) q[2];
sx q[2];
rz(-0.54577959) q[2];
sx q[2];
rz(-1.3679999) q[2];
rz(0.3041501) q[3];
sx q[3];
rz(-1.5757898) q[3];
sx q[3];
rz(-2.4450541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2279385) q[0];
sx q[0];
rz(-2.9253687) q[0];
sx q[0];
rz(-2.6044593) q[0];
rz(0.74983239) q[1];
sx q[1];
rz(-2.0593819) q[1];
sx q[1];
rz(0.73712635) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0219624) q[0];
sx q[0];
rz(-3.0263794) q[0];
sx q[0];
rz(2.553279) q[0];
rz(-pi) q[1];
rz(-1.7456878) q[2];
sx q[2];
rz(-1.6564257) q[2];
sx q[2];
rz(-2.0599287) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3143718) q[1];
sx q[1];
rz(-1.4431568) q[1];
sx q[1];
rz(0.26083827) q[1];
rz(-pi) q[2];
rz(1.6629822) q[3];
sx q[3];
rz(-1.8927729) q[3];
sx q[3];
rz(0.76065242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4445112) q[2];
sx q[2];
rz(-2.7883174) q[2];
sx q[2];
rz(-0.14661655) q[2];
rz(-1.692449) q[3];
sx q[3];
rz(-1.861898) q[3];
sx q[3];
rz(-0.52391887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.523282) q[0];
sx q[0];
rz(-1.016541) q[0];
sx q[0];
rz(1.7200394) q[0];
rz(1.8824185) q[1];
sx q[1];
rz(-0.6764532) q[1];
sx q[1];
rz(0.4037942) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.291639) q[0];
sx q[0];
rz(-1.8516415) q[0];
sx q[0];
rz(2.3333059) q[0];
rz(-2.8844599) q[2];
sx q[2];
rz(-1.7499353) q[2];
sx q[2];
rz(0.15878294) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.48133367) q[1];
sx q[1];
rz(-2.1161882) q[1];
sx q[1];
rz(2.27687) q[1];
x q[2];
rz(-2.4533692) q[3];
sx q[3];
rz(-1.7587708) q[3];
sx q[3];
rz(1.926762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28973618) q[2];
sx q[2];
rz(-1.7278262) q[2];
sx q[2];
rz(-2.9750321) q[2];
rz(1.6377431) q[3];
sx q[3];
rz(-0.47347355) q[3];
sx q[3];
rz(3.04305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.903776) q[0];
sx q[0];
rz(-0.38919583) q[0];
sx q[0];
rz(2.26407) q[0];
rz(0.96218836) q[1];
sx q[1];
rz(-2.018237) q[1];
sx q[1];
rz(-2.7872564) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0838881) q[0];
sx q[0];
rz(-1.5885) q[0];
sx q[0];
rz(0.23772765) q[0];
x q[1];
rz(1.7155859) q[2];
sx q[2];
rz(-1.5939404) q[2];
sx q[2];
rz(-1.0095846) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6207032) q[1];
sx q[1];
rz(-2.1701309) q[1];
sx q[1];
rz(0.66461892) q[1];
rz(-pi) q[2];
rz(-2.6316492) q[3];
sx q[3];
rz(-1.5048319) q[3];
sx q[3];
rz(3.0086781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80498901) q[2];
sx q[2];
rz(-0.73770928) q[2];
sx q[2];
rz(-1.0450411) q[2];
rz(0.40237829) q[3];
sx q[3];
rz(-1.9741524) q[3];
sx q[3];
rz(1.4761402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43110111) q[0];
sx q[0];
rz(-3.0248108) q[0];
sx q[0];
rz(2.2905599) q[0];
rz(-2.7878413) q[1];
sx q[1];
rz(-1.4101135) q[1];
sx q[1];
rz(-2.2869349) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86164323) q[0];
sx q[0];
rz(-0.33726528) q[0];
sx q[0];
rz(0.57798903) q[0];
rz(-pi) q[1];
rz(-0.24979892) q[2];
sx q[2];
rz(-1.7160176) q[2];
sx q[2];
rz(1.9171361) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9766337) q[1];
sx q[1];
rz(-1.7384496) q[1];
sx q[1];
rz(2.555545) q[1];
rz(0.1885957) q[3];
sx q[3];
rz(-2.7321987) q[3];
sx q[3];
rz(-0.90367095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9510368) q[2];
sx q[2];
rz(-0.64728105) q[2];
sx q[2];
rz(2.5212042) q[2];
rz(-0.22647151) q[3];
sx q[3];
rz(-1.3969996) q[3];
sx q[3];
rz(2.5605719) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9613551) q[0];
sx q[0];
rz(-1.9397475) q[0];
sx q[0];
rz(-0.015137976) q[0];
rz(1.0221647) q[1];
sx q[1];
rz(-2.731555) q[1];
sx q[1];
rz(1.2000363) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7040492) q[0];
sx q[0];
rz(-1.3945197) q[0];
sx q[0];
rz(0.86221077) q[0];
x q[1];
rz(1.5115509) q[2];
sx q[2];
rz(-0.97840913) q[2];
sx q[2];
rz(-0.5991583) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4440587) q[1];
sx q[1];
rz(-1.6007533) q[1];
sx q[1];
rz(-1.0590963) q[1];
rz(-0.8127549) q[3];
sx q[3];
rz(-2.4852607) q[3];
sx q[3];
rz(1.9539273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6156561) q[2];
sx q[2];
rz(-2.2286712) q[2];
sx q[2];
rz(2.9021662) q[2];
rz(0.058852363) q[3];
sx q[3];
rz(-2.3299496) q[3];
sx q[3];
rz(2.8656901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9489768) q[0];
sx q[0];
rz(-1.9276351) q[0];
sx q[0];
rz(0.97491997) q[0];
rz(-0.67543593) q[1];
sx q[1];
rz(-0.91921872) q[1];
sx q[1];
rz(-2.1598037) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0115765) q[0];
sx q[0];
rz(-2.1476319) q[0];
sx q[0];
rz(-0.31188282) q[0];
x q[1];
rz(-2.7300225) q[2];
sx q[2];
rz(-1.5082422) q[2];
sx q[2];
rz(2.6175509) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2416545) q[1];
sx q[1];
rz(-1.3119427) q[1];
sx q[1];
rz(-0.83251013) q[1];
rz(-pi) q[2];
rz(-1.124106) q[3];
sx q[3];
rz(-1.5367931) q[3];
sx q[3];
rz(-1.7336577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5741817) q[2];
sx q[2];
rz(-0.97102204) q[2];
sx q[2];
rz(1.9166454) q[2];
rz(2.775906) q[3];
sx q[3];
rz(-0.94958011) q[3];
sx q[3];
rz(-1.1398116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.097261978) q[0];
sx q[0];
rz(-1.8987569) q[0];
sx q[0];
rz(1.4415997) q[0];
rz(3.0615831) q[1];
sx q[1];
rz(-1.7623822) q[1];
sx q[1];
rz(-0.0026127876) q[1];
rz(1.4314797) q[2];
sx q[2];
rz(-2.0124544) q[2];
sx q[2];
rz(-2.6228946) q[2];
rz(1.5936218) q[3];
sx q[3];
rz(-2.0232087) q[3];
sx q[3];
rz(3.051331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
