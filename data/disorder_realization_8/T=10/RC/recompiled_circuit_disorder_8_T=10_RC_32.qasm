OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.33558694) q[0];
sx q[0];
rz(-2.196329) q[0];
sx q[0];
rz(0.52559108) q[0];
rz(0.2431915) q[1];
sx q[1];
rz(4.3742124) q[1];
sx q[1];
rz(10.32962) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5358937) q[0];
sx q[0];
rz(-2.4568757) q[0];
sx q[0];
rz(2.2936355) q[0];
rz(-1.3741671) q[2];
sx q[2];
rz(-1.7927205) q[2];
sx q[2];
rz(-2.1264399) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2658087) q[1];
sx q[1];
rz(-1.8059397) q[1];
sx q[1];
rz(1.9782515) q[1];
x q[2];
rz(1.2600793) q[3];
sx q[3];
rz(-1.7508535) q[3];
sx q[3];
rz(2.3678399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20377542) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(-0.096244372) q[2];
rz(1.0359267) q[3];
sx q[3];
rz(-0.38714287) q[3];
sx q[3];
rz(2.9878785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7006943) q[0];
sx q[0];
rz(-2.7504524) q[0];
sx q[0];
rz(0.76517117) q[0];
rz(-1.2922497) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(-0.66295019) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5040341) q[0];
sx q[0];
rz(-3.0507037) q[0];
sx q[0];
rz(-0.81764098) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2234369) q[2];
sx q[2];
rz(-1.895972) q[2];
sx q[2];
rz(2.6174389) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1124277) q[1];
sx q[1];
rz(-1.6532835) q[1];
sx q[1];
rz(-2.6998181) q[1];
x q[2];
rz(0.48682537) q[3];
sx q[3];
rz(-1.7142222) q[3];
sx q[3];
rz(-1.6008582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.092676) q[2];
sx q[2];
rz(-2.0930223) q[2];
sx q[2];
rz(2.8125787) q[2];
rz(0.66550231) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(-1.8238508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5927785) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(2.9192525) q[0];
rz(1.0173343) q[1];
sx q[1];
rz(-2.4203114) q[1];
sx q[1];
rz(0.51868784) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.737239) q[0];
sx q[0];
rz(-2.3271932) q[0];
sx q[0];
rz(1.0712207) q[0];
x q[1];
rz(2.8027595) q[2];
sx q[2];
rz(-2.860184) q[2];
sx q[2];
rz(-1.1178521) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0842523) q[1];
sx q[1];
rz(-2.3358314) q[1];
sx q[1];
rz(1.3303824) q[1];
rz(-pi) q[2];
rz(0.054611562) q[3];
sx q[3];
rz(-1.5715989) q[3];
sx q[3];
rz(2.4459248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8609994) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(-0.068543531) q[2];
rz(-2.5391501) q[3];
sx q[3];
rz(-2.3790363) q[3];
sx q[3];
rz(0.13901916) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.65072) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(2.9673476) q[0];
rz(2.6113367) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(0.51309103) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78631567) q[0];
sx q[0];
rz(-2.516541) q[0];
sx q[0];
rz(0.13537188) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3344343) q[2];
sx q[2];
rz(-2.0370738) q[2];
sx q[2];
rz(1.5392787) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9070223) q[1];
sx q[1];
rz(-1.7090194) q[1];
sx q[1];
rz(-0.45781086) q[1];
rz(-1.7170834) q[3];
sx q[3];
rz(-0.20855599) q[3];
sx q[3];
rz(-1.8640765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6667368) q[2];
sx q[2];
rz(-0.36617827) q[2];
sx q[2];
rz(-2.9117057) q[2];
rz(0.41904467) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(0.45927799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4693562) q[0];
sx q[0];
rz(-0.75119632) q[0];
sx q[0];
rz(2.4647734) q[0];
rz(2.6485486) q[1];
sx q[1];
rz(-0.9489916) q[1];
sx q[1];
rz(-2.5255323) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2992633) q[0];
sx q[0];
rz(-3.0780601) q[0];
sx q[0];
rz(-2.1701943) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30454238) q[2];
sx q[2];
rz(-2.7603622) q[2];
sx q[2];
rz(0.234137) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1787781) q[1];
sx q[1];
rz(-1.5652834) q[1];
sx q[1];
rz(2.5639736) q[1];
x q[2];
rz(1.1957937) q[3];
sx q[3];
rz(-2.4295394) q[3];
sx q[3];
rz(-2.860481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4074576) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(-0.25536728) q[2];
rz(1.6051965) q[3];
sx q[3];
rz(-0.95241773) q[3];
sx q[3];
rz(-0.77409625) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7191294) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(2.6690924) q[0];
rz(0.52608144) q[1];
sx q[1];
rz(-2.9317347) q[1];
sx q[1];
rz(-0.88476673) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5793943) q[0];
sx q[0];
rz(-1.644855) q[0];
sx q[0];
rz(-1.9522569) q[0];
x q[1];
rz(-0.21005819) q[2];
sx q[2];
rz(-1.5569599) q[2];
sx q[2];
rz(0.42523281) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5342418) q[1];
sx q[1];
rz(-2.631175) q[1];
sx q[1];
rz(-0.95153248) q[1];
x q[2];
rz(-0.24598908) q[3];
sx q[3];
rz(-2.0912598) q[3];
sx q[3];
rz(2.1465079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3970967) q[2];
sx q[2];
rz(-0.8049736) q[2];
sx q[2];
rz(0.3113783) q[2];
rz(1.3686251) q[3];
sx q[3];
rz(-0.45752782) q[3];
sx q[3];
rz(0.5293203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7235274) q[0];
sx q[0];
rz(-2.8630246) q[0];
sx q[0];
rz(-3.080522) q[0];
rz(-3.1014077) q[1];
sx q[1];
rz(-1.9804852) q[1];
sx q[1];
rz(0.73289245) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3467348) q[0];
sx q[0];
rz(-1.96694) q[0];
sx q[0];
rz(-0.57647716) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.272846) q[2];
sx q[2];
rz(-1.1940496) q[2];
sx q[2];
rz(-1.6422611) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0098795) q[1];
sx q[1];
rz(-0.94505802) q[1];
sx q[1];
rz(-0.0080607944) q[1];
rz(-pi) q[2];
rz(-1.2039127) q[3];
sx q[3];
rz(-0.11493472) q[3];
sx q[3];
rz(-0.68655187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0733033) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(0.51458365) q[2];
rz(1.9040646) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(3.085389) q[0];
sx q[0];
rz(-1.6563002) q[0];
sx q[0];
rz(2.432166) q[0];
rz(1.6363232) q[1];
sx q[1];
rz(-1.0737597) q[1];
sx q[1];
rz(0.27871305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3364209) q[0];
sx q[0];
rz(-1.780691) q[0];
sx q[0];
rz(-1.7091442) q[0];
x q[1];
rz(2.7630373) q[2];
sx q[2];
rz(-2.5494529) q[2];
sx q[2];
rz(2.2373667) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0540788) q[1];
sx q[1];
rz(-1.8808108) q[1];
sx q[1];
rz(0.22884303) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2392427) q[3];
sx q[3];
rz(-2.9058876) q[3];
sx q[3];
rz(1.2329668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8401106) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(-0.056079496) q[2];
rz(-2.2864443) q[3];
sx q[3];
rz(-2.6931098) q[3];
sx q[3];
rz(-2.7364031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99496019) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(0.96889281) q[0];
rz(0.47337198) q[1];
sx q[1];
rz(-2.3469766) q[1];
sx q[1];
rz(1.999058) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0363077) q[0];
sx q[0];
rz(-1.4799911) q[0];
sx q[0];
rz(1.8663976) q[0];
x q[1];
rz(1.5683453) q[2];
sx q[2];
rz(-0.8674538) q[2];
sx q[2];
rz(0.0991115) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0513623) q[1];
sx q[1];
rz(-1.6962595) q[1];
sx q[1];
rz(-0.54162234) q[1];
rz(-1.3167131) q[3];
sx q[3];
rz(-2.4554606) q[3];
sx q[3];
rz(-2.7568267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.896495) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(1.0207821) q[2];
rz(0.3237237) q[3];
sx q[3];
rz(-2.3886069) q[3];
sx q[3];
rz(-3.1304205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97994119) q[0];
sx q[0];
rz(-3.1136944) q[0];
sx q[0];
rz(2.4401869) q[0];
rz(-0.91570634) q[1];
sx q[1];
rz(-2.1332824) q[1];
sx q[1];
rz(1.9030301) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0093875) q[0];
sx q[0];
rz(-1.3552249) q[0];
sx q[0];
rz(2.0995887) q[0];
rz(-0.054758666) q[2];
sx q[2];
rz(-0.75373947) q[2];
sx q[2];
rz(-2.4278305) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.065579942) q[1];
sx q[1];
rz(-0.91846839) q[1];
sx q[1];
rz(-0.14821649) q[1];
rz(0.978312) q[3];
sx q[3];
rz(-1.6654135) q[3];
sx q[3];
rz(2.832151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4548268) q[2];
sx q[2];
rz(-1.0964311) q[2];
sx q[2];
rz(-3.0977541) q[2];
rz(-1.94058) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(-1.003585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4713521) q[0];
sx q[0];
rz(-2.4155407) q[0];
sx q[0];
rz(1.7759905) q[0];
rz(-3.1162221) q[1];
sx q[1];
rz(-1.8352958) q[1];
sx q[1];
rz(-1.8713554) q[1];
rz(2.2123443) q[2];
sx q[2];
rz(-0.48454787) q[2];
sx q[2];
rz(1.5018644) q[2];
rz(-1.745789) q[3];
sx q[3];
rz(-2.3621515) q[3];
sx q[3];
rz(3.0099517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
