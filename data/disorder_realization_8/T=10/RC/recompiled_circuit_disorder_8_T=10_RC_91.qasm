OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8060057) q[0];
sx q[0];
rz(-0.94526362) q[0];
sx q[0];
rz(-0.52559108) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(2.2367509) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63431595) q[0];
sx q[0];
rz(-2.0024558) q[0];
sx q[0];
rz(2.1202203) q[0];
rz(0.71360795) q[2];
sx q[2];
rz(-2.8461694) q[2];
sx q[2];
rz(-1.7507391) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.20476725) q[1];
sx q[1];
rz(-1.1751886) q[1];
sx q[1];
rz(0.25524615) q[1];
rz(1.2600793) q[3];
sx q[3];
rz(-1.7508535) q[3];
sx q[3];
rz(2.3678399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.20377542) q[2];
sx q[2];
rz(-1.4062466) q[2];
sx q[2];
rz(0.096244372) q[2];
rz(-2.105666) q[3];
sx q[3];
rz(-0.38714287) q[3];
sx q[3];
rz(-0.15371418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44089833) q[0];
sx q[0];
rz(-2.7504524) q[0];
sx q[0];
rz(-0.76517117) q[0];
rz(1.8493429) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(-0.66295019) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11766079) q[0];
sx q[0];
rz(-1.6370602) q[0];
sx q[0];
rz(-3.0793385) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0775954) q[2];
sx q[2];
rz(-2.4231744) q[2];
sx q[2];
rz(-1.6990627) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.029164974) q[1];
sx q[1];
rz(-1.6532835) q[1];
sx q[1];
rz(-2.6998181) q[1];
rz(-1.7327659) q[3];
sx q[3];
rz(-1.0893981) q[3];
sx q[3];
rz(3.0961406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.092676) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(0.32901397) q[2];
rz(2.4760903) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(1.8238508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5488141) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(2.9192525) q[0];
rz(-1.0173343) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(-2.6229048) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9496574) q[0];
sx q[0];
rz(-1.2149095) q[0];
sx q[0];
rz(-0.8215254) q[0];
x q[1];
rz(-2.8027595) q[2];
sx q[2];
rz(-0.28140861) q[2];
sx q[2];
rz(-1.1178521) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7171214) q[1];
sx q[1];
rz(-2.3470504) q[1];
sx q[1];
rz(-2.8984927) q[1];
rz(-pi) q[2];
rz(3.1268901) q[3];
sx q[3];
rz(-3.0869752) q[3];
sx q[3];
rz(-0.86044776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8609994) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(0.068543531) q[2];
rz(-0.60244256) q[3];
sx q[3];
rz(-2.3790363) q[3];
sx q[3];
rz(-0.13901916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.65072) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(0.17424507) q[0];
rz(-2.6113367) q[1];
sx q[1];
rz(-1.6590051) q[1];
sx q[1];
rz(-2.6285016) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67447353) q[0];
sx q[0];
rz(-1.6498483) q[0];
sx q[0];
rz(-0.62069686) q[0];
x q[1];
rz(0.60855234) q[2];
sx q[2];
rz(-2.2366479) q[2];
sx q[2];
rz(2.7666639) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.40401134) q[1];
sx q[1];
rz(-1.1176795) q[1];
sx q[1];
rz(1.7246507) q[1];
rz(-pi) q[2];
rz(-1.7771878) q[3];
sx q[3];
rz(-1.6009814) q[3];
sx q[3];
rz(-2.9914732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6667368) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(2.9117057) q[2];
rz(-0.41904467) q[3];
sx q[3];
rz(-1.7925526) q[3];
sx q[3];
rz(-2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6722365) q[0];
sx q[0];
rz(-2.3903963) q[0];
sx q[0];
rz(-2.4647734) q[0];
rz(0.49304402) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(0.61606032) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6989243) q[0];
sx q[0];
rz(-1.6232423) q[0];
sx q[0];
rz(-0.035874493) q[0];
x q[1];
rz(-1.6904171) q[2];
sx q[2];
rz(-1.9336485) q[2];
sx q[2];
rz(-2.5809443) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.38356009) q[1];
sx q[1];
rz(-0.57764232) q[1];
sx q[1];
rz(-0.010096117) q[1];
rz(-pi) q[2];
rz(1.945799) q[3];
sx q[3];
rz(-2.4295394) q[3];
sx q[3];
rz(2.860481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.73413509) q[2];
sx q[2];
rz(-2.5871758) q[2];
sx q[2];
rz(0.25536728) q[2];
rz(1.5363961) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(-0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7191294) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(-2.6690924) q[0];
rz(-0.52608144) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(-0.88476673) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5793943) q[0];
sx q[0];
rz(-1.644855) q[0];
sx q[0];
rz(-1.1893358) q[0];
rz(-2.9315345) q[2];
sx q[2];
rz(-1.5569599) q[2];
sx q[2];
rz(2.7163598) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51984519) q[1];
sx q[1];
rz(-1.2832844) q[1];
sx q[1];
rz(-1.9985755) q[1];
x q[2];
rz(2.8956036) q[3];
sx q[3];
rz(-1.0503328) q[3];
sx q[3];
rz(0.99508475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.74449599) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(-0.3113783) q[2];
rz(-1.3686251) q[3];
sx q[3];
rz(-0.45752782) q[3];
sx q[3];
rz(-0.5293203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7235274) q[0];
sx q[0];
rz(-2.8630246) q[0];
sx q[0];
rz(-0.061070651) q[0];
rz(-3.1014077) q[1];
sx q[1];
rz(-1.9804852) q[1];
sx q[1];
rz(-2.4087002) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6717364) q[0];
sx q[0];
rz(-2.0977019) q[0];
sx q[0];
rz(-2.0335474) q[0];
rz(-pi) q[1];
rz(-0.47809017) q[2];
sx q[2];
rz(-0.92667898) q[2];
sx q[2];
rz(-0.23020506) q[2];
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
rz(-2.1965346) q[1];
sx q[1];
rz(0.0080607944) q[1];
rz(1.2039127) q[3];
sx q[3];
rz(-3.0266579) q[3];
sx q[3];
rz(-0.68655187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0733033) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(2.627009) q[2];
rz(1.9040646) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(-0.54491836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.085389) q[0];
sx q[0];
rz(-1.6563002) q[0];
sx q[0];
rz(2.432166) q[0];
rz(1.5052694) q[1];
sx q[1];
rz(-2.067833) q[1];
sx q[1];
rz(-2.8628796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3942791) q[0];
sx q[0];
rz(-2.8907667) q[0];
sx q[0];
rz(-2.5670811) q[0];
x q[1];
rz(-1.8144572) q[2];
sx q[2];
rz(-2.116034) q[2];
sx q[2];
rz(1.3512163) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0875138) q[1];
sx q[1];
rz(-1.8808108) q[1];
sx q[1];
rz(-0.22884303) q[1];
rz(-pi) q[2];
rz(1.3474943) q[3];
sx q[3];
rz(-1.4947065) q[3];
sx q[3];
rz(-0.014811024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.30148208) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(3.0855132) q[2];
rz(-2.2864443) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(-0.40518951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99496019) q[0];
sx q[0];
rz(-0.17630795) q[0];
sx q[0];
rz(-2.1726998) q[0];
rz(-0.47337198) q[1];
sx q[1];
rz(-2.3469766) q[1];
sx q[1];
rz(1.1425346) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56209598) q[0];
sx q[0];
rz(-1.8651433) q[0];
sx q[0];
rz(3.0466945) q[0];
rz(-0.70334401) q[2];
sx q[2];
rz(-1.568927) q[2];
sx q[2];
rz(1.4700996) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4559608) q[1];
sx q[1];
rz(-0.5545485) q[1];
sx q[1];
rz(2.901652) q[1];
x q[2];
rz(-2.9386018) q[3];
sx q[3];
rz(-2.2309125) q[3];
sx q[3];
rz(-0.70860329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.896495) q[2];
sx q[2];
rz(-1.921804) q[2];
sx q[2];
rz(-2.1208105) q[2];
rz(2.8178689) q[3];
sx q[3];
rz(-0.75298572) q[3];
sx q[3];
rz(-3.1304205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(-2.1616515) q[0];
sx q[0];
rz(-3.1136944) q[0];
sx q[0];
rz(-2.4401869) q[0];
rz(-2.2258863) q[1];
sx q[1];
rz(-2.1332824) q[1];
sx q[1];
rz(-1.9030301) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4558444) q[0];
sx q[0];
rz(-1.0554753) q[0];
sx q[0];
rz(-0.2483764) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5194703) q[2];
sx q[2];
rz(-0.81846279) q[2];
sx q[2];
rz(0.78879702) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0760127) q[1];
sx q[1];
rz(-0.91846839) q[1];
sx q[1];
rz(-0.14821649) q[1];
x q[2];
rz(0.11390399) q[3];
sx q[3];
rz(-2.1602727) q[3];
sx q[3];
rz(-1.9437499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.68676585) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(-3.0977541) q[2];
rz(1.2010126) q[3];
sx q[3];
rz(-2.4062556) q[3];
sx q[3];
rz(-2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4713521) q[0];
sx q[0];
rz(-0.72605194) q[0];
sx q[0];
rz(-1.3656021) q[0];
rz(3.1162221) q[1];
sx q[1];
rz(-1.3062968) q[1];
sx q[1];
rz(1.2702373) q[1];
rz(-0.30517898) q[2];
sx q[2];
rz(-1.1883493) q[2];
sx q[2];
rz(-2.3408008) q[2];
rz(0.17037114) q[3];
sx q[3];
rz(-2.3352565) q[3];
sx q[3];
rz(-0.37526216) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
