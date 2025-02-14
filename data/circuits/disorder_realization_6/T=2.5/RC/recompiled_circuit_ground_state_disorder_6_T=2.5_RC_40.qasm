OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4226469) q[0];
sx q[0];
rz(-0.5889686) q[0];
sx q[0];
rz(0.11864057) q[0];
rz(2.150382) q[1];
sx q[1];
rz(-1.041643) q[1];
sx q[1];
rz(0.51377327) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27072752) q[0];
sx q[0];
rz(-2.0692192) q[0];
sx q[0];
rz(-0.99683572) q[0];
rz(3.0209885) q[2];
sx q[2];
rz(-2.6636811) q[2];
sx q[2];
rz(1.8574886) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39756718) q[1];
sx q[1];
rz(-1.7850707) q[1];
sx q[1];
rz(-2.812723) q[1];
x q[2];
rz(-0.051187201) q[3];
sx q[3];
rz(-0.59268206) q[3];
sx q[3];
rz(0.41646233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4673246) q[2];
sx q[2];
rz(-1.6215723) q[2];
sx q[2];
rz(2.8094214) q[2];
rz(-0.25003555) q[3];
sx q[3];
rz(-1.9146405) q[3];
sx q[3];
rz(-1.8488319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2501204) q[0];
sx q[0];
rz(-1.8880867) q[0];
sx q[0];
rz(2.6943595) q[0];
rz(-2.1121292) q[1];
sx q[1];
rz(-0.58164683) q[1];
sx q[1];
rz(-3.0404125) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6529236) q[0];
sx q[0];
rz(-0.61572585) q[0];
sx q[0];
rz(-2.4020345) q[0];
rz(-pi) q[1];
rz(-0.20002983) q[2];
sx q[2];
rz(-1.0009655) q[2];
sx q[2];
rz(2.9404158) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9969719) q[1];
sx q[1];
rz(-2.640814) q[1];
sx q[1];
rz(-2.5624496) q[1];
x q[2];
rz(-2.2096735) q[3];
sx q[3];
rz(-1.6611691) q[3];
sx q[3];
rz(2.3109009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0589361) q[2];
sx q[2];
rz(-2.3857748) q[2];
sx q[2];
rz(-1.5400881) q[2];
rz(0.61795175) q[3];
sx q[3];
rz(-2.5585585) q[3];
sx q[3];
rz(-1.9790953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59442941) q[0];
sx q[0];
rz(-1.0048486) q[0];
sx q[0];
rz(0.51602236) q[0];
rz(-2.0891321) q[1];
sx q[1];
rz(-2.2155589) q[1];
sx q[1];
rz(-1.1722391) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2207551) q[0];
sx q[0];
rz(-1.1467548) q[0];
sx q[0];
rz(1.5829257) q[0];
x q[1];
rz(1.486023) q[2];
sx q[2];
rz(-0.71226487) q[2];
sx q[2];
rz(-0.9916457) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.13048417) q[1];
sx q[1];
rz(-1.3846894) q[1];
sx q[1];
rz(-2.0687201) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1276305) q[3];
sx q[3];
rz(-0.80028557) q[3];
sx q[3];
rz(1.518702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8664794) q[2];
sx q[2];
rz(-1.2629513) q[2];
sx q[2];
rz(-2.4647554) q[2];
rz(-0.46946851) q[3];
sx q[3];
rz(-0.89478409) q[3];
sx q[3];
rz(1.2137871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6770099) q[0];
sx q[0];
rz(-1.0053758) q[0];
sx q[0];
rz(-1.9418035) q[0];
rz(-2.5489573) q[1];
sx q[1];
rz(-1.0721782) q[1];
sx q[1];
rz(-1.4592272) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98226794) q[0];
sx q[0];
rz(-1.766681) q[0];
sx q[0];
rz(-0.59024337) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8267385) q[2];
sx q[2];
rz(-2.3530966) q[2];
sx q[2];
rz(2.7389604) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.063733405) q[1];
sx q[1];
rz(-0.95653906) q[1];
sx q[1];
rz(2.8251404) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68379559) q[3];
sx q[3];
rz(-1.8701564) q[3];
sx q[3];
rz(1.9366858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49372855) q[2];
sx q[2];
rz(-2.04546) q[2];
sx q[2];
rz(0.68191051) q[2];
rz(0.41334263) q[3];
sx q[3];
rz(-1.5324493) q[3];
sx q[3];
rz(-1.1558862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2979564) q[0];
sx q[0];
rz(-1.4763259) q[0];
sx q[0];
rz(2.8635039) q[0];
rz(0.9043215) q[1];
sx q[1];
rz(-1.6878637) q[1];
sx q[1];
rz(0.98731891) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8191746) q[0];
sx q[0];
rz(-1.8083982) q[0];
sx q[0];
rz(-0.95652076) q[0];
rz(-3.004641) q[2];
sx q[2];
rz(-1.5986575) q[2];
sx q[2];
rz(2.7951023) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1121171) q[1];
sx q[1];
rz(-1.5095081) q[1];
sx q[1];
rz(-1.960683) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0342399) q[3];
sx q[3];
rz(-1.8066386) q[3];
sx q[3];
rz(2.4135116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.014293369) q[2];
sx q[2];
rz(-2.1743446) q[2];
sx q[2];
rz(-2.3822752) q[2];
rz(1.6589818) q[3];
sx q[3];
rz(-1.9448152) q[3];
sx q[3];
rz(-0.35879859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9229014) q[0];
sx q[0];
rz(-2.3221115) q[0];
sx q[0];
rz(-2.4500093) q[0];
rz(0.85451952) q[1];
sx q[1];
rz(-0.71047345) q[1];
sx q[1];
rz(-0.42133322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3940659) q[0];
sx q[0];
rz(-2.0025064) q[0];
sx q[0];
rz(-0.075448087) q[0];
x q[1];
rz(2.8745804) q[2];
sx q[2];
rz(-1.8521554) q[2];
sx q[2];
rz(-0.54131484) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5798963) q[1];
sx q[1];
rz(-1.4868127) q[1];
sx q[1];
rz(-0.3616228) q[1];
rz(-pi) q[2];
rz(2.5365127) q[3];
sx q[3];
rz(-1.7972094) q[3];
sx q[3];
rz(-0.070158557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.93115807) q[2];
sx q[2];
rz(-1.4782108) q[2];
sx q[2];
rz(-2.9388536) q[2];
rz(-1.5431131) q[3];
sx q[3];
rz(-0.32035443) q[3];
sx q[3];
rz(-1.4690442) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28432524) q[0];
sx q[0];
rz(-1.7131282) q[0];
sx q[0];
rz(2.7427234) q[0];
rz(1.9786037) q[1];
sx q[1];
rz(-2.3154924) q[1];
sx q[1];
rz(2.861048) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5359176) q[0];
sx q[0];
rz(-1.1725764) q[0];
sx q[0];
rz(-1.1541973) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0899827) q[2];
sx q[2];
rz(-0.74724846) q[2];
sx q[2];
rz(2.0488536) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.043634) q[1];
sx q[1];
rz(-2.6065738) q[1];
sx q[1];
rz(2.19997) q[1];
rz(-pi) q[2];
rz(-2.1637502) q[3];
sx q[3];
rz(-0.8675608) q[3];
sx q[3];
rz(3.0271526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3205388) q[2];
sx q[2];
rz(-2.1619449) q[2];
sx q[2];
rz(1.4873571) q[2];
rz(0.94333831) q[3];
sx q[3];
rz(-2.2127559) q[3];
sx q[3];
rz(-1.1613891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5281552) q[0];
sx q[0];
rz(-1.4035839) q[0];
sx q[0];
rz(0.50233895) q[0];
rz(0.54652864) q[1];
sx q[1];
rz(-1.8859452) q[1];
sx q[1];
rz(-2.6796403) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34794688) q[0];
sx q[0];
rz(-1.0805305) q[0];
sx q[0];
rz(-2.8099174) q[0];
x q[1];
rz(-0.047872825) q[2];
sx q[2];
rz(-0.92000735) q[2];
sx q[2];
rz(-1.4818459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5344055) q[1];
sx q[1];
rz(-2.7173373) q[1];
sx q[1];
rz(-1.3035745) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1199014) q[3];
sx q[3];
rz(-2.6110161) q[3];
sx q[3];
rz(-2.2013045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.97303566) q[2];
sx q[2];
rz(-2.904197) q[2];
sx q[2];
rz(-0.17042223) q[2];
rz(2.6863344) q[3];
sx q[3];
rz(-1.596343) q[3];
sx q[3];
rz(-2.4236603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97487226) q[0];
sx q[0];
rz(-1.6396739) q[0];
sx q[0];
rz(0.18044743) q[0];
rz(2.6249053) q[1];
sx q[1];
rz(-1.5770715) q[1];
sx q[1];
rz(2.7076941) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.949373) q[0];
sx q[0];
rz(-1.1447971) q[0];
sx q[0];
rz(2.1550687) q[0];
x q[1];
rz(-1.1230311) q[2];
sx q[2];
rz(-1.8862714) q[2];
sx q[2];
rz(1.4184784) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5916768) q[1];
sx q[1];
rz(-2.2001451) q[1];
sx q[1];
rz(2.9931328) q[1];
x q[2];
rz(-0.29124864) q[3];
sx q[3];
rz(-0.90121239) q[3];
sx q[3];
rz(-0.88676363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.58069289) q[2];
sx q[2];
rz(-2.7036724) q[2];
sx q[2];
rz(1.3194579) q[2];
rz(-0.27160078) q[3];
sx q[3];
rz(-1.8235794) q[3];
sx q[3];
rz(1.1118838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6962947) q[0];
sx q[0];
rz(-2.2309208) q[0];
sx q[0];
rz(1.4060422) q[0];
rz(1.7156853) q[1];
sx q[1];
rz(-1.1049263) q[1];
sx q[1];
rz(-1.8941194) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6794327) q[0];
sx q[0];
rz(-1.832167) q[0];
sx q[0];
rz(1.0699468) q[0];
x q[1];
rz(-1.6883759) q[2];
sx q[2];
rz(-1.6248068) q[2];
sx q[2];
rz(0.62979311) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6599802) q[1];
sx q[1];
rz(-1.2289189) q[1];
sx q[1];
rz(-1.5274164) q[1];
x q[2];
rz(-2.5188451) q[3];
sx q[3];
rz(-1.0120262) q[3];
sx q[3];
rz(1.8468685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4282816) q[2];
sx q[2];
rz(-2.9094628) q[2];
sx q[2];
rz(3.0408119) q[2];
rz(-0.0095602592) q[3];
sx q[3];
rz(-1.5567501) q[3];
sx q[3];
rz(-1.6077707) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99209256) q[0];
sx q[0];
rz(-1.7763573) q[0];
sx q[0];
rz(-1.1471164) q[0];
rz(-0.59303444) q[1];
sx q[1];
rz(-1.4321764) q[1];
sx q[1];
rz(2.166116) q[1];
rz(-3.0708457) q[2];
sx q[2];
rz(-1.9683217) q[2];
sx q[2];
rz(-0.34839658) q[2];
rz(1.5884052) q[3];
sx q[3];
rz(-1.3497769) q[3];
sx q[3];
rz(-2.6365437) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
