OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(-0.57305133) q[0];
sx q[0];
rz(-2.2990062) q[0];
rz(-1.0358345) q[1];
sx q[1];
rz(-2.0422715) q[1];
sx q[1];
rz(1.6834747) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1928756) q[0];
sx q[0];
rz(-0.84515453) q[0];
sx q[0];
rz(1.8564419) q[0];
rz(-pi) q[1];
rz(0.61383944) q[2];
sx q[2];
rz(-1.5547353) q[2];
sx q[2];
rz(-1.2889372) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7588501) q[1];
sx q[1];
rz(-1.9470125) q[1];
sx q[1];
rz(1.43169) q[1];
rz(-pi) q[2];
rz(0.17957844) q[3];
sx q[3];
rz(-0.60185963) q[3];
sx q[3];
rz(1.7367401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52790102) q[2];
sx q[2];
rz(-2.1353022) q[2];
sx q[2];
rz(2.9620985) q[2];
rz(1.2256631) q[3];
sx q[3];
rz(-1.7951199) q[3];
sx q[3];
rz(-2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.74801385) q[0];
sx q[0];
rz(-0.8809692) q[0];
sx q[0];
rz(2.8161312) q[0];
rz(-1.356396) q[1];
sx q[1];
rz(-2.0929095) q[1];
sx q[1];
rz(1.9869841) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5529454) q[0];
sx q[0];
rz(-1.5879022) q[0];
sx q[0];
rz(3.1228035) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7484659) q[2];
sx q[2];
rz(-2.1596585) q[2];
sx q[2];
rz(1.4002422) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3719912) q[1];
sx q[1];
rz(-2.3725315) q[1];
sx q[1];
rz(-0.12188697) q[1];
x q[2];
rz(-1.3080018) q[3];
sx q[3];
rz(-1.3920708) q[3];
sx q[3];
rz(-2.5457515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6894199) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(0.88341218) q[2];
rz(0.47131053) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(2.3538891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.8283591) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(1.5154243) q[0];
rz(-2.5405163) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(2.0498958) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1900345) q[0];
sx q[0];
rz(-1.0150195) q[0];
sx q[0];
rz(-0.65727289) q[0];
rz(-1.0513564) q[2];
sx q[2];
rz(-0.72548496) q[2];
sx q[2];
rz(1.5067593) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.162902) q[1];
sx q[1];
rz(-2.2802417) q[1];
sx q[1];
rz(-0.49180007) q[1];
x q[2];
rz(3.0136209) q[3];
sx q[3];
rz(-1.5068753) q[3];
sx q[3];
rz(-1.909006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(2.2606405) q[2];
rz(1.7679924) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(-1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83051935) q[0];
sx q[0];
rz(-1.3922465) q[0];
sx q[0];
rz(-2.7048892) q[0];
rz(-2.9084335) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(-2.8312347) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.045517) q[0];
sx q[0];
rz(-1.1928416) q[0];
sx q[0];
rz(1.7128574) q[0];
x q[1];
rz(-0.15375806) q[2];
sx q[2];
rz(-2.4506844) q[2];
sx q[2];
rz(1.5916057) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0536641) q[1];
sx q[1];
rz(-1.9287319) q[1];
sx q[1];
rz(3.1217561) q[1];
rz(-3.1060018) q[3];
sx q[3];
rz(-1.4471874) q[3];
sx q[3];
rz(1.3328758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0115396) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(2.0641573) q[2];
rz(-0.056190101) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(-1.594054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9524277) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(2.8919343) q[0];
rz(-1.5769618) q[1];
sx q[1];
rz(-2.3639634) q[1];
sx q[1];
rz(0.87019428) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4207626) q[0];
sx q[0];
rz(-1.8846858) q[0];
sx q[0];
rz(1.785196) q[0];
rz(-2.444961) q[2];
sx q[2];
rz(-1.7156892) q[2];
sx q[2];
rz(0.87755132) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.732547) q[1];
sx q[1];
rz(-0.88102075) q[1];
sx q[1];
rz(1.1613261) q[1];
rz(-pi) q[2];
rz(-1.8364041) q[3];
sx q[3];
rz(-0.57816539) q[3];
sx q[3];
rz(-0.8997013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8683118) q[2];
sx q[2];
rz(-1.3262649) q[2];
sx q[2];
rz(-0.67374054) q[2];
rz(-2.8379748) q[3];
sx q[3];
rz(-1.2250591) q[3];
sx q[3];
rz(-1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7917787) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(0.2579903) q[0];
rz(0.42516431) q[1];
sx q[1];
rz(-0.95562569) q[1];
sx q[1];
rz(-1.649883) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.314635) q[0];
sx q[0];
rz(-1.1420113) q[0];
sx q[0];
rz(1.6786806) q[0];
rz(-pi) q[1];
rz(2.6368124) q[2];
sx q[2];
rz(-0.74540388) q[2];
sx q[2];
rz(0.2573075) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7266453) q[1];
sx q[1];
rz(-2.2989095) q[1];
sx q[1];
rz(2.3689518) q[1];
rz(-pi) q[2];
rz(-1.8317354) q[3];
sx q[3];
rz(-0.3762227) q[3];
sx q[3];
rz(-0.46686831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1288746) q[2];
sx q[2];
rz(-0.96367633) q[2];
sx q[2];
rz(-3.0498665) q[2];
rz(2.2979459) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(0.89404026) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82350746) q[0];
sx q[0];
rz(-1.894269) q[0];
sx q[0];
rz(-2.7303625) q[0];
rz(2.2757018) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(-0.033989865) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86173979) q[0];
sx q[0];
rz(-1.732677) q[0];
sx q[0];
rz(2.3539761) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88862822) q[2];
sx q[2];
rz(-0.76105984) q[2];
sx q[2];
rz(1.6737446) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6778292) q[1];
sx q[1];
rz(-0.67968183) q[1];
sx q[1];
rz(1.6512647) q[1];
x q[2];
rz(-2.5266685) q[3];
sx q[3];
rz(-1.4498386) q[3];
sx q[3];
rz(1.3161591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6035446) q[2];
sx q[2];
rz(-0.59843439) q[2];
sx q[2];
rz(-0.87654385) q[2];
rz(-2.792568) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-2.9984737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.8975163) q[0];
sx q[0];
rz(-1.4325457) q[0];
sx q[0];
rz(-0.39392719) q[0];
rz(-0.36755964) q[1];
sx q[1];
rz(-1.3840679) q[1];
sx q[1];
rz(-1.4454909) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1720393) q[0];
sx q[0];
rz(-1.5759828) q[0];
sx q[0];
rz(2.0041549) q[0];
rz(-pi) q[1];
rz(2.5567899) q[2];
sx q[2];
rz(-1.0324761) q[2];
sx q[2];
rz(-0.98758299) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8135012) q[1];
sx q[1];
rz(-2.5655167) q[1];
sx q[1];
rz(2.8495795) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1962542) q[3];
sx q[3];
rz(-2.1350386) q[3];
sx q[3];
rz(-2.4953147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8470856) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(0.40714804) q[2];
rz(1.5173222) q[3];
sx q[3];
rz(-1.9842691) q[3];
sx q[3];
rz(-0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3354934) q[0];
sx q[0];
rz(-0.51500106) q[0];
sx q[0];
rz(-1.8898213) q[0];
rz(-2.4720526) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(-2.8318185) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4315902) q[0];
sx q[0];
rz(-0.54034034) q[0];
sx q[0];
rz(-2.8938328) q[0];
rz(-pi) q[1];
rz(1.7477112) q[2];
sx q[2];
rz(-1.5801016) q[2];
sx q[2];
rz(-0.44405802) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.27948353) q[1];
sx q[1];
rz(-0.21394193) q[1];
sx q[1];
rz(1.2941542) q[1];
rz(-pi) q[2];
rz(-2.7731032) q[3];
sx q[3];
rz(-2.512305) q[3];
sx q[3];
rz(0.027241782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4391675) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(1.2072198) q[2];
rz(2.1045945) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(2.4859378) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050215125) q[0];
sx q[0];
rz(-1.8176879) q[0];
sx q[0];
rz(1.9357095) q[0];
rz(-2.5559015) q[1];
sx q[1];
rz(-2.060545) q[1];
sx q[1];
rz(-1.6419798) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4831055) q[0];
sx q[0];
rz(-0.30880901) q[0];
sx q[0];
rz(-1.9103861) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0707804) q[2];
sx q[2];
rz(-2.376308) q[2];
sx q[2];
rz(-2.8634957) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75913945) q[1];
sx q[1];
rz(-2.0331953) q[1];
sx q[1];
rz(2.9049302) q[1];
rz(-pi) q[2];
rz(-0.0025000574) q[3];
sx q[3];
rz(-1.3979988) q[3];
sx q[3];
rz(0.97027422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2575834) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(0.94669) q[2];
rz(2.7729014) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(-2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35836999) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(3.070667) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(0.46438607) q[2];
sx q[2];
rz(-1.1190363) q[2];
sx q[2];
rz(-2.6418532) q[2];
rz(2.0675038) q[3];
sx q[3];
rz(-0.91377331) q[3];
sx q[3];
rz(2.5627315) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];