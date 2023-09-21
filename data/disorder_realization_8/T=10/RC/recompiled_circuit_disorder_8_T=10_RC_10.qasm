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
rz(2.6160016) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(-0.90484172) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5072767) q[0];
sx q[0];
rz(-2.0024558) q[0];
sx q[0];
rz(-1.0213724) q[0];
rz(-2.4279847) q[2];
sx q[2];
rz(-0.2954233) q[2];
sx q[2];
rz(1.7507391) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.875784) q[1];
sx q[1];
rz(-1.8059397) q[1];
sx q[1];
rz(1.1633412) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0337898) q[3];
sx q[3];
rz(-0.35764965) q[3];
sx q[3];
rz(1.8358177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9378172) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(0.096244372) q[2];
rz(-1.0359267) q[3];
sx q[3];
rz(-2.7544498) q[3];
sx q[3];
rz(-0.15371418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7006943) q[0];
sx q[0];
rz(-2.7504524) q[0];
sx q[0];
rz(2.3764215) q[0];
rz(1.2922497) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(0.66295019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0239319) q[0];
sx q[0];
rz(-1.5045325) q[0];
sx q[0];
rz(0.06225417) q[0];
rz(0.91815572) q[2];
sx q[2];
rz(-1.2456206) q[2];
sx q[2];
rz(-2.6174389) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4274802) q[1];
sx q[1];
rz(-0.44891) q[1];
sx q[1];
rz(0.19101363) q[1];
rz(-2.6547673) q[3];
sx q[3];
rz(-1.7142222) q[3];
sx q[3];
rz(1.5407345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.048916653) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(0.32901397) q[2];
rz(0.66550231) q[3];
sx q[3];
rz(-2.9232959) q[3];
sx q[3];
rz(1.8238508) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5488141) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(-2.9192525) q[0];
rz(-2.1242583) q[1];
sx q[1];
rz(-2.4203114) q[1];
sx q[1];
rz(-2.6229048) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.737239) q[0];
sx q[0];
rz(-2.3271932) q[0];
sx q[0];
rz(2.0703719) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33883314) q[2];
sx q[2];
rz(-2.860184) q[2];
sx q[2];
rz(-2.0237405) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0573404) q[1];
sx q[1];
rz(-0.80576128) q[1];
sx q[1];
rz(1.3303824) q[1];
rz(-pi) q[2];
rz(3.1268901) q[3];
sx q[3];
rz(-0.054617453) q[3];
sx q[3];
rz(-2.2811449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8609994) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(-3.0730491) q[2];
rz(-2.5391501) q[3];
sx q[3];
rz(-0.76255637) q[3];
sx q[3];
rz(3.0025735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4908726) q[0];
sx q[0];
rz(-0.77712178) q[0];
sx q[0];
rz(-2.9673476) q[0];
rz(-2.6113367) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(-0.51309103) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4671191) q[0];
sx q[0];
rz(-1.6498483) q[0];
sx q[0];
rz(-2.5208958) q[0];
x q[1];
rz(-0.80715837) q[2];
sx q[2];
rz(-2.0370738) q[2];
sx q[2];
rz(1.5392787) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2345703) q[1];
sx q[1];
rz(-1.7090194) q[1];
sx q[1];
rz(-2.6837818) q[1];
rz(-pi) q[2];
rz(-1.3644049) q[3];
sx q[3];
rz(-1.5406113) q[3];
sx q[3];
rz(-2.9914732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6667368) q[2];
sx q[2];
rz(-0.36617827) q[2];
sx q[2];
rz(-0.22988698) q[2];
rz(-0.41904467) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6722365) q[0];
sx q[0];
rz(-2.3903963) q[0];
sx q[0];
rz(-2.4647734) q[0];
rz(2.6485486) q[1];
sx q[1];
rz(-0.9489916) q[1];
sx q[1];
rz(0.61606032) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2992633) q[0];
sx q[0];
rz(-0.063532524) q[0];
sx q[0];
rz(0.97139831) q[0];
x q[1];
rz(-2.8370503) q[2];
sx q[2];
rz(-2.7603622) q[2];
sx q[2];
rz(-2.9074557) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9628145) q[1];
sx q[1];
rz(-1.5763092) q[1];
sx q[1];
rz(-2.5639736) q[1];
rz(-1.945799) q[3];
sx q[3];
rz(-0.71205322) q[3];
sx q[3];
rz(2.860481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73413509) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(2.8862254) q[2];
rz(1.5363961) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(-0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42246321) q[0];
sx q[0];
rz(-0.96619773) q[0];
sx q[0];
rz(0.47250026) q[0];
rz(-0.52608144) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(-0.88476673) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.103325) q[0];
sx q[0];
rz(-1.9511576) q[0];
sx q[0];
rz(-0.079770712) q[0];
rz(0.066263513) q[2];
sx q[2];
rz(-2.931086) q[2];
sx q[2];
rz(-1.2103684) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51984519) q[1];
sx q[1];
rz(-1.8583082) q[1];
sx q[1];
rz(-1.9985755) q[1];
rz(-2.8956036) q[3];
sx q[3];
rz(-1.0503328) q[3];
sx q[3];
rz(2.1465079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74449599) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(0.3113783) q[2];
rz(-1.7729676) q[3];
sx q[3];
rz(-2.6840648) q[3];
sx q[3];
rz(-0.5293203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41806528) q[0];
sx q[0];
rz(-0.27856809) q[0];
sx q[0];
rz(0.061070651) q[0];
rz(3.1014077) q[1];
sx q[1];
rz(-1.9804852) q[1];
sx q[1];
rz(2.4087002) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6717364) q[0];
sx q[0];
rz(-2.0977019) q[0];
sx q[0];
rz(-2.0335474) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0211208) q[2];
sx q[2];
rz(-2.3602544) q[2];
sx q[2];
rz(2.6598425) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1317132) q[1];
sx q[1];
rz(-0.94505802) q[1];
sx q[1];
rz(0.0080607944) q[1];
x q[2];
rz(1.4634499) q[3];
sx q[3];
rz(-1.5296474) q[3];
sx q[3];
rz(-0.51957182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0733033) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(0.51458365) q[2];
rz(1.2375281) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(-2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056203689) q[0];
sx q[0];
rz(-1.6563002) q[0];
sx q[0];
rz(2.432166) q[0];
rz(-1.5052694) q[1];
sx q[1];
rz(-1.0737597) q[1];
sx q[1];
rz(0.27871305) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80517171) q[0];
sx q[0];
rz(-1.3609017) q[0];
sx q[0];
rz(-1.4324485) q[0];
rz(-0.3785554) q[2];
sx q[2];
rz(-0.59213973) q[2];
sx q[2];
rz(-2.2373667) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41234327) q[1];
sx q[1];
rz(-1.7885498) q[1];
sx q[1];
rz(-1.8885683) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.078019402) q[3];
sx q[3];
rz(-1.3481513) q[3];
sx q[3];
rz(-1.5732461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8401106) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(-0.056079496) q[2];
rz(-0.85514832) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(0.40518951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1466325) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(2.1726998) q[0];
rz(0.47337198) q[1];
sx q[1];
rz(-0.7946161) q[1];
sx q[1];
rz(1.1425346) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5794967) q[0];
sx q[0];
rz(-1.8651433) q[0];
sx q[0];
rz(-3.0466945) q[0];
x q[1];
rz(-0.70334401) q[2];
sx q[2];
rz(-1.568927) q[2];
sx q[2];
rz(1.4700996) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4559608) q[1];
sx q[1];
rz(-0.5545485) q[1];
sx q[1];
rz(-2.901652) q[1];
rz(-pi) q[2];
rz(-1.8248796) q[3];
sx q[3];
rz(-0.68613201) q[3];
sx q[3];
rz(0.38476598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2450976) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(2.1208105) q[2];
rz(0.3237237) q[3];
sx q[3];
rz(-2.3886069) q[3];
sx q[3];
rz(-3.1304205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1616515) q[0];
sx q[0];
rz(-0.027898235) q[0];
sx q[0];
rz(0.7014057) q[0];
rz(2.2258863) q[1];
sx q[1];
rz(-1.0083102) q[1];
sx q[1];
rz(1.2385626) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21047132) q[0];
sx q[0];
rz(-0.56715542) q[0];
sx q[0];
rz(1.9803067) q[0];
x q[1];
rz(3.086834) q[2];
sx q[2];
rz(-2.3878532) q[2];
sx q[2];
rz(2.4278305) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.065579942) q[1];
sx q[1];
rz(-2.2231243) q[1];
sx q[1];
rz(-0.14821649) q[1];
rz(-pi) q[2];
rz(-0.978312) q[3];
sx q[3];
rz(-1.6654135) q[3];
sx q[3];
rz(0.30944165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.68676585) q[2];
sx q[2];
rz(-1.0964311) q[2];
sx q[2];
rz(-3.0977541) q[2];
rz(1.94058) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(1.003585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-0.92924835) q[2];
sx q[2];
rz(-0.48454787) q[2];
sx q[2];
rz(1.5018644) q[2];
rz(-2.3425441) q[3];
sx q[3];
rz(-1.6934762) q[3];
sx q[3];
rz(1.3140524) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
