OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.61252874) q[0];
sx q[0];
rz(3.7747369) q[0];
sx q[0];
rz(9.7882191) q[0];
rz(2.5545622) q[1];
sx q[1];
rz(-0.5293923) q[1];
sx q[1];
rz(-1.5076948) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76985659) q[0];
sx q[0];
rz(-1.6210067) q[0];
sx q[0];
rz(1.4771853) q[0];
rz(-pi) q[1];
rz(0.71552566) q[2];
sx q[2];
rz(-1.9811077) q[2];
sx q[2];
rz(-2.7844051) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.22836338) q[1];
sx q[1];
rz(-2.1268561) q[1];
sx q[1];
rz(-0.25198752) q[1];
rz(-2.5938183) q[3];
sx q[3];
rz(-2.2507077) q[3];
sx q[3];
rz(2.821673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2734566) q[2];
sx q[2];
rz(-3.0507705) q[2];
sx q[2];
rz(2.691213) q[2];
rz(0.069933794) q[3];
sx q[3];
rz(-2.0262597) q[3];
sx q[3];
rz(-0.44225881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0792585) q[0];
sx q[0];
rz(-0.82759768) q[0];
sx q[0];
rz(2.8424679) q[0];
rz(-0.66080719) q[1];
sx q[1];
rz(-1.1102763) q[1];
sx q[1];
rz(0.92346907) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7851334) q[0];
sx q[0];
rz(-1.74455) q[0];
sx q[0];
rz(-1.2453402) q[0];
x q[1];
rz(-1.7910576) q[2];
sx q[2];
rz(-1.8008322) q[2];
sx q[2];
rz(-2.736614) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0638173) q[1];
sx q[1];
rz(-2.7690921) q[1];
sx q[1];
rz(-0.17613303) q[1];
rz(-pi) q[2];
rz(0.98387444) q[3];
sx q[3];
rz(-1.3345342) q[3];
sx q[3];
rz(-0.4224225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6580761) q[2];
sx q[2];
rz(-0.4008171) q[2];
sx q[2];
rz(0.083902396) q[2];
rz(0.99533254) q[3];
sx q[3];
rz(-1.1570802) q[3];
sx q[3];
rz(1.6307433) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8453269) q[0];
sx q[0];
rz(-1.2355329) q[0];
sx q[0];
rz(2.1790225) q[0];
rz(2.8341017) q[1];
sx q[1];
rz(-0.88491416) q[1];
sx q[1];
rz(0.21893758) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013857582) q[0];
sx q[0];
rz(-1.9163736) q[0];
sx q[0];
rz(1.4968977) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23425183) q[2];
sx q[2];
rz(-2.2356389) q[2];
sx q[2];
rz(-1.6981036) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1477222) q[1];
sx q[1];
rz(-1.1712554) q[1];
sx q[1];
rz(-0.0030722105) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2454065) q[3];
sx q[3];
rz(-2.0749669) q[3];
sx q[3];
rz(-1.4880796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4585877) q[2];
sx q[2];
rz(-1.747749) q[2];
sx q[2];
rz(-1.5031507) q[2];
rz(2.891244) q[3];
sx q[3];
rz(-0.72385794) q[3];
sx q[3];
rz(2.943128) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4918764) q[0];
sx q[0];
rz(-0.32418495) q[0];
sx q[0];
rz(-2.8247996) q[0];
rz(3.1253539) q[1];
sx q[1];
rz(-2.3490678) q[1];
sx q[1];
rz(0.24831717) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1386928) q[0];
sx q[0];
rz(-1.0637245) q[0];
sx q[0];
rz(2.4207628) q[0];
x q[1];
rz(-1.4961719) q[2];
sx q[2];
rz(-1.1047102) q[2];
sx q[2];
rz(-1.3206732) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.59245306) q[1];
sx q[1];
rz(-2.7902979) q[1];
sx q[1];
rz(-0.37254803) q[1];
rz(-pi) q[2];
rz(-0.19222741) q[3];
sx q[3];
rz(-1.3557443) q[3];
sx q[3];
rz(0.59753093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2636286) q[2];
sx q[2];
rz(-3.0402277) q[2];
sx q[2];
rz(-1.144009) q[2];
rz(2.3794203) q[3];
sx q[3];
rz(-0.52505571) q[3];
sx q[3];
rz(-2.150056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85631973) q[0];
sx q[0];
rz(-0.090949051) q[0];
sx q[0];
rz(-2.0891304) q[0];
rz(0.4902803) q[1];
sx q[1];
rz(-1.0888381) q[1];
sx q[1];
rz(3.0163684) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019524272) q[0];
sx q[0];
rz(-2.6112289) q[0];
sx q[0];
rz(-2.2838462) q[0];
rz(-pi) q[1];
rz(2.7020821) q[2];
sx q[2];
rz(-2.548303) q[2];
sx q[2];
rz(0.35542929) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0908806) q[1];
sx q[1];
rz(-1.4987603) q[1];
sx q[1];
rz(-1.3204095) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85876561) q[3];
sx q[3];
rz(-1.4031583) q[3];
sx q[3];
rz(-1.110838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6584955) q[2];
sx q[2];
rz(-0.53762895) q[2];
sx q[2];
rz(1.5888692) q[2];
rz(-0.056564864) q[3];
sx q[3];
rz(-1.5141234) q[3];
sx q[3];
rz(-2.8642843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5618943) q[0];
sx q[0];
rz(-2.6772006) q[0];
sx q[0];
rz(2.7138746) q[0];
rz(-2.600421) q[1];
sx q[1];
rz(-2.7310557) q[1];
sx q[1];
rz(-1.3637095) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.246733) q[0];
sx q[0];
rz(-1.680516) q[0];
sx q[0];
rz(2.1734957) q[0];
rz(-0.12597106) q[2];
sx q[2];
rz(-0.4301542) q[2];
sx q[2];
rz(-0.90829721) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8783166) q[1];
sx q[1];
rz(-1.8748302) q[1];
sx q[1];
rz(2.2869907) q[1];
rz(1.1579944) q[3];
sx q[3];
rz(-1.5808453) q[3];
sx q[3];
rz(-1.5665986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.042818729) q[2];
sx q[2];
rz(-2.6218178) q[2];
sx q[2];
rz(0.40936142) q[2];
rz(-2.9218819) q[3];
sx q[3];
rz(-1.5922056) q[3];
sx q[3];
rz(1.6519914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4396502) q[0];
sx q[0];
rz(-0.42124978) q[0];
sx q[0];
rz(0.14081328) q[0];
rz(2.6048062) q[1];
sx q[1];
rz(-1.355492) q[1];
sx q[1];
rz(0.37667325) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55467419) q[0];
sx q[0];
rz(-1.7307052) q[0];
sx q[0];
rz(-1.317509) q[0];
x q[1];
rz(-1.9697625) q[2];
sx q[2];
rz(-2.2062183) q[2];
sx q[2];
rz(2.2081326) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46493775) q[1];
sx q[1];
rz(-2.3466517) q[1];
sx q[1];
rz(1.388768) q[1];
x q[2];
rz(-2.7262732) q[3];
sx q[3];
rz(-1.7139846) q[3];
sx q[3];
rz(0.48986942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.10304) q[2];
sx q[2];
rz(-1.975375) q[2];
sx q[2];
rz(-0.80123365) q[2];
rz(-1.0771105) q[3];
sx q[3];
rz(-0.25682768) q[3];
sx q[3];
rz(-2.9991155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.899026) q[0];
sx q[0];
rz(-2.9809256) q[0];
sx q[0];
rz(2.2005079) q[0];
rz(2.2451027) q[1];
sx q[1];
rz(-1.4184378) q[1];
sx q[1];
rz(-2.3268708) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0705256) q[0];
sx q[0];
rz(-1.1909426) q[0];
sx q[0];
rz(-1.1761355) q[0];
rz(-pi) q[1];
rz(-0.17085393) q[2];
sx q[2];
rz(-1.5207982) q[2];
sx q[2];
rz(1.6544261) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.22413177) q[1];
sx q[1];
rz(-0.96012174) q[1];
sx q[1];
rz(2.7236205) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1188468) q[3];
sx q[3];
rz(-2.4999385) q[3];
sx q[3];
rz(1.1430571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5251289) q[2];
sx q[2];
rz(-0.63462555) q[2];
sx q[2];
rz(-2.5647707) q[2];
rz(0.47354627) q[3];
sx q[3];
rz(-1.9927931) q[3];
sx q[3];
rz(3.098251) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7035377) q[0];
sx q[0];
rz(-0.75551581) q[0];
sx q[0];
rz(3.0536861) q[0];
rz(1.5880623) q[1];
sx q[1];
rz(-0.5843662) q[1];
sx q[1];
rz(-2.8839674) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15730298) q[0];
sx q[0];
rz(-2.7098795) q[0];
sx q[0];
rz(1.2832416) q[0];
x q[1];
rz(-0.070775971) q[2];
sx q[2];
rz(-0.94460605) q[2];
sx q[2];
rz(-0.51838057) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.96366072) q[1];
sx q[1];
rz(-2.6448634) q[1];
sx q[1];
rz(-0.50951783) q[1];
rz(0.99434234) q[3];
sx q[3];
rz(-1.0691158) q[3];
sx q[3];
rz(-2.8714406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.075500175) q[2];
sx q[2];
rz(-1.4614212) q[2];
sx q[2];
rz(3.063805) q[2];
rz(1.0545688) q[3];
sx q[3];
rz(-2.4013897) q[3];
sx q[3];
rz(-2.4906702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8973549) q[0];
sx q[0];
rz(-0.033441823) q[0];
sx q[0];
rz(0.43575409) q[0];
rz(-3.1089973) q[1];
sx q[1];
rz(-0.73740149) q[1];
sx q[1];
rz(-1.3479412) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4530058) q[0];
sx q[0];
rz(-2.3212715) q[0];
sx q[0];
rz(1.561101) q[0];
rz(-pi) q[1];
rz(-1.8135494) q[2];
sx q[2];
rz(-1.0565783) q[2];
sx q[2];
rz(0.16788858) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4159591) q[1];
sx q[1];
rz(-1.6068649) q[1];
sx q[1];
rz(-3.1276324) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8043933) q[3];
sx q[3];
rz(-1.8775619) q[3];
sx q[3];
rz(-1.5535958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9255479) q[2];
sx q[2];
rz(-2.1074769) q[2];
sx q[2];
rz(-0.75630581) q[2];
rz(-3.0594399) q[3];
sx q[3];
rz(-2.7492456) q[3];
sx q[3];
rz(0.69123417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3638851) q[0];
sx q[0];
rz(-2.4019699) q[0];
sx q[0];
rz(-0.85516047) q[0];
rz(-2.0289519) q[1];
sx q[1];
rz(-1.1323977) q[1];
sx q[1];
rz(-0.36066537) q[1];
rz(-1.821221) q[2];
sx q[2];
rz(-1.741521) q[2];
sx q[2];
rz(2.3157816) q[2];
rz(1.011547) q[3];
sx q[3];
rz(-1.9725298) q[3];
sx q[3];
rz(-1.6172668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
