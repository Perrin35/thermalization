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
rz(2.8528557) q[0];
sx q[0];
rz(-0.69536916) q[0];
sx q[0];
rz(-2.8737972) q[0];
rz(0.42203045) q[1];
sx q[1];
rz(-2.2208417) q[1];
sx q[1];
rz(1.2738127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0604297) q[0];
sx q[0];
rz(-1.3295242) q[0];
sx q[0];
rz(-3.0953498) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.535475) q[2];
sx q[2];
rz(-2.173758) q[2];
sx q[2];
rz(-0.22113344) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.63034025) q[1];
sx q[1];
rz(-2.4773438) q[1];
sx q[1];
rz(-0.11099191) q[1];
x q[2];
rz(1.8146562) q[3];
sx q[3];
rz(-0.14774665) q[3];
sx q[3];
rz(0.53099957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5219118) q[2];
sx q[2];
rz(-2.5002067) q[2];
sx q[2];
rz(0.042595159) q[2];
rz(-2.8604782) q[3];
sx q[3];
rz(-1.5646076) q[3];
sx q[3];
rz(0.59578305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5113145) q[0];
sx q[0];
rz(-1.1742641) q[0];
sx q[0];
rz(2.0654772) q[0];
rz(1.0307182) q[1];
sx q[1];
rz(-2.0298256) q[1];
sx q[1];
rz(-1.8310742) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1640747) q[0];
sx q[0];
rz(-0.65854544) q[0];
sx q[0];
rz(1.6048628) q[0];
rz(2.7380472) q[2];
sx q[2];
rz(-0.87658823) q[2];
sx q[2];
rz(-2.780404) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.84634631) q[1];
sx q[1];
rz(-0.2057067) q[1];
sx q[1];
rz(0.61726112) q[1];
rz(2.5111273) q[3];
sx q[3];
rz(-0.63172904) q[3];
sx q[3];
rz(-0.28030685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.92404667) q[2];
sx q[2];
rz(-0.7520389) q[2];
sx q[2];
rz(2.1853866) q[2];
rz(1.3077959) q[3];
sx q[3];
rz(-0.81495133) q[3];
sx q[3];
rz(3.0202878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.2995375) q[0];
sx q[0];
rz(-2.6326023) q[0];
sx q[0];
rz(-3.1053542) q[0];
rz(0.80129519) q[1];
sx q[1];
rz(-1.5472629) q[1];
sx q[1];
rz(1.3005728) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0220003) q[0];
sx q[0];
rz(-1.730607) q[0];
sx q[0];
rz(-0.10800604) q[0];
rz(-pi) q[1];
rz(0.97359263) q[2];
sx q[2];
rz(-1.5351211) q[2];
sx q[2];
rz(-3.1247471) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3727468) q[1];
sx q[1];
rz(-2.0018491) q[1];
sx q[1];
rz(0.94668862) q[1];
rz(2.4518029) q[3];
sx q[3];
rz(-1.2827323) q[3];
sx q[3];
rz(2.7369839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3761882) q[2];
sx q[2];
rz(-1.3724962) q[2];
sx q[2];
rz(0.21793951) q[2];
rz(2.229522) q[3];
sx q[3];
rz(-1.4927161) q[3];
sx q[3];
rz(-0.85165858) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4172149) q[0];
sx q[0];
rz(-2.0700924) q[0];
sx q[0];
rz(-2.3260314) q[0];
rz(1.0527481) q[1];
sx q[1];
rz(-2.2595854) q[1];
sx q[1];
rz(1.7900593) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4824351) q[0];
sx q[0];
rz(-2.0391984) q[0];
sx q[0];
rz(1.4841561) q[0];
rz(-pi) q[1];
rz(1.4830681) q[2];
sx q[2];
rz(-2.2854439) q[2];
sx q[2];
rz(1.8217979) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.38511577) q[1];
sx q[1];
rz(-2.6653892) q[1];
sx q[1];
rz(1.6333196) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4249486) q[3];
sx q[3];
rz(-1.9633246) q[3];
sx q[3];
rz(-1.6300622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1286596) q[2];
sx q[2];
rz(-1.9482875) q[2];
sx q[2];
rz(2.7093757) q[2];
rz(0.84248078) q[3];
sx q[3];
rz(-1.0336927) q[3];
sx q[3];
rz(2.1459818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.1603482) q[0];
sx q[0];
rz(-2.6762185) q[0];
sx q[0];
rz(-0.55111849) q[0];
rz(0.61800686) q[1];
sx q[1];
rz(-1.2851241) q[1];
sx q[1];
rz(-2.0726223) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81449303) q[0];
sx q[0];
rz(-2.2393423) q[0];
sx q[0];
rz(1.8446246) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48241291) q[2];
sx q[2];
rz(-1.97792) q[2];
sx q[2];
rz(-0.77821748) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7643508) q[1];
sx q[1];
rz(-2.3016986) q[1];
sx q[1];
rz(-0.16907558) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.426151) q[3];
sx q[3];
rz(-2.0791884) q[3];
sx q[3];
rz(-1.7252462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7908287) q[2];
sx q[2];
rz(-2.5574234) q[2];
sx q[2];
rz(0.9551777) q[2];
rz(2.5480934) q[3];
sx q[3];
rz(-2.1953526) q[3];
sx q[3];
rz(1.745863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-1.3747568) q[0];
sx q[0];
rz(-3.0992442) q[0];
sx q[0];
rz(1.7497077) q[0];
rz(1.998924) q[1];
sx q[1];
rz(-1.7997768) q[1];
sx q[1];
rz(1.3124189) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8439633) q[0];
sx q[0];
rz(-2.2052551) q[0];
sx q[0];
rz(-2.5353801) q[0];
x q[1];
rz(-0.8719123) q[2];
sx q[2];
rz(-1.0191227) q[2];
sx q[2];
rz(2.0786503) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23199546) q[1];
sx q[1];
rz(-1.0568406) q[1];
sx q[1];
rz(1.5088874) q[1];
rz(1.7465318) q[3];
sx q[3];
rz(-0.32009691) q[3];
sx q[3];
rz(-1.8863581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.352508) q[2];
sx q[2];
rz(-0.75560537) q[2];
sx q[2];
rz(1.4136723) q[2];
rz(-0.46755725) q[3];
sx q[3];
rz(-0.65842015) q[3];
sx q[3];
rz(0.55268923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5107875) q[0];
sx q[0];
rz(-0.063022114) q[0];
sx q[0];
rz(2.4600273) q[0];
rz(1.1850146) q[1];
sx q[1];
rz(-0.76528913) q[1];
sx q[1];
rz(-0.78651816) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9174355) q[0];
sx q[0];
rz(-2.6154714) q[0];
sx q[0];
rz(0.71147646) q[0];
x q[1];
rz(-0.94633905) q[2];
sx q[2];
rz(-1.1558487) q[2];
sx q[2];
rz(0.13042658) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9481755) q[1];
sx q[1];
rz(-1.6778291) q[1];
sx q[1];
rz(2.2501037) q[1];
rz(-pi) q[2];
rz(-0.83966484) q[3];
sx q[3];
rz(-1.6014301) q[3];
sx q[3];
rz(-3.125906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.53175348) q[2];
sx q[2];
rz(-2.9015151) q[2];
sx q[2];
rz(-1.5717724) q[2];
rz(1.8543367) q[3];
sx q[3];
rz(-1.5232892) q[3];
sx q[3];
rz(0.94304812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59931961) q[0];
sx q[0];
rz(-0.71745187) q[0];
sx q[0];
rz(-1.188311) q[0];
rz(-0.020847281) q[1];
sx q[1];
rz(-0.7531082) q[1];
sx q[1];
rz(0.59897024) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63510977) q[0];
sx q[0];
rz(-2.0744893) q[0];
sx q[0];
rz(0.11258188) q[0];
x q[1];
rz(-1.125419) q[2];
sx q[2];
rz(-0.76596224) q[2];
sx q[2];
rz(2.9474023) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4074583) q[1];
sx q[1];
rz(-0.86454138) q[1];
sx q[1];
rz(-2.4962673) q[1];
rz(-2.7008301) q[3];
sx q[3];
rz(-0.89019934) q[3];
sx q[3];
rz(1.1632533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.282436) q[2];
sx q[2];
rz(-1.2500637) q[2];
sx q[2];
rz(0.51521987) q[2];
rz(-0.53269261) q[3];
sx q[3];
rz(-2.0566514) q[3];
sx q[3];
rz(-2.2044619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(3.1135947) q[0];
sx q[0];
rz(-0.6830712) q[0];
sx q[0];
rz(-2.7711476) q[0];
rz(1.0796374) q[1];
sx q[1];
rz(-0.41079435) q[1];
sx q[1];
rz(0.058578514) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3315225) q[0];
sx q[0];
rz(-2.5096623) q[0];
sx q[0];
rz(2.8721316) q[0];
rz(-0.04444261) q[2];
sx q[2];
rz(-1.343633) q[2];
sx q[2];
rz(0.73611605) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9875659) q[1];
sx q[1];
rz(-2.3515764) q[1];
sx q[1];
rz(-1.436019) q[1];
rz(-pi) q[2];
rz(0.58038099) q[3];
sx q[3];
rz(-2.2182052) q[3];
sx q[3];
rz(3.1064432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2728682) q[2];
sx q[2];
rz(-2.3385907) q[2];
sx q[2];
rz(2.8105984) q[2];
rz(-0.013414772) q[3];
sx q[3];
rz(-0.9534854) q[3];
sx q[3];
rz(-2.0851871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5633504) q[0];
sx q[0];
rz(-2.712482) q[0];
sx q[0];
rz(0.44813928) q[0];
rz(-3.0478802) q[1];
sx q[1];
rz(-0.30148503) q[1];
sx q[1];
rz(1.2154481) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5543723) q[0];
sx q[0];
rz(-2.0036812) q[0];
sx q[0];
rz(-3.109193) q[0];
rz(-pi) q[1];
rz(-1.1555668) q[2];
sx q[2];
rz(-2.3390963) q[2];
sx q[2];
rz(1.323483) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7772953) q[1];
sx q[1];
rz(-0.60579311) q[1];
sx q[1];
rz(-1.08849) q[1];
x q[2];
rz(-1.2566725) q[3];
sx q[3];
rz(-0.41659714) q[3];
sx q[3];
rz(2.9521717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1756246) q[2];
sx q[2];
rz(-1.5529239) q[2];
sx q[2];
rz(-0.69941163) q[2];
rz(-1.2753963) q[3];
sx q[3];
rz(-1.9857152) q[3];
sx q[3];
rz(-1.2709966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4324343) q[0];
sx q[0];
rz(-2.2812738) q[0];
sx q[0];
rz(-1.9914837) q[0];
rz(-1.3954096) q[1];
sx q[1];
rz(-1.2184873) q[1];
sx q[1];
rz(1.7582735) q[1];
rz(1.4906314) q[2];
sx q[2];
rz(-2.2350428) q[2];
sx q[2];
rz(2.5000917) q[2];
rz(1.4518573) q[3];
sx q[3];
rz(-1.107286) q[3];
sx q[3];
rz(-2.7881691) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
