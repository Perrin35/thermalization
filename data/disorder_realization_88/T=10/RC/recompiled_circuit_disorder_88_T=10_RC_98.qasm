OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.24703439) q[0];
sx q[0];
rz(-1.8859099) q[0];
sx q[0];
rz(0.32796252) q[0];
rz(-0.2344996) q[1];
sx q[1];
rz(3.3426715) q[1];
sx q[1];
rz(9.3333416) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8592005) q[0];
sx q[0];
rz(-0.99635591) q[0];
sx q[0];
rz(1.0919071) q[0];
x q[1];
rz(-0.47317998) q[2];
sx q[2];
rz(-0.28943974) q[2];
sx q[2];
rz(-2.6602886) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4788471) q[1];
sx q[1];
rz(-1.6531684) q[1];
sx q[1];
rz(-1.3326416) q[1];
rz(-pi) q[2];
rz(-2.923008) q[3];
sx q[3];
rz(-0.9871452) q[3];
sx q[3];
rz(-1.0047131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4771007) q[2];
sx q[2];
rz(-0.97350073) q[2];
sx q[2];
rz(1.1260024) q[2];
rz(2.8664355) q[3];
sx q[3];
rz(-2.5313009) q[3];
sx q[3];
rz(0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1317516) q[0];
sx q[0];
rz(-2.4480902) q[0];
sx q[0];
rz(-2.4480208) q[0];
rz(2.0454848) q[1];
sx q[1];
rz(-0.98384905) q[1];
sx q[1];
rz(-0.19031659) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099146518) q[0];
sx q[0];
rz(-0.45396807) q[0];
sx q[0];
rz(1.1048261) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1778189) q[2];
sx q[2];
rz(-2.071638) q[2];
sx q[2];
rz(0.21265342) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9840934) q[1];
sx q[1];
rz(-1.5930953) q[1];
sx q[1];
rz(-0.84866546) q[1];
rz(-0.72782794) q[3];
sx q[3];
rz(-1.2065294) q[3];
sx q[3];
rz(1.4621853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6341614) q[2];
sx q[2];
rz(-2.4298411) q[2];
sx q[2];
rz(-1.921839) q[2];
rz(-2.9988585) q[3];
sx q[3];
rz(-2.1295363) q[3];
sx q[3];
rz(2.232961) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74293566) q[0];
sx q[0];
rz(-2.0516899) q[0];
sx q[0];
rz(0.8202585) q[0];
rz(-0.29207686) q[1];
sx q[1];
rz(-1.074011) q[1];
sx q[1];
rz(1.8935727) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0322198) q[0];
sx q[0];
rz(-1.6633699) q[0];
sx q[0];
rz(-2.5412482) q[0];
rz(1.5065932) q[2];
sx q[2];
rz(-1.3639796) q[2];
sx q[2];
rz(-1.9931672) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4856845) q[1];
sx q[1];
rz(-2.471172) q[1];
sx q[1];
rz(3.0337435) q[1];
x q[2];
rz(-2.8415801) q[3];
sx q[3];
rz(-2.1999151) q[3];
sx q[3];
rz(2.574207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8212006) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(-1.9281663) q[2];
rz(2.976867) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40959013) q[0];
sx q[0];
rz(-1.2816757) q[0];
sx q[0];
rz(-2.9555292) q[0];
rz(-0.20446725) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(1.8444555) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1870677) q[0];
sx q[0];
rz(-1.6141119) q[0];
sx q[0];
rz(1.4403507) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9739431) q[2];
sx q[2];
rz(-1.285458) q[2];
sx q[2];
rz(-2.167785) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0190977) q[1];
sx q[1];
rz(-0.88224733) q[1];
sx q[1];
rz(0.81329878) q[1];
rz(-0.013838776) q[3];
sx q[3];
rz(-1.3912364) q[3];
sx q[3];
rz(1.9765215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2356448) q[2];
sx q[2];
rz(-2.3489958) q[2];
sx q[2];
rz(2.2223991) q[2];
rz(2.8202608) q[3];
sx q[3];
rz(-1.0682169) q[3];
sx q[3];
rz(-1.8937768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17386757) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(-0.87819535) q[0];
rz(1.8163266) q[1];
sx q[1];
rz(-1.4122496) q[1];
sx q[1];
rz(1.7153046) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25709897) q[0];
sx q[0];
rz(-1.4887267) q[0];
sx q[0];
rz(0.9135855) q[0];
rz(-pi) q[1];
rz(-1.2582785) q[2];
sx q[2];
rz(-1.3442355) q[2];
sx q[2];
rz(0.3414008) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9702455) q[1];
sx q[1];
rz(-2.4023224) q[1];
sx q[1];
rz(1.7319748) q[1];
rz(-0.77107314) q[3];
sx q[3];
rz(-1.8043451) q[3];
sx q[3];
rz(-2.0296823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8695801) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(-0.041794725) q[2];
rz(-3.0801008) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(-0.4367691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.3549266) q[0];
sx q[0];
rz(-3.0556398) q[0];
sx q[0];
rz(-2.1110995) q[0];
rz(-0.73973918) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(-0.57156634) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.556658) q[0];
sx q[0];
rz(-0.99781636) q[0];
sx q[0];
rz(-0.95227382) q[0];
rz(-pi) q[1];
rz(-0.90494855) q[2];
sx q[2];
rz(-1.4706503) q[2];
sx q[2];
rz(-0.03749321) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.9532277) q[1];
sx q[1];
rz(-2.0247012) q[1];
sx q[1];
rz(-2.8737349) q[1];
rz(-pi) q[2];
rz(-0.84047517) q[3];
sx q[3];
rz(-0.90208902) q[3];
sx q[3];
rz(-2.3088282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.968154) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(0.56419939) q[2];
rz(-3.0155904) q[3];
sx q[3];
rz(-1.6832451) q[3];
sx q[3];
rz(-0.29461598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-1.0903704) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(-2.4293161) q[0];
rz(-2.6157216) q[1];
sx q[1];
rz(-0.41627517) q[1];
sx q[1];
rz(0.66551048) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8004566) q[0];
sx q[0];
rz(-1.3285713) q[0];
sx q[0];
rz(3.1405297) q[0];
rz(1.3940784) q[2];
sx q[2];
rz(-0.70509796) q[2];
sx q[2];
rz(1.2202028) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4588786) q[1];
sx q[1];
rz(-2.6280118) q[1];
sx q[1];
rz(-2.0290124) q[1];
rz(-pi) q[2];
rz(-1.3286367) q[3];
sx q[3];
rz(-1.7112964) q[3];
sx q[3];
rz(-1.5797918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15726382) q[2];
sx q[2];
rz(-0.66005808) q[2];
sx q[2];
rz(1.4228014) q[2];
rz(-0.11519365) q[3];
sx q[3];
rz(-0.51518232) q[3];
sx q[3];
rz(-0.19259024) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0916864) q[0];
sx q[0];
rz(-1.1356857) q[0];
sx q[0];
rz(0.19083047) q[0];
rz(0.62675369) q[1];
sx q[1];
rz(-2.125506) q[1];
sx q[1];
rz(-2.802882) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3074293) q[0];
sx q[0];
rz(-1.9142262) q[0];
sx q[0];
rz(2.9817581) q[0];
rz(-pi) q[1];
rz(1.3391101) q[2];
sx q[2];
rz(-1.1294239) q[2];
sx q[2];
rz(2.1053134) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0847816) q[1];
sx q[1];
rz(-1.4598795) q[1];
sx q[1];
rz(2.0463498) q[1];
rz(-pi) q[2];
rz(-0.2770822) q[3];
sx q[3];
rz(-1.8859366) q[3];
sx q[3];
rz(-1.3354697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4954341) q[2];
sx q[2];
rz(-2.50714) q[2];
sx q[2];
rz(0.70043606) q[2];
rz(0.8979848) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(-0.17351304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.609628) q[0];
sx q[0];
rz(-2.6599929) q[0];
sx q[0];
rz(-0.65834808) q[0];
rz(0.61093962) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(-3.0019965) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0848207) q[0];
sx q[0];
rz(-3.0773101) q[0];
sx q[0];
rz(-1.8875185) q[0];
x q[1];
rz(0.058768674) q[2];
sx q[2];
rz(-2.3261056) q[2];
sx q[2];
rz(-1.8975443) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.82397205) q[1];
sx q[1];
rz(-2.6156153) q[1];
sx q[1];
rz(3.0082263) q[1];
x q[2];
rz(1.6406052) q[3];
sx q[3];
rz(-0.6456635) q[3];
sx q[3];
rz(-2.0696236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5247941) q[2];
sx q[2];
rz(-2.1308265) q[2];
sx q[2];
rz(-0.33111462) q[2];
rz(0.75774276) q[3];
sx q[3];
rz(-0.38882935) q[3];
sx q[3];
rz(-0.087879114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.2963294) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(2.9560126) q[0];
rz(2.045385) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(1.4846444) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34996841) q[0];
sx q[0];
rz(-1.1317283) q[0];
sx q[0];
rz(-0.29859782) q[0];
x q[1];
rz(0.89016685) q[2];
sx q[2];
rz(-0.3496799) q[2];
sx q[2];
rz(-0.65469757) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6757322) q[1];
sx q[1];
rz(-1.670174) q[1];
sx q[1];
rz(-1.9468716) q[1];
rz(-pi) q[2];
rz(-3.1366523) q[3];
sx q[3];
rz(-0.86588174) q[3];
sx q[3];
rz(0.73965328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2075656) q[2];
sx q[2];
rz(-0.88946122) q[2];
sx q[2];
rz(-2.5893842) q[2];
rz(0.77783716) q[3];
sx q[3];
rz(-2.2838897) q[3];
sx q[3];
rz(-2.623693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8778397) q[0];
sx q[0];
rz(-1.5120266) q[0];
sx q[0];
rz(1.7396447) q[0];
rz(-0.18763018) q[1];
sx q[1];
rz(-1.3627121) q[1];
sx q[1];
rz(-0.77308853) q[1];
rz(-1.7170231) q[2];
sx q[2];
rz(-1.29371) q[2];
sx q[2];
rz(0.84809662) q[2];
rz(2.4545112) q[3];
sx q[3];
rz(-1.0071181) q[3];
sx q[3];
rz(-1.3345171) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
