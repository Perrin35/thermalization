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
rz(4.3972754) q[0];
sx q[0];
rz(9.7527405) q[0];
rz(2.9070931) q[1];
sx q[1];
rz(-0.20107888) q[1];
sx q[1];
rz(0.091436401) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5782195) q[0];
sx q[0];
rz(-1.1735998) q[0];
sx q[0];
rz(-0.63011516) q[0];
rz(-pi) q[1];
rz(1.435906) q[2];
sx q[2];
rz(-1.8276668) q[2];
sx q[2];
rz(-2.169662) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.92802231) q[1];
sx q[1];
rz(-1.808128) q[1];
sx q[1];
rz(3.0568394) q[1];
rz(-pi) q[2];
rz(-2.1655733) q[3];
sx q[3];
rz(-1.3888437) q[3];
sx q[3];
rz(0.44427696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4771007) q[2];
sx q[2];
rz(-2.1680919) q[2];
sx q[2];
rz(2.0155902) q[2];
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
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098410957) q[0];
sx q[0];
rz(-2.4480902) q[0];
sx q[0];
rz(2.4480208) q[0];
rz(2.0454848) q[1];
sx q[1];
rz(-0.98384905) q[1];
sx q[1];
rz(2.9512761) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8961401) q[0];
sx q[0];
rz(-1.3724694) q[0];
sx q[0];
rz(-1.1597) q[0];
rz(-0.61043592) q[2];
sx q[2];
rz(-0.62610859) q[2];
sx q[2];
rz(-0.49952835) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1574993) q[1];
sx q[1];
rz(-1.5484973) q[1];
sx q[1];
rz(-0.84866546) q[1];
rz(0.5204366) q[3];
sx q[3];
rz(-2.3428829) q[3];
sx q[3];
rz(-0.48898104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.50743121) q[2];
sx q[2];
rz(-2.4298411) q[2];
sx q[2];
rz(1.2197536) q[2];
rz(-2.9988585) q[3];
sx q[3];
rz(-2.1295363) q[3];
sx q[3];
rz(-0.90863168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74293566) q[0];
sx q[0];
rz(-2.0516899) q[0];
sx q[0];
rz(0.8202585) q[0];
rz(-2.8495158) q[1];
sx q[1];
rz(-1.074011) q[1];
sx q[1];
rz(1.2480199) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67278636) q[0];
sx q[0];
rz(-0.60657036) q[0];
sx q[0];
rz(0.16288217) q[0];
x q[1];
rz(1.5065932) q[2];
sx q[2];
rz(-1.777613) q[2];
sx q[2];
rz(-1.1484255) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.9995211) q[1];
sx q[1];
rz(-1.6377249) q[1];
sx q[1];
rz(0.66758572) q[1];
rz(-pi) q[2];
rz(-2.8415801) q[3];
sx q[3];
rz(-2.1999151) q[3];
sx q[3];
rz(-0.56738561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.32039207) q[2];
sx q[2];
rz(-2.0911066) q[2];
sx q[2];
rz(-1.2134264) q[2];
rz(-0.16472566) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40959013) q[0];
sx q[0];
rz(-1.2816757) q[0];
sx q[0];
rz(0.18606342) q[0];
rz(0.20446725) q[1];
sx q[1];
rz(-0.47195131) q[1];
sx q[1];
rz(-1.2971372) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0766749) q[0];
sx q[0];
rz(-0.13741048) q[0];
sx q[0];
rz(1.8924367) q[0];
rz(-2.9739431) q[2];
sx q[2];
rz(-1.285458) q[2];
sx q[2];
rz(2.167785) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2343826) q[1];
sx q[1];
rz(-2.130059) q[1];
sx q[1];
rz(2.2940966) q[1];
x q[2];
rz(1.6468871) q[3];
sx q[3];
rz(-0.18008672) q[3];
sx q[3];
rz(2.0538581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.90594784) q[2];
sx q[2];
rz(-2.3489958) q[2];
sx q[2];
rz(2.2223991) q[2];
rz(-0.32133189) q[3];
sx q[3];
rz(-2.0733757) q[3];
sx q[3];
rz(-1.2478158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17386757) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(2.2633973) q[0];
rz(-1.8163266) q[1];
sx q[1];
rz(-1.7293431) q[1];
sx q[1];
rz(-1.426288) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4195376) q[0];
sx q[0];
rz(-2.4800322) q[0];
sx q[0];
rz(-1.4369591) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8833141) q[2];
sx q[2];
rz(-1.7973571) q[2];
sx q[2];
rz(2.8001919) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51902276) q[1];
sx q[1];
rz(-1.6791324) q[1];
sx q[1];
rz(-0.83801724) q[1];
rz(-pi) q[2];
rz(-2.3705195) q[3];
sx q[3];
rz(-1.3372476) q[3];
sx q[3];
rz(1.1119103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2720126) q[2];
sx q[2];
rz(-11*pi/13) q[2];
sx q[2];
rz(3.0997979) q[2];
rz(0.061491866) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(2.7048236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7866661) q[0];
sx q[0];
rz(-0.085952856) q[0];
sx q[0];
rz(2.1110995) q[0];
rz(-0.73973918) q[1];
sx q[1];
rz(-1.5286427) q[1];
sx q[1];
rz(-2.5700263) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33486903) q[0];
sx q[0];
rz(-0.81672251) q[0];
sx q[0];
rz(2.4094894) q[0];
rz(-pi) q[1];
rz(2.2366441) q[2];
sx q[2];
rz(-1.6709423) q[2];
sx q[2];
rz(0.03749321) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.73733444) q[1];
sx q[1];
rz(-1.3306276) q[1];
sx q[1];
rz(2.0391697) q[1];
rz(-0.84047517) q[3];
sx q[3];
rz(-0.90208902) q[3];
sx q[3];
rz(-2.3088282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.17343865) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(2.5773933) q[2];
rz(-3.0155904) q[3];
sx q[3];
rz(-1.6832451) q[3];
sx q[3];
rz(-0.29461598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0512222) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(2.4293161) q[0];
rz(0.5258711) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(2.4760822) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3715079) q[0];
sx q[0];
rz(-1.5718282) q[0];
sx q[0];
rz(1.3285711) q[0];
rz(-0.14851103) q[2];
sx q[2];
rz(-0.87887895) q[2];
sx q[2];
rz(-1.4505381) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1674541) q[1];
sx q[1];
rz(-2.0270837) q[1];
sx q[1];
rz(0.2445226) q[1];
rz(-2.1036759) q[3];
sx q[3];
rz(-0.27927342) q[3];
sx q[3];
rz(2.6168407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.15726382) q[2];
sx q[2];
rz(-0.66005808) q[2];
sx q[2];
rz(-1.4228014) q[2];
rz(-0.11519365) q[3];
sx q[3];
rz(-0.51518232) q[3];
sx q[3];
rz(-0.19259024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0916864) q[0];
sx q[0];
rz(-1.1356857) q[0];
sx q[0];
rz(2.9507622) q[0];
rz(2.514839) q[1];
sx q[1];
rz(-2.125506) q[1];
sx q[1];
rz(-0.33871067) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7539464) q[0];
sx q[0];
rz(-0.37746143) q[0];
sx q[0];
rz(-1.9895372) q[0];
rz(-0.45192265) q[2];
sx q[2];
rz(-1.3616614) q[2];
sx q[2];
rz(-0.6349596) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.415886) q[1];
sx q[1];
rz(-2.6542414) q[1];
sx q[1];
rz(1.332167) q[1];
rz(0.87266463) q[3];
sx q[3];
rz(-0.41655311) q[3];
sx q[3];
rz(1.0636898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64615858) q[2];
sx q[2];
rz(-2.50714) q[2];
sx q[2];
rz(2.4411566) q[2];
rz(-2.2436079) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53196466) q[0];
sx q[0];
rz(-2.6599929) q[0];
sx q[0];
rz(-0.65834808) q[0];
rz(0.61093962) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(-3.0019965) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0848207) q[0];
sx q[0];
rz(-3.0773101) q[0];
sx q[0];
rz(1.2540741) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.082824) q[2];
sx q[2];
rz(-0.81548703) q[2];
sx q[2];
rz(1.8975443) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3176206) q[1];
sx q[1];
rz(-0.52597731) q[1];
sx q[1];
rz(3.0082263) q[1];
rz(-pi) q[2];
rz(3.0890907) q[3];
sx q[3];
rz(-0.92696654) q[3];
sx q[3];
rz(-2.1569463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5247941) q[2];
sx q[2];
rz(-2.1308265) q[2];
sx q[2];
rz(2.810478) q[2];
rz(-2.3838499) q[3];
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
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2963294) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(-2.9560126) q[0];
rz(1.0962076) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(1.6569482) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0508704) q[0];
sx q[0];
rz(-1.8403247) q[0];
sx q[0];
rz(-2.0275293) q[0];
rz(1.2946285) q[2];
sx q[2];
rz(-1.3535) q[2];
sx q[2];
rz(-1.5751788) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1440925) q[1];
sx q[1];
rz(-1.9449241) q[1];
sx q[1];
rz(-3.0348026) q[1];
rz(-3.1366523) q[3];
sx q[3];
rz(-0.86588174) q[3];
sx q[3];
rz(0.73965328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.93402702) q[2];
sx q[2];
rz(-0.88946122) q[2];
sx q[2];
rz(-2.5893842) q[2];
rz(-2.3637555) q[3];
sx q[3];
rz(-2.2838897) q[3];
sx q[3];
rz(-2.623693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8778397) q[0];
sx q[0];
rz(-1.5120266) q[0];
sx q[0];
rz(1.7396447) q[0];
rz(2.9539625) q[1];
sx q[1];
rz(-1.3627121) q[1];
sx q[1];
rz(-0.77308853) q[1];
rz(-1.7170231) q[2];
sx q[2];
rz(-1.29371) q[2];
sx q[2];
rz(0.84809662) q[2];
rz(-2.256176) q[3];
sx q[3];
rz(-2.1366742) q[3];
sx q[3];
rz(-2.4921806) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];