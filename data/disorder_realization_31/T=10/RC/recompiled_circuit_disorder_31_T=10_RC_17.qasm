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
rz(2.5685413) q[0];
sx q[0];
rz(11.723784) q[0];
rz(-1.0358345) q[1];
sx q[1];
rz(-2.0422715) q[1];
sx q[1];
rz(1.6834747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18544491) q[0];
sx q[0];
rz(-1.3584134) q[0];
sx q[0];
rz(0.74622112) q[0];
x q[1];
rz(-0.027878472) q[2];
sx q[2];
rz(-0.61402245) q[2];
sx q[2];
rz(-0.30464722) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.018651389) q[1];
sx q[1];
rz(-0.39995799) q[1];
sx q[1];
rz(-0.33756983) q[1];
x q[2];
rz(2.5472766) q[3];
sx q[3];
rz(-1.4694957) q[3];
sx q[3];
rz(-0.017410226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.52790102) q[2];
sx q[2];
rz(-2.1353022) q[2];
sx q[2];
rz(0.17949417) q[2];
rz(-1.2256631) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(0.82204449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74801385) q[0];
sx q[0];
rz(-0.8809692) q[0];
sx q[0];
rz(-2.8161312) q[0];
rz(-1.356396) q[1];
sx q[1];
rz(-2.0929095) q[1];
sx q[1];
rz(1.9869841) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5529454) q[0];
sx q[0];
rz(-1.5879022) q[0];
sx q[0];
rz(-3.1228035) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1968368) q[2];
sx q[2];
rz(-1.2465887) q[2];
sx q[2];
rz(2.7446483) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6007538) q[1];
sx q[1];
rz(-0.80889091) q[1];
sx q[1];
rz(-1.6879338) q[1];
rz(1.8335908) q[3];
sx q[3];
rz(-1.3920708) q[3];
sx q[3];
rz(-2.5457515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6894199) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(2.2581805) q[2];
rz(0.47131053) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8283591) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(1.5154243) q[0];
rz(-0.60107636) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(1.0916969) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1611623) q[0];
sx q[0];
rz(-0.83320252) q[0];
sx q[0];
rz(-2.3479793) q[0];
x q[1];
rz(2.0902363) q[2];
sx q[2];
rz(-0.72548496) q[2];
sx q[2];
rz(1.5067593) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8852383) q[1];
sx q[1];
rz(-1.2043722) q[1];
sx q[1];
rz(0.79856915) q[1];
rz(-pi) q[2];
rz(1.5063498) q[3];
sx q[3];
rz(-1.6985053) q[3];
sx q[3];
rz(-0.32999048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(-2.2606405) q[2];
rz(-1.3736003) q[3];
sx q[3];
rz(-1.6146086) q[3];
sx q[3];
rz(-2.1239471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83051935) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(-0.4367035) q[0];
rz(2.9084335) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(2.8312347) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.719602) q[0];
sx q[0];
rz(-1.4388226) q[0];
sx q[0];
rz(2.7601526) q[0];
rz(-pi) q[1];
rz(1.4448302) q[2];
sx q[2];
rz(-2.2519886) q[2];
sx q[2];
rz(-1.3931526) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1102317) q[1];
sx q[1];
rz(-0.35846113) q[1];
sx q[1];
rz(-1.5178174) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1060018) q[3];
sx q[3];
rz(-1.6944052) q[3];
sx q[3];
rz(-1.3328758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13005304) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(-2.0641573) q[2];
rz(0.056190101) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(-1.594054) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9524277) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(-0.24965832) q[0];
rz(-1.5646308) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(0.87019428) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72083005) q[0];
sx q[0];
rz(-1.8846858) q[0];
sx q[0];
rz(1.785196) q[0];
rz(-pi) q[1];
x q[1];
rz(2.444961) q[2];
sx q[2];
rz(-1.4259035) q[2];
sx q[2];
rz(0.87755132) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.40904564) q[1];
sx q[1];
rz(-2.2605719) q[1];
sx q[1];
rz(-1.9802666) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9719475) q[3];
sx q[3];
rz(-1.0153474) q[3];
sx q[3];
rz(-0.58562216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8683118) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(-0.67374054) q[2];
rz(-2.8379748) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(1.822086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34981397) q[0];
sx q[0];
rz(-2.204201) q[0];
sx q[0];
rz(-0.2579903) q[0];
rz(0.42516431) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(-1.4917096) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0598037) q[0];
sx q[0];
rz(-2.7002618) q[0];
sx q[0];
rz(0.23131891) q[0];
rz(-pi) q[1];
rz(0.67955534) q[2];
sx q[2];
rz(-1.2365885) q[2];
sx q[2];
rz(-2.213775) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3867492) q[1];
sx q[1];
rz(-1.0068839) q[1];
sx q[1];
rz(-2.2350603) q[1];
rz(-pi) q[2];
rz(-1.9353988) q[3];
sx q[3];
rz(-1.4758665) q[3];
sx q[3];
rz(-1.7942384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.012718) q[2];
sx q[2];
rz(-0.96367633) q[2];
sx q[2];
rz(3.0498665) q[2];
rz(-0.84364676) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(-2.2475524) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3180852) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(-2.7303625) q[0];
rz(-0.86589083) q[1];
sx q[1];
rz(-0.31232467) q[1];
sx q[1];
rz(-3.1076028) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995178) q[0];
sx q[0];
rz(-0.80054987) q[0];
sx q[0];
rz(2.9151025) q[0];
x q[1];
rz(-2.2074239) q[2];
sx q[2];
rz(-1.120943) q[2];
sx q[2];
rz(0.63461441) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6778292) q[1];
sx q[1];
rz(-0.67968183) q[1];
sx q[1];
rz(-1.6512647) q[1];
rz(-pi) q[2];
rz(2.9339318) q[3];
sx q[3];
rz(-2.5163979) q[3];
sx q[3];
rz(3.0561662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5380481) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(-0.87654385) q[2];
rz(0.34902469) q[3];
sx q[3];
rz(-1.9411496) q[3];
sx q[3];
rz(-0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2440764) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(0.39392719) q[0];
rz(2.774033) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(-1.6961018) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60364265) q[0];
sx q[0];
rz(-1.1374439) q[0];
sx q[0];
rz(-3.135878) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3169575) q[2];
sx q[2];
rz(-0.7728918) q[2];
sx q[2];
rz(3.0658714) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6519424) q[1];
sx q[1];
rz(-1.4133269) q[1];
sx q[1];
rz(0.55649346) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66283488) q[3];
sx q[3];
rz(-2.088306) q[3];
sx q[3];
rz(-0.55596065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8470856) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(-0.40714804) q[2];
rz(-1.6242705) q[3];
sx q[3];
rz(-1.9842691) q[3];
sx q[3];
rz(-0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.3354934) q[0];
sx q[0];
rz(-0.51500106) q[0];
sx q[0];
rz(1.8898213) q[0];
rz(-0.66954008) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(2.8318185) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71000242) q[0];
sx q[0];
rz(-0.54034034) q[0];
sx q[0];
rz(0.24775981) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0094527761) q[2];
sx q[2];
rz(-1.3938892) q[2];
sx q[2];
rz(-2.0131907) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.138315) q[1];
sx q[1];
rz(-1.365108) q[1];
sx q[1];
rz(-3.0823207) q[1];
rz(-pi) q[2];
rz(1.8272607) q[3];
sx q[3];
rz(-2.1520352) q[3];
sx q[3];
rz(-0.47282156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.70242515) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(-1.2072198) q[2];
rz(1.0369982) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(0.65565482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0913775) q[0];
sx q[0];
rz(-1.8176879) q[0];
sx q[0];
rz(-1.2058831) q[0];
rz(-0.58569113) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(-1.6419798) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9046281) q[0];
sx q[0];
rz(-1.6722073) q[0];
sx q[0];
rz(1.8629835) q[0];
rz(-pi) q[1];
rz(1.502938) q[2];
sx q[2];
rz(-0.80791622) q[2];
sx q[2];
rz(-2.9615336) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2227576) q[1];
sx q[1];
rz(-1.7822052) q[1];
sx q[1];
rz(-2.0445776) q[1];
rz(-pi) q[2];
rz(0.0025000574) q[3];
sx q[3];
rz(-1.3979988) q[3];
sx q[3];
rz(-0.97027422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.88400921) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(-0.94669) q[2];
rz(-2.7729014) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(-0.45599109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(0.35836999) q[0];
sx q[0];
rz(-1.1483648) q[0];
sx q[0];
rz(-0.4233465) q[0];
rz(0.070925698) q[1];
sx q[1];
rz(-1.6880886) q[1];
sx q[1];
rz(-0.26500519) q[1];
rz(2.0680188) q[2];
sx q[2];
rz(-1.9855269) q[2];
sx q[2];
rz(-1.2863458) q[2];
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
