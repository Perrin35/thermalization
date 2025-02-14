OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.13088432) q[0];
sx q[0];
rz(-0.67357981) q[0];
sx q[0];
rz(2.3186865) q[0];
rz(1.4505439) q[1];
sx q[1];
rz(-0.6757285) q[1];
sx q[1];
rz(-2.9339209) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0595221) q[0];
sx q[0];
rz(-2.2769109) q[0];
sx q[0];
rz(2.6816899) q[0];
x q[1];
rz(-1.0124341) q[2];
sx q[2];
rz(-1.2744546) q[2];
sx q[2];
rz(1.7425962) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.365578) q[1];
sx q[1];
rz(-1.6322337) q[1];
sx q[1];
rz(0.12428026) q[1];
rz(-pi) q[2];
rz(-1.9413596) q[3];
sx q[3];
rz(-1.5041122) q[3];
sx q[3];
rz(-2.5197864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9660008) q[2];
sx q[2];
rz(-1.5660183) q[2];
sx q[2];
rz(-1.9123745) q[2];
rz(1.2980596) q[3];
sx q[3];
rz(-1.7652054) q[3];
sx q[3];
rz(-1.1840597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9640279) q[0];
sx q[0];
rz(-0.25968817) q[0];
sx q[0];
rz(2.2143256) q[0];
rz(2.941653) q[1];
sx q[1];
rz(-1.4727458) q[1];
sx q[1];
rz(2.1520069) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4662113) q[0];
sx q[0];
rz(-0.37801925) q[0];
sx q[0];
rz(1.2751352) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0485905) q[2];
sx q[2];
rz(-0.43605294) q[2];
sx q[2];
rz(2.2047037) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4633219) q[1];
sx q[1];
rz(-0.090677977) q[1];
sx q[1];
rz(2.1342127) q[1];
rz(-pi) q[2];
rz(0.58122509) q[3];
sx q[3];
rz(-0.94454256) q[3];
sx q[3];
rz(-0.083409781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2135311) q[2];
sx q[2];
rz(-1.7240883) q[2];
sx q[2];
rz(-2.7878063) q[2];
rz(0.95156041) q[3];
sx q[3];
rz(-0.74717251) q[3];
sx q[3];
rz(-2.814754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3554409) q[0];
sx q[0];
rz(-0.94796258) q[0];
sx q[0];
rz(-1.0107262) q[0];
rz(-3.082869) q[1];
sx q[1];
rz(-1.6207691) q[1];
sx q[1];
rz(-2.7322863) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3407281) q[0];
sx q[0];
rz(-2.0052064) q[0];
sx q[0];
rz(-0.77462642) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1597761) q[2];
sx q[2];
rz(-1.7226698) q[2];
sx q[2];
rz(2.0645666) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61534897) q[1];
sx q[1];
rz(-1.0305291) q[1];
sx q[1];
rz(-1.3413642) q[1];
rz(-pi) q[2];
rz(-2.6326551) q[3];
sx q[3];
rz(-1.5889611) q[3];
sx q[3];
rz(3.09336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.056444082) q[2];
sx q[2];
rz(-1.2480382) q[2];
sx q[2];
rz(-0.15963456) q[2];
rz(1.2498445) q[3];
sx q[3];
rz(-1.4188473) q[3];
sx q[3];
rz(-1.0861402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1551664) q[0];
sx q[0];
rz(-2.7041589) q[0];
sx q[0];
rz(2.0470108) q[0];
rz(1.9704341) q[1];
sx q[1];
rz(-2.0234225) q[1];
sx q[1];
rz(-1.0135244) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4552559) q[0];
sx q[0];
rz(-1.5323973) q[0];
sx q[0];
rz(3.0976035) q[0];
x q[1];
rz(-3.1379478) q[2];
sx q[2];
rz(-1.0763775) q[2];
sx q[2];
rz(1.0629176) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.62134939) q[1];
sx q[1];
rz(-1.3101556) q[1];
sx q[1];
rz(1.6127178) q[1];
x q[2];
rz(-1.3427686) q[3];
sx q[3];
rz(-2.6963391) q[3];
sx q[3];
rz(2.5955615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1363498) q[2];
sx q[2];
rz(-1.760027) q[2];
sx q[2];
rz(-1.8278149) q[2];
rz(2.2037196) q[3];
sx q[3];
rz(-1.8504986) q[3];
sx q[3];
rz(-0.098793678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90792847) q[0];
sx q[0];
rz(-1.6133244) q[0];
sx q[0];
rz(-2.0966356) q[0];
rz(0.16432556) q[1];
sx q[1];
rz(-1.3427837) q[1];
sx q[1];
rz(-0.16990653) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92238802) q[0];
sx q[0];
rz(-1.4996075) q[0];
sx q[0];
rz(1.7573331) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3218574) q[2];
sx q[2];
rz(-0.46870527) q[2];
sx q[2];
rz(3.1196032) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.53105866) q[1];
sx q[1];
rz(-1.8137168) q[1];
sx q[1];
rz(0.90621913) q[1];
rz(-pi) q[2];
rz(1.4395797) q[3];
sx q[3];
rz(-0.84892143) q[3];
sx q[3];
rz(-0.038824507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56875151) q[2];
sx q[2];
rz(-2.4390287) q[2];
sx q[2];
rz(1.0315726) q[2];
rz(-2.0555563) q[3];
sx q[3];
rz(-0.7998172) q[3];
sx q[3];
rz(-1.731855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80424911) q[0];
sx q[0];
rz(-2.0954837) q[0];
sx q[0];
rz(1.7403437) q[0];
rz(2.7119472) q[1];
sx q[1];
rz(-0.88637561) q[1];
sx q[1];
rz(-1.3053798) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84112924) q[0];
sx q[0];
rz(-0.80781898) q[0];
sx q[0];
rz(-0.43231583) q[0];
rz(-pi) q[1];
rz(3.1370509) q[2];
sx q[2];
rz(-1.1318739) q[2];
sx q[2];
rz(-2.2786841) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4652047) q[1];
sx q[1];
rz(-0.25034764) q[1];
sx q[1];
rz(-0.12488229) q[1];
rz(0.15302739) q[3];
sx q[3];
rz(-1.8108441) q[3];
sx q[3];
rz(1.9757063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1244916) q[2];
sx q[2];
rz(-1.2064826) q[2];
sx q[2];
rz(0.98480946) q[2];
rz(1.5752327) q[3];
sx q[3];
rz(-1.5896348) q[3];
sx q[3];
rz(-2.8563833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.259909) q[0];
sx q[0];
rz(-0.30170983) q[0];
sx q[0];
rz(1.3551706) q[0];
rz(0.024070865) q[1];
sx q[1];
rz(-1.5211952) q[1];
sx q[1];
rz(2.7640061) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1771072) q[0];
sx q[0];
rz(-1.4146574) q[0];
sx q[0];
rz(2.4859758) q[0];
rz(-pi) q[1];
rz(0.88084282) q[2];
sx q[2];
rz(-0.27327785) q[2];
sx q[2];
rz(2.5246594) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3096034) q[1];
sx q[1];
rz(-1.1909232) q[1];
sx q[1];
rz(-1.989945) q[1];
rz(-pi) q[2];
rz(2.2144775) q[3];
sx q[3];
rz(-0.52442951) q[3];
sx q[3];
rz(2.7853109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5923803) q[2];
sx q[2];
rz(-0.33752957) q[2];
sx q[2];
rz(-0.40360061) q[2];
rz(0.18925439) q[3];
sx q[3];
rz(-2.0892102) q[3];
sx q[3];
rz(2.6846867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52983487) q[0];
sx q[0];
rz(-0.54946041) q[0];
sx q[0];
rz(-0.31103617) q[0];
rz(-0.95651904) q[1];
sx q[1];
rz(-1.6638959) q[1];
sx q[1];
rz(-0.32435736) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8984025) q[0];
sx q[0];
rz(-1.0694188) q[0];
sx q[0];
rz(1.7883975) q[0];
rz(3.0359984) q[2];
sx q[2];
rz(-0.51266951) q[2];
sx q[2];
rz(-1.6216506) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.16087577) q[1];
sx q[1];
rz(-2.389713) q[1];
sx q[1];
rz(-0.94161011) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89357425) q[3];
sx q[3];
rz(-1.8704318) q[3];
sx q[3];
rz(-1.1232337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4057464) q[2];
sx q[2];
rz(-0.79068557) q[2];
sx q[2];
rz(2.9131367) q[2];
rz(-2.6089846) q[3];
sx q[3];
rz(-2.3753128) q[3];
sx q[3];
rz(-0.84958357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.547895) q[0];
sx q[0];
rz(-0.51012796) q[0];
sx q[0];
rz(-1.8713895) q[0];
rz(-0.58865976) q[1];
sx q[1];
rz(-1.9344067) q[1];
sx q[1];
rz(0.38280815) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13148334) q[0];
sx q[0];
rz(-1.3416222) q[0];
sx q[0];
rz(3.1160627) q[0];
x q[1];
rz(-0.50254681) q[2];
sx q[2];
rz(-0.48641962) q[2];
sx q[2];
rz(-0.18196276) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7748357) q[1];
sx q[1];
rz(-2.2765571) q[1];
sx q[1];
rz(-1.7362167) q[1];
rz(2.0476663) q[3];
sx q[3];
rz(-2.5989957) q[3];
sx q[3];
rz(0.93520704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1560912) q[2];
sx q[2];
rz(-1.5727377) q[2];
sx q[2];
rz(0.4001948) q[2];
rz(-2.8943446) q[3];
sx q[3];
rz(-1.6797545) q[3];
sx q[3];
rz(-2.7483773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8409214) q[0];
sx q[0];
rz(-1.8074169) q[0];
sx q[0];
rz(2.1260496) q[0];
rz(-2.4726942) q[1];
sx q[1];
rz(-1.4543507) q[1];
sx q[1];
rz(1.6355754) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4168303) q[0];
sx q[0];
rz(-0.93092954) q[0];
sx q[0];
rz(-0.93389966) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4898275) q[2];
sx q[2];
rz(-0.6547857) q[2];
sx q[2];
rz(-2.1534065) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2336894) q[1];
sx q[1];
rz(-1.091106) q[1];
sx q[1];
rz(3.0976899) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2478836) q[3];
sx q[3];
rz(-2.0314616) q[3];
sx q[3];
rz(3.0627444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.16613913) q[2];
sx q[2];
rz(-2.0659122) q[2];
sx q[2];
rz(1.5258741) q[2];
rz(-0.29159355) q[3];
sx q[3];
rz(-0.44277954) q[3];
sx q[3];
rz(2.7899138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6780554) q[0];
sx q[0];
rz(-0.96374496) q[0];
sx q[0];
rz(-2.5441334) q[0];
rz(-0.28868227) q[1];
sx q[1];
rz(-2.3434227) q[1];
sx q[1];
rz(-3.1176288) q[1];
rz(-0.63715061) q[2];
sx q[2];
rz(-0.49497866) q[2];
sx q[2];
rz(-0.53866932) q[2];
rz(1.3194094) q[3];
sx q[3];
rz(-2.0181927) q[3];
sx q[3];
rz(-0.041139091) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
