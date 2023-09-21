OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.55943638) q[0];
sx q[0];
rz(-2.5506033) q[0];
sx q[0];
rz(-0.58340573) q[0];
rz(-0.18435873) q[1];
sx q[1];
rz(-2.1579722) q[1];
sx q[1];
rz(0.89259994) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4757724) q[0];
sx q[0];
rz(-1.2736397) q[0];
sx q[0];
rz(-2.1509403) q[0];
rz(-pi) q[1];
rz(1.0640261) q[2];
sx q[2];
rz(-2.1760586) q[2];
sx q[2];
rz(-1.93047) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0776699) q[1];
sx q[1];
rz(-2.3706672) q[1];
sx q[1];
rz(-2.0248807) q[1];
rz(-pi) q[2];
rz(-0.75099545) q[3];
sx q[3];
rz(-1.5392116) q[3];
sx q[3];
rz(-0.75814523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2261752) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(-1.9809451) q[2];
rz(0.21696572) q[3];
sx q[3];
rz(-0.522885) q[3];
sx q[3];
rz(-2.0863566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0579257) q[0];
sx q[0];
rz(-0.96494976) q[0];
sx q[0];
rz(2.5657186) q[0];
rz(1.2469762) q[1];
sx q[1];
rz(-1.8449102) q[1];
sx q[1];
rz(-1.1670246) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.089433) q[0];
sx q[0];
rz(-2.8773327) q[0];
sx q[0];
rz(-1.3454076) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1115233) q[2];
sx q[2];
rz(-1.7833372) q[2];
sx q[2];
rz(-0.18268798) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6845219) q[1];
sx q[1];
rz(-1.5091389) q[1];
sx q[1];
rz(2.3350299) q[1];
rz(-pi) q[2];
rz(-0.59252177) q[3];
sx q[3];
rz(-2.3791109) q[3];
sx q[3];
rz(-0.034686397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.033826753) q[2];
sx q[2];
rz(-1.1788538) q[2];
sx q[2];
rz(2.9272184) q[2];
rz(-0.073444627) q[3];
sx q[3];
rz(-2.6918604) q[3];
sx q[3];
rz(-2.8607821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6630994) q[0];
sx q[0];
rz(-2.5254624) q[0];
sx q[0];
rz(-1.2269155) q[0];
rz(2.7413209) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(2.1267557) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0600216) q[0];
sx q[0];
rz(-2.2131753) q[0];
sx q[0];
rz(-1.6409671) q[0];
x q[1];
rz(-1.048676) q[2];
sx q[2];
rz(-2.6819957) q[2];
sx q[2];
rz(-2.2082579) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.2753607) q[1];
sx q[1];
rz(-2.1664157) q[1];
sx q[1];
rz(-3.0443405) q[1];
x q[2];
rz(-1.8869927) q[3];
sx q[3];
rz(-2.4305775) q[3];
sx q[3];
rz(-0.3871813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0456475) q[2];
sx q[2];
rz(-2.084705) q[2];
sx q[2];
rz(2.2423559) q[2];
rz(-0.69747654) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(-0.60788679) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65524453) q[0];
sx q[0];
rz(-1.0826033) q[0];
sx q[0];
rz(-0.43193257) q[0];
rz(-0.63255429) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(-2.5057709) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0974554) q[0];
sx q[0];
rz(-2.6090528) q[0];
sx q[0];
rz(1.2116648) q[0];
rz(-pi) q[1];
rz(-0.16480883) q[2];
sx q[2];
rz(-2.2909082) q[2];
sx q[2];
rz(-2.1317496) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7446049) q[1];
sx q[1];
rz(-2.3481391) q[1];
sx q[1];
rz(2.8810487) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28663978) q[3];
sx q[3];
rz(-1.7754103) q[3];
sx q[3];
rz(-0.74433792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2137961) q[2];
sx q[2];
rz(-0.14363229) q[2];
sx q[2];
rz(-0.68871838) q[2];
rz(0.33411807) q[3];
sx q[3];
rz(-1.9709316) q[3];
sx q[3];
rz(-0.14373246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14389811) q[0];
sx q[0];
rz(-2.4425638) q[0];
sx q[0];
rz(1.0744263) q[0];
rz(0.74514666) q[1];
sx q[1];
rz(-1.4829758) q[1];
sx q[1];
rz(2.863046) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2715831) q[0];
sx q[0];
rz(-0.9013114) q[0];
sx q[0];
rz(2.5894126) q[0];
rz(-pi) q[1];
rz(-0.61168806) q[2];
sx q[2];
rz(-1.4809161) q[2];
sx q[2];
rz(-0.22345605) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4651471) q[1];
sx q[1];
rz(-2.4009581) q[1];
sx q[1];
rz(-2.1315104) q[1];
x q[2];
rz(2.1513125) q[3];
sx q[3];
rz(-2.2970082) q[3];
sx q[3];
rz(-2.5648404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6986065) q[2];
sx q[2];
rz(-1.7437982) q[2];
sx q[2];
rz(2.8473575) q[2];
rz(3.0596628) q[3];
sx q[3];
rz(-2.6223845) q[3];
sx q[3];
rz(-3.1053655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.0329523) q[0];
sx q[0];
rz(-2.2886798) q[0];
sx q[0];
rz(3.1325353) q[0];
rz(-2.5065705) q[1];
sx q[1];
rz(-0.68990866) q[1];
sx q[1];
rz(0.10805282) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2436284) q[0];
sx q[0];
rz(-1.122323) q[0];
sx q[0];
rz(1.2178221) q[0];
rz(-pi) q[1];
rz(2.8073505) q[2];
sx q[2];
rz(-0.40194449) q[2];
sx q[2];
rz(-2.3964756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6335771) q[1];
sx q[1];
rz(-1.0883691) q[1];
sx q[1];
rz(2.0111994) q[1];
rz(-pi) q[2];
rz(2.1170656) q[3];
sx q[3];
rz(-1.3427991) q[3];
sx q[3];
rz(-2.8834164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.99284995) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(1.8704869) q[2];
rz(3.0631915) q[3];
sx q[3];
rz(-1.4784808) q[3];
sx q[3];
rz(1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25061297) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(-0.65761956) q[0];
rz(1.3972067) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(0.89362842) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4511787) q[0];
sx q[0];
rz(-0.72890857) q[0];
sx q[0];
rz(-2.6364987) q[0];
rz(-2.0486084) q[2];
sx q[2];
rz(-1.3078948) q[2];
sx q[2];
rz(-2.773657) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0221755) q[1];
sx q[1];
rz(-0.70964538) q[1];
sx q[1];
rz(0.31125734) q[1];
x q[2];
rz(-0.75974792) q[3];
sx q[3];
rz(-1.6413416) q[3];
sx q[3];
rz(2.7618046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.39067337) q[2];
sx q[2];
rz(-2.1942287) q[2];
sx q[2];
rz(2.612109) q[2];
rz(2.6654065) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(-2.7990394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3787518) q[0];
sx q[0];
rz(-1.5126001) q[0];
sx q[0];
rz(-1.09028) q[0];
rz(0.11225637) q[1];
sx q[1];
rz(-1.1039762) q[1];
sx q[1];
rz(-1.1539248) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33015108) q[0];
sx q[0];
rz(-0.92130565) q[0];
sx q[0];
rz(1.5280456) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5904028) q[2];
sx q[2];
rz(-1.4442208) q[2];
sx q[2];
rz(-1.7984496) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.47626074) q[1];
sx q[1];
rz(-1.6325103) q[1];
sx q[1];
rz(2.3803821) q[1];
x q[2];
rz(-0.38903799) q[3];
sx q[3];
rz(-0.93718796) q[3];
sx q[3];
rz(2.6694359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.551679) q[2];
sx q[2];
rz(-0.38107291) q[2];
sx q[2];
rz(-2.7611458) q[2];
rz(-2.0137265) q[3];
sx q[3];
rz(-1.7815855) q[3];
sx q[3];
rz(-1.5765223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34148759) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(2.5019116) q[0];
rz(1.9027963) q[1];
sx q[1];
rz(-1.4150554) q[1];
sx q[1];
rz(-1.9715086) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1146678) q[0];
sx q[0];
rz(-1.2171193) q[0];
sx q[0];
rz(2.7915519) q[0];
x q[1];
rz(1.8413999) q[2];
sx q[2];
rz(-0.81297183) q[2];
sx q[2];
rz(1.1635309) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1994836) q[1];
sx q[1];
rz(-1.0892727) q[1];
sx q[1];
rz(1.7755309) q[1];
x q[2];
rz(1.7438668) q[3];
sx q[3];
rz(-2.462938) q[3];
sx q[3];
rz(0.21733397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3141979) q[2];
sx q[2];
rz(-0.65076995) q[2];
sx q[2];
rz(2.7588552) q[2];
rz(-2.2132204) q[3];
sx q[3];
rz(-1.174077) q[3];
sx q[3];
rz(-2.476957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9160354) q[0];
sx q[0];
rz(-1.6397497) q[0];
sx q[0];
rz(-2.8826707) q[0];
rz(2.4312773) q[1];
sx q[1];
rz(-1.0537035) q[1];
sx q[1];
rz(0.47992596) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5537162) q[0];
sx q[0];
rz(-1.9037316) q[0];
sx q[0];
rz(-2.2399708) q[0];
rz(-pi) q[1];
rz(-0.9058814) q[2];
sx q[2];
rz(-1.3833628) q[2];
sx q[2];
rz(0.71042176) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6595315) q[1];
sx q[1];
rz(-1.2188984) q[1];
sx q[1];
rz(-0.54606502) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94060937) q[3];
sx q[3];
rz(-0.77946957) q[3];
sx q[3];
rz(1.1704695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.84247983) q[2];
sx q[2];
rz(-2.2128236) q[2];
sx q[2];
rz(-1.3170362) q[2];
rz(1.2420098) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(2.4035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.15923545) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(-0.97933979) q[1];
sx q[1];
rz(-2.0420488) q[1];
sx q[1];
rz(2.2534823) q[1];
rz(1.7469035) q[2];
sx q[2];
rz(-2.0460143) q[2];
sx q[2];
rz(2.569414) q[2];
rz(2.8638774) q[3];
sx q[3];
rz(-1.3125827) q[3];
sx q[3];
rz(-1.7607341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
