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
rz(-0.32514206) q[0];
sx q[0];
rz(-2.4976375) q[0];
sx q[0];
rz(2.08026) q[0];
rz(-pi) q[1];
rz(-2.0775665) q[2];
sx q[2];
rz(-0.96553409) q[2];
sx q[2];
rz(1.93047) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0776699) q[1];
sx q[1];
rz(-0.7709255) q[1];
sx q[1];
rz(-2.0248807) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0953232) q[3];
sx q[3];
rz(-2.390063) q[3];
sx q[3];
rz(-0.7788333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2261752) q[2];
sx q[2];
rz(-1.6025275) q[2];
sx q[2];
rz(1.1606476) q[2];
rz(-0.21696572) q[3];
sx q[3];
rz(-0.522885) q[3];
sx q[3];
rz(2.0863566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-0.57587409) q[0];
rz(1.8946164) q[1];
sx q[1];
rz(-1.2966825) q[1];
sx q[1];
rz(-1.1670246) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052159667) q[0];
sx q[0];
rz(-2.8773327) q[0];
sx q[0];
rz(1.7961851) q[0];
rz(-2.8950047) q[2];
sx q[2];
rz(-1.0435259) q[2];
sx q[2];
rz(-1.2621244) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0868623) q[1];
sx q[1];
rz(-0.80838258) q[1];
sx q[1];
rz(-0.085309172) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59252177) q[3];
sx q[3];
rz(-0.7624818) q[3];
sx q[3];
rz(-3.1069063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.033826753) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(2.9272184) q[2];
rz(-3.068148) q[3];
sx q[3];
rz(-0.44973222) q[3];
sx q[3];
rz(0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6630994) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(-1.9146772) q[0];
rz(-0.40027174) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(2.1267557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0815711) q[0];
sx q[0];
rz(-0.9284174) q[0];
sx q[0];
rz(-1.6409671) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8995908) q[2];
sx q[2];
rz(-1.1761464) q[2];
sx q[2];
rz(1.637527) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.240757) q[1];
sx q[1];
rz(-1.6512617) q[1];
sx q[1];
rz(2.1686173) q[1];
x q[2];
rz(0.88481836) q[3];
sx q[3];
rz(-1.775145) q[3];
sx q[3];
rz(0.94061461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0959452) q[2];
sx q[2];
rz(-1.0568876) q[2];
sx q[2];
rz(-0.8992368) q[2];
rz(-0.69747654) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(2.5337059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65524453) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(-0.43193257) q[0];
rz(0.63255429) q[1];
sx q[1];
rz(-0.41699854) q[1];
sx q[1];
rz(0.63582173) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83946562) q[0];
sx q[0];
rz(-1.7502022) q[0];
sx q[0];
rz(1.0666215) q[0];
rz(-1.7556778) q[2];
sx q[2];
rz(-0.73542483) q[2];
sx q[2];
rz(0.76279574) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7829261) q[1];
sx q[1];
rz(-1.7554605) q[1];
sx q[1];
rz(-0.77628805) q[1];
rz(-pi) q[2];
rz(2.5084247) q[3];
sx q[3];
rz(-0.35053262) q[3];
sx q[3];
rz(-1.4299973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2137961) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(-2.4528743) q[2];
rz(2.8074746) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(2.9978602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9976945) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(1.0744263) q[0];
rz(2.396446) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(-0.27854663) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51019788) q[0];
sx q[0];
rz(-0.83980951) q[0];
sx q[0];
rz(2.1561119) q[0];
rz(-2.5299046) q[2];
sx q[2];
rz(-1.4809161) q[2];
sx q[2];
rz(0.22345605) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3814195) q[1];
sx q[1];
rz(-0.96254327) q[1];
sx q[1];
rz(0.45254032) q[1];
rz(-pi) q[2];
rz(-0.99028011) q[3];
sx q[3];
rz(-2.2970082) q[3];
sx q[3];
rz(0.57675225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6986065) q[2];
sx q[2];
rz(-1.7437982) q[2];
sx q[2];
rz(2.8473575) q[2];
rz(3.0596628) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(3.1053655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1086403) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(3.1325353) q[0];
rz(2.5065705) q[1];
sx q[1];
rz(-0.68990866) q[1];
sx q[1];
rz(3.0335398) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6271583) q[0];
sx q[0];
rz(-1.8875727) q[0];
sx q[0];
rz(0.47382521) q[0];
x q[1];
rz(0.33424218) q[2];
sx q[2];
rz(-0.40194449) q[2];
sx q[2];
rz(-0.7451171) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15247004) q[1];
sx q[1];
rz(-1.9580541) q[1];
sx q[1];
rz(-0.5247922) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8764358) q[3];
sx q[3];
rz(-1.0401871) q[3];
sx q[3];
rz(-1.6924072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1487427) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(-1.8704869) q[2];
rz(-3.0631915) q[3];
sx q[3];
rz(-1.4784808) q[3];
sx q[3];
rz(-1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8909797) q[0];
sx q[0];
rz(-1.5654634) q[0];
sx q[0];
rz(2.4839731) q[0];
rz(-1.744386) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(0.89362842) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2715627) q[0];
sx q[0];
rz(-1.8989519) q[0];
sx q[0];
rz(-0.66332711) q[0];
rz(2.0486084) q[2];
sx q[2];
rz(-1.8336979) q[2];
sx q[2];
rz(-2.773657) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.5205198) q[1];
sx q[1];
rz(-0.90172651) q[1];
sx q[1];
rz(-1.8280162) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10223933) q[3];
sx q[3];
rz(-0.7623626) q[3];
sx q[3];
rz(1.2650714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7509193) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(-2.612109) q[2];
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
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3787518) q[0];
sx q[0];
rz(-1.6289926) q[0];
sx q[0];
rz(-2.0513127) q[0];
rz(-3.0293363) q[1];
sx q[1];
rz(-2.0376164) q[1];
sx q[1];
rz(-1.9876678) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9268122) q[0];
sx q[0];
rz(-1.6048389) q[0];
sx q[0];
rz(2.4916617) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9032193) q[2];
sx q[2];
rz(-2.5775238) q[2];
sx q[2];
rz(-2.7114045) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1057508) q[1];
sx q[1];
rz(-2.3301947) q[1];
sx q[1];
rz(1.6559385) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38903799) q[3];
sx q[3];
rz(-0.93718796) q[3];
sx q[3];
rz(0.47215677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5899137) q[2];
sx q[2];
rz(-0.38107291) q[2];
sx q[2];
rz(-2.7611458) q[2];
rz(2.0137265) q[3];
sx q[3];
rz(-1.7815855) q[3];
sx q[3];
rz(-1.5650704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34148759) q[0];
sx q[0];
rz(-1.2416168) q[0];
sx q[0];
rz(-2.5019116) q[0];
rz(-1.9027963) q[1];
sx q[1];
rz(-1.4150554) q[1];
sx q[1];
rz(1.9715086) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3568048) q[0];
sx q[0];
rz(-2.6492282) q[0];
sx q[0];
rz(0.82226336) q[0];
rz(-2.8662888) q[2];
sx q[2];
rz(-0.79553662) q[2];
sx q[2];
rz(-0.77992935) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.942109) q[1];
sx q[1];
rz(-2.0523199) q[1];
sx q[1];
rz(-1.3660618) q[1];
rz(-0.13799237) q[3];
sx q[3];
rz(-2.237461) q[3];
sx q[3];
rz(0.43825144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3141979) q[2];
sx q[2];
rz(-2.4908227) q[2];
sx q[2];
rz(0.38273746) q[2];
rz(0.9283723) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(-0.66463566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9160354) q[0];
sx q[0];
rz(-1.6397497) q[0];
sx q[0];
rz(2.8826707) q[0];
rz(-2.4312773) q[1];
sx q[1];
rz(-2.0878891) q[1];
sx q[1];
rz(-2.6616667) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7326638) q[0];
sx q[0];
rz(-0.73584475) q[0];
sx q[0];
rz(1.0622513) q[0];
rz(-pi) q[1];
rz(-2.2357113) q[2];
sx q[2];
rz(-1.3833628) q[2];
sx q[2];
rz(2.4311709) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.5727947) q[1];
sx q[1];
rz(-0.63981445) q[1];
sx q[1];
rz(-2.5261643) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89703538) q[3];
sx q[3];
rz(-1.1437136) q[3];
sx q[3];
rz(0.87891146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2991128) q[2];
sx q[2];
rz(-2.2128236) q[2];
sx q[2];
rz(1.8245565) q[2];
rz(1.8995829) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(-2.4035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15923545) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(0.97933979) q[1];
sx q[1];
rz(-1.0995438) q[1];
sx q[1];
rz(-0.88811036) q[1];
rz(-0.48158823) q[2];
sx q[2];
rz(-1.7272186) q[2];
sx q[2];
rz(-2.061736) q[2];
rz(-0.76673037) q[3];
sx q[3];
rz(-0.37692108) q[3];
sx q[3];
rz(2.2212096) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
