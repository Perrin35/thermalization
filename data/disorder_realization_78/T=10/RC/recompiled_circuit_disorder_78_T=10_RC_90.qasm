OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5821563) q[0];
sx q[0];
rz(-0.59098935) q[0];
sx q[0];
rz(0.58340573) q[0];
rz(2.9572339) q[1];
sx q[1];
rz(-0.98362041) q[1];
sx q[1];
rz(-0.89259994) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8164506) q[0];
sx q[0];
rz(-0.6439552) q[0];
sx q[0];
rz(-2.08026) q[0];
rz(-pi) q[1];
rz(0.61172723) q[2];
sx q[2];
rz(-0.76843843) q[2];
sx q[2];
rz(-0.4380463) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.8436369) q[1];
sx q[1];
rz(-1.2601818) q[1];
sx q[1];
rz(-0.85308869) q[1];
rz(-pi) q[2];
rz(-1.6139908) q[3];
sx q[3];
rz(-0.82026635) q[3];
sx q[3];
rz(-0.84212069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9154174) q[2];
sx q[2];
rz(-1.6025275) q[2];
sx q[2];
rz(1.9809451) q[2];
rz(-0.21696572) q[3];
sx q[3];
rz(-0.522885) q[3];
sx q[3];
rz(-1.0552361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.083667) q[0];
sx q[0];
rz(-2.1766429) q[0];
sx q[0];
rz(-0.57587409) q[0];
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
rz(0.28536797) q[0];
sx q[0];
rz(-1.8282187) q[0];
sx q[0];
rz(3.0811937) q[0];
rz(-1.1738271) q[2];
sx q[2];
rz(-2.5644828) q[2];
sx q[2];
rz(1.4156262) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0547304) q[1];
sx q[1];
rz(-2.3332101) q[1];
sx q[1];
rz(-3.0562835) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4715273) q[3];
sx q[3];
rz(-1.1747922) q[3];
sx q[3];
rz(1.083064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.033826753) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(-2.9272184) q[2];
rz(3.068148) q[3];
sx q[3];
rz(-2.6918604) q[3];
sx q[3];
rz(-2.8607821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630994) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(1.9146772) q[0];
rz(2.7413209) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(2.1267557) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44714156) q[0];
sx q[0];
rz(-1.6269636) q[0];
sx q[0];
rz(0.64356128) q[0];
rz(-pi) q[1];
rz(-2.8995908) q[2];
sx q[2];
rz(-1.9654462) q[2];
sx q[2];
rz(1.637527) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.44753578) q[1];
sx q[1];
rz(-2.5390365) q[1];
sx q[1];
rz(1.4284929) q[1];
x q[2];
rz(0.88481836) q[3];
sx q[3];
rz(-1.775145) q[3];
sx q[3];
rz(0.94061461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0456475) q[2];
sx q[2];
rz(-2.084705) q[2];
sx q[2];
rz(0.8992368) q[2];
rz(-0.69747654) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(-0.60788679) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65524453) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(-2.7096601) q[0];
rz(0.63255429) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(-0.63582173) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.302127) q[0];
sx q[0];
rz(-1.3913904) q[0];
sx q[0];
rz(2.0749712) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16480883) q[2];
sx q[2];
rz(-2.2909082) q[2];
sx q[2];
rz(1.009843) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7446049) q[1];
sx q[1];
rz(-0.79345353) q[1];
sx q[1];
rz(0.26054392) q[1];
rz(-1.7838578) q[3];
sx q[3];
rz(-1.8512929) q[3];
sx q[3];
rz(2.3749542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2137961) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(0.68871838) q[2];
rz(-0.33411807) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(2.9978602) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14389811) q[0];
sx q[0];
rz(-2.4425638) q[0];
sx q[0];
rz(-2.0671663) q[0];
rz(2.396446) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(-0.27854663) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8700096) q[0];
sx q[0];
rz(-0.9013114) q[0];
sx q[0];
rz(-0.55218009) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4611545) q[2];
sx q[2];
rz(-2.1796558) q[2];
sx q[2];
rz(-1.4102175) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.67644557) q[1];
sx q[1];
rz(-2.4009581) q[1];
sx q[1];
rz(-2.1315104) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5883702) q[3];
sx q[3];
rz(-0.8953989) q[3];
sx q[3];
rz(-1.355987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6986065) q[2];
sx q[2];
rz(-1.7437982) q[2];
sx q[2];
rz(-0.29423514) q[2];
rz(0.081929835) q[3];
sx q[3];
rz(-2.6223845) q[3];
sx q[3];
rz(3.1053655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1086403) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(-3.1325353) q[0];
rz(0.63502216) q[1];
sx q[1];
rz(-0.68990866) q[1];
sx q[1];
rz(-3.0335398) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2436284) q[0];
sx q[0];
rz(-2.0192696) q[0];
sx q[0];
rz(1.2178221) q[0];
x q[1];
rz(-1.7093541) q[2];
sx q[2];
rz(-1.9493305) q[2];
sx q[2];
rz(1.1059424) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5080155) q[1];
sx q[1];
rz(-1.0883691) q[1];
sx q[1];
rz(1.1303933) q[1];
rz(-pi) q[2];
rz(2.8764358) q[3];
sx q[3];
rz(-1.0401871) q[3];
sx q[3];
rz(1.6924072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99284995) q[2];
sx q[2];
rz(-1.7603346) q[2];
sx q[2];
rz(1.2711058) q[2];
rz(3.0631915) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(-1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25061297) q[0];
sx q[0];
rz(-1.5654634) q[0];
sx q[0];
rz(-2.4839731) q[0];
rz(1.744386) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(-0.89362842) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2715627) q[0];
sx q[0];
rz(-1.2426408) q[0];
sx q[0];
rz(0.66332711) q[0];
x q[1];
rz(-0.29427476) q[2];
sx q[2];
rz(-2.0308959) q[2];
sx q[2];
rz(-1.0690881) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0221755) q[1];
sx q[1];
rz(-0.70964538) q[1];
sx q[1];
rz(0.31125734) q[1];
rz(-pi) q[2];
rz(1.4736389) q[3];
sx q[3];
rz(-2.328184) q[3];
sx q[3];
rz(-1.1241476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7509193) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(0.52948362) q[2];
rz(-2.6654065) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(2.7990394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7628409) q[0];
sx q[0];
rz(-1.6289926) q[0];
sx q[0];
rz(-1.09028) q[0];
rz(3.0293363) q[1];
sx q[1];
rz(-2.0376164) q[1];
sx q[1];
rz(1.9876678) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.740828) q[0];
sx q[0];
rz(-2.4908998) q[0];
sx q[0];
rz(3.0853737) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4225142) q[2];
sx q[2];
rz(-2.1170756) q[2];
sx q[2];
rz(-0.15020457) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47626074) q[1];
sx q[1];
rz(-1.5090824) q[1];
sx q[1];
rz(-2.3803821) q[1];
x q[2];
rz(-1.0941986) q[3];
sx q[3];
rz(-0.72924858) q[3];
sx q[3];
rz(0.13343982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.551679) q[2];
sx q[2];
rz(-0.38107291) q[2];
sx q[2];
rz(-2.7611458) q[2];
rz(1.1278661) q[3];
sx q[3];
rz(-1.3600072) q[3];
sx q[3];
rz(-1.5650704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8001051) q[0];
sx q[0];
rz(-1.2416168) q[0];
sx q[0];
rz(-0.63968101) q[0];
rz(1.2387964) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(-1.9715086) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6696475) q[0];
sx q[0];
rz(-1.8983316) q[0];
sx q[0];
rz(1.1963084) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3001928) q[2];
sx q[2];
rz(-0.81297183) q[2];
sx q[2];
rz(-1.1635309) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.942109) q[1];
sx q[1];
rz(-1.0892727) q[1];
sx q[1];
rz(-1.7755309) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2421078) q[3];
sx q[3];
rz(-1.4624819) q[3];
sx q[3];
rz(1.9233821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8273948) q[2];
sx q[2];
rz(-0.65076995) q[2];
sx q[2];
rz(-2.7588552) q[2];
rz(0.9283723) q[3];
sx q[3];
rz(-1.174077) q[3];
sx q[3];
rz(-2.476957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-2.9160354) q[0];
sx q[0];
rz(-1.6397497) q[0];
sx q[0];
rz(-2.8826707) q[0];
rz(-0.71031538) q[1];
sx q[1];
rz(-1.0537035) q[1];
sx q[1];
rz(-2.6616667) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23586789) q[0];
sx q[0];
rz(-0.94434443) q[0];
sx q[0];
rz(-0.41525526) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8690228) q[2];
sx q[2];
rz(-2.4546461) q[2];
sx q[2];
rz(0.62703122) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2952134) q[1];
sx q[1];
rz(-1.0615674) q[1];
sx q[1];
rz(1.9766115) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2445573) q[3];
sx q[3];
rz(-1.1437136) q[3];
sx q[3];
rz(-2.2626812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2991128) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(1.3170362) q[2];
rz(-1.2420098) q[3];
sx q[3];
rz(-2.1879991) q[3];
sx q[3];
rz(2.4035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15923545) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(-2.1622529) q[1];
sx q[1];
rz(-1.0995438) q[1];
sx q[1];
rz(-0.88811036) q[1];
rz(2.813415) q[2];
sx q[2];
rz(-2.6371418) q[2];
sx q[2];
rz(-0.20124659) q[2];
rz(-0.27771523) q[3];
sx q[3];
rz(-1.3125827) q[3];
sx q[3];
rz(-1.7607341) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
