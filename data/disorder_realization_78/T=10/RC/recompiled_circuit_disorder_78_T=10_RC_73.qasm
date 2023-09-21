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
rz(3.732582) q[0];
sx q[0];
rz(8.8413722) q[0];
rz(-0.18435873) q[1];
sx q[1];
rz(4.1252131) q[1];
sx q[1];
rz(10.317378) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8569782) q[0];
sx q[0];
rz(-1.0191139) q[0];
sx q[0];
rz(-2.7906228) q[0];
rz(-pi) q[1];
rz(2.0775665) q[2];
sx q[2];
rz(-0.96553409) q[2];
sx q[2];
rz(1.2111226) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0639227) q[1];
sx q[1];
rz(-0.7709255) q[1];
sx q[1];
rz(2.0248807) q[1];
rz(-3.0953232) q[3];
sx q[3];
rz(-0.75152961) q[3];
sx q[3];
rz(2.3627594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9154174) q[2];
sx q[2];
rz(-1.6025275) q[2];
sx q[2];
rz(1.1606476) q[2];
rz(-0.21696572) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(-2.0863566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.0579257) q[0];
sx q[0];
rz(-0.96494976) q[0];
sx q[0];
rz(0.57587409) q[0];
rz(-1.8946164) q[1];
sx q[1];
rz(-1.2966825) q[1];
sx q[1];
rz(-1.974568) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8407699) q[0];
sx q[0];
rz(-1.5123899) q[0];
sx q[0];
rz(-1.8286684) q[0];
x q[1];
rz(-1.1738271) q[2];
sx q[2];
rz(-2.5644828) q[2];
sx q[2];
rz(-1.7259665) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.45707073) q[1];
sx q[1];
rz(-1.6324537) q[1];
sx q[1];
rz(0.80656273) q[1];
x q[2];
rz(-0.59252177) q[3];
sx q[3];
rz(-2.3791109) q[3];
sx q[3];
rz(3.1069063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1077659) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(-2.9272184) q[2];
rz(3.068148) q[3];
sx q[3];
rz(-2.6918604) q[3];
sx q[3];
rz(0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6630994) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(1.2269155) q[0];
rz(-0.40027174) q[1];
sx q[1];
rz(-1.8881256) q[1];
sx q[1];
rz(-2.1267557) q[1];
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
rz(-1.1655408) q[2];
sx q[2];
rz(-1.793867) q[2];
sx q[2];
rz(-0.16135339) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.240757) q[1];
sx q[1];
rz(-1.490331) q[1];
sx q[1];
rz(-0.97297538) q[1];
x q[2];
rz(2.879911) q[3];
sx q[3];
rz(-0.90173429) q[3];
sx q[3];
rz(-0.79479587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0456475) q[2];
sx q[2];
rz(-2.084705) q[2];
sx q[2];
rz(2.2423559) q[2];
rz(-0.69747654) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(2.5337059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65524453) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(-0.43193257) q[0];
rz(-2.5090384) q[1];
sx q[1];
rz(-0.41699854) q[1];
sx q[1];
rz(-2.5057709) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5084002) q[0];
sx q[0];
rz(-2.0661372) q[0];
sx q[0];
rz(0.20423996) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7556778) q[2];
sx q[2];
rz(-2.4061678) q[2];
sx q[2];
rz(2.3787969) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1078474) q[1];
sx q[1];
rz(-0.81106942) q[1];
sx q[1];
rz(1.3147522) q[1];
rz(-1.7838578) q[3];
sx q[3];
rz(-1.8512929) q[3];
sx q[3];
rz(-0.76663843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2137961) q[2];
sx q[2];
rz(-0.14363229) q[2];
sx q[2];
rz(-2.4528743) q[2];
rz(-0.33411807) q[3];
sx q[3];
rz(-1.9709316) q[3];
sx q[3];
rz(0.14373246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9976945) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(-2.0671663) q[0];
rz(-2.396446) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(0.27854663) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2715831) q[0];
sx q[0];
rz(-2.2402813) q[0];
sx q[0];
rz(2.5894126) q[0];
rz(-pi) q[1];
rz(2.9859221) q[2];
sx q[2];
rz(-2.5241733) q[2];
sx q[2];
rz(1.9215259) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67644557) q[1];
sx q[1];
rz(-2.4009581) q[1];
sx q[1];
rz(1.0100823) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3260818) q[3];
sx q[3];
rz(-1.1482571) q[3];
sx q[3];
rz(0.58327196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6022588) q[0];
sx q[0];
rz(-0.56319153) q[0];
sx q[0];
rz(-2.5186033) q[0];
rz(-pi) q[1];
rz(-2.8073505) q[2];
sx q[2];
rz(-2.7396482) q[2];
sx q[2];
rz(-2.3964756) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.84032816) q[1];
sx q[1];
rz(-2.5003308) q[1];
sx q[1];
rz(-2.458359) q[1];
rz(-pi) q[2];
rz(-0.26515682) q[3];
sx q[3];
rz(-1.0401871) q[3];
sx q[3];
rz(-1.4491855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1487427) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(1.8704869) q[2];
rz(-0.078401119) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(-1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(2.8909797) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(-0.65761956) q[0];
rz(1.744386) q[1];
sx q[1];
rz(-1.9669292) q[1];
sx q[1];
rz(2.2479642) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2715627) q[0];
sx q[0];
rz(-1.2426408) q[0];
sx q[0];
rz(2.4782655) q[0];
rz(-2.8473179) q[2];
sx q[2];
rz(-2.0308959) q[2];
sx q[2];
rz(1.0690881) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.5205198) q[1];
sx q[1];
rz(-0.90172651) q[1];
sx q[1];
rz(-1.8280162) q[1];
x q[2];
rz(0.10223933) q[3];
sx q[3];
rz(-0.7623626) q[3];
sx q[3];
rz(1.2650714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.39067337) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(2.612109) q[2];
rz(2.6654065) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(0.34255323) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7628409) q[0];
sx q[0];
rz(-1.5126001) q[0];
sx q[0];
rz(2.0513127) q[0];
rz(-0.11225637) q[1];
sx q[1];
rz(-2.0376164) q[1];
sx q[1];
rz(1.9876678) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.740828) q[0];
sx q[0];
rz(-2.4908998) q[0];
sx q[0];
rz(-0.056218938) q[0];
rz(-pi) q[1];
rz(0.23837337) q[2];
sx q[2];
rz(-2.5775238) q[2];
sx q[2];
rz(-0.43018815) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1591783) q[1];
sx q[1];
rz(-2.3783861) q[1];
sx q[1];
rz(0.089341954) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0941986) q[3];
sx q[3];
rz(-2.4123441) q[3];
sx q[3];
rz(-0.13343982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.551679) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(-2.7611458) q[2];
rz(-1.1278661) q[3];
sx q[3];
rz(-1.7815855) q[3];
sx q[3];
rz(-1.5650704) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8001051) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(2.5019116) q[0];
rz(1.9027963) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(1.9715086) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026924883) q[0];
sx q[0];
rz(-1.2171193) q[0];
sx q[0];
rz(-2.7915519) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27530382) q[2];
sx q[2];
rz(-0.79553662) q[2];
sx q[2];
rz(0.77992935) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3635892) q[1];
sx q[1];
rz(-0.52007857) q[1];
sx q[1];
rz(0.37104721) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3977259) q[3];
sx q[3];
rz(-2.462938) q[3];
sx q[3];
rz(-0.21733397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8273948) q[2];
sx q[2];
rz(-2.4908227) q[2];
sx q[2];
rz(2.7588552) q[2];
rz(2.2132204) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(-2.476957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2255573) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(-0.25892192) q[0];
rz(-0.71031538) q[1];
sx q[1];
rz(-2.0878891) q[1];
sx q[1];
rz(-0.47992596) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23586789) q[0];
sx q[0];
rz(-2.1972482) q[0];
sx q[0];
rz(-2.7263374) q[0];
rz(-pi) q[1];
rz(0.9058814) q[2];
sx q[2];
rz(-1.3833628) q[2];
sx q[2];
rz(2.4311709) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.568798) q[1];
sx q[1];
rz(-0.63981445) q[1];
sx q[1];
rz(2.5261643) q[1];
rz(-pi) q[2];
rz(0.94060937) q[3];
sx q[3];
rz(-0.77946957) q[3];
sx q[3];
rz(-1.1704695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2991128) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(1.8245565) q[2];
rz(1.8995829) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(0.73808134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.15923545) q[0];
sx q[0];
rz(-0.032082162) q[0];
sx q[0];
rz(-1.4627009) q[0];
rz(-2.1622529) q[1];
sx q[1];
rz(-1.0995438) q[1];
sx q[1];
rz(-0.88811036) q[1];
rz(2.6600044) q[2];
sx q[2];
rz(-1.7272186) q[2];
sx q[2];
rz(-2.061736) q[2];
rz(0.27771523) q[3];
sx q[3];
rz(-1.82901) q[3];
sx q[3];
rz(1.3808586) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];