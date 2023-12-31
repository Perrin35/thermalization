OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8945583) q[0];
sx q[0];
rz(-1.2556827) q[0];
sx q[0];
rz(-0.32796252) q[0];
rz(-0.2344996) q[1];
sx q[1];
rz(3.3426715) q[1];
sx q[1];
rz(9.3333416) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2823921) q[0];
sx q[0];
rz(-0.99635591) q[0];
sx q[0];
rz(-1.0919071) q[0];
rz(0.25912164) q[2];
sx q[2];
rz(-1.4403575) q[2];
sx q[2];
rz(-2.5082617) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5603197) q[1];
sx q[1];
rz(-0.25174192) q[1];
sx q[1];
rz(-1.9074349) q[1];
rz(-pi) q[2];
rz(2.1655733) q[3];
sx q[3];
rz(-1.3888437) q[3];
sx q[3];
rz(2.6973157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4771007) q[2];
sx q[2];
rz(-0.97350073) q[2];
sx q[2];
rz(1.1260024) q[2];
rz(-2.8664355) q[3];
sx q[3];
rz(-2.5313009) q[3];
sx q[3];
rz(2.2385105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098410957) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(-2.4480208) q[0];
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
rz(-1.2454525) q[0];
sx q[0];
rz(-1.3724694) q[0];
sx q[0];
rz(1.1597) q[0];
x q[1];
rz(-2.6066577) q[2];
sx q[2];
rz(-1.2282279) q[2];
sx q[2];
rz(1.554622) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.43860501) q[1];
sx q[1];
rz(-0.72241306) q[1];
sx q[1];
rz(-1.5370675) q[1];
x q[2];
rz(0.5204366) q[3];
sx q[3];
rz(-2.3428829) q[3];
sx q[3];
rz(-0.48898104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.50743121) q[2];
sx q[2];
rz(-2.4298411) q[2];
sx q[2];
rz(1.2197536) q[2];
rz(0.1427342) q[3];
sx q[3];
rz(-1.0120564) q[3];
sx q[3];
rz(-2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74293566) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(-2.3213342) q[0];
rz(0.29207686) q[1];
sx q[1];
rz(-1.074011) q[1];
sx q[1];
rz(1.2480199) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1093729) q[0];
sx q[0];
rz(-1.6633699) q[0];
sx q[0];
rz(2.5412482) q[0];
x q[1];
rz(-1.6349995) q[2];
sx q[2];
rz(-1.3639796) q[2];
sx q[2];
rz(1.1484255) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4856845) q[1];
sx q[1];
rz(-2.471172) q[1];
sx q[1];
rz(-0.1078492) q[1];
rz(-pi) q[2];
rz(-1.9565342) q[3];
sx q[3];
rz(-0.68813656) q[3];
sx q[3];
rz(-3.0582173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8212006) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(1.2134264) q[2];
rz(0.16472566) q[3];
sx q[3];
rz(-0.89403331) q[3];
sx q[3];
rz(-0.081610672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7320025) q[0];
sx q[0];
rz(-1.2816757) q[0];
sx q[0];
rz(-2.9555292) q[0];
rz(-2.9371254) q[1];
sx q[1];
rz(-0.47195131) q[1];
sx q[1];
rz(-1.2971372) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7521833) q[0];
sx q[0];
rz(-1.4404738) q[0];
sx q[0];
rz(3.0979063) q[0];
rz(-pi) q[1];
rz(-1.859971) q[2];
sx q[2];
rz(-1.7316069) q[2];
sx q[2];
rz(-0.5493872) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2343826) q[1];
sx q[1];
rz(-1.0115336) q[1];
sx q[1];
rz(2.2940966) q[1];
x q[2];
rz(0.013838776) q[3];
sx q[3];
rz(-1.3912364) q[3];
sx q[3];
rz(1.1650712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2356448) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(-2.2223991) q[2];
rz(0.32133189) q[3];
sx q[3];
rz(-2.0733757) q[3];
sx q[3];
rz(1.2478158) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17386757) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(-2.2633973) q[0];
rz(-1.8163266) q[1];
sx q[1];
rz(-1.7293431) q[1];
sx q[1];
rz(1.7153046) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.250524) q[0];
sx q[0];
rz(-0.91618012) q[0];
sx q[0];
rz(3.0380681) q[0];
rz(-1.8833141) q[2];
sx q[2];
rz(-1.3442355) q[2];
sx q[2];
rz(2.8001919) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9702455) q[1];
sx q[1];
rz(-0.73927021) q[1];
sx q[1];
rz(-1.4096178) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32894965) q[3];
sx q[3];
rz(-2.3429686) q[3];
sx q[3];
rz(0.22508276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8695801) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(0.041794725) q[2];
rz(-0.061491866) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(-2.7048236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3549266) q[0];
sx q[0];
rz(-3.0556398) q[0];
sx q[0];
rz(2.1110995) q[0];
rz(0.73973918) q[1];
sx q[1];
rz(-1.5286427) q[1];
sx q[1];
rz(2.5700263) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8067236) q[0];
sx q[0];
rz(-2.3248701) q[0];
sx q[0];
rz(2.4094894) q[0];
rz(-pi) q[1];
rz(2.2366441) q[2];
sx q[2];
rz(-1.6709423) q[2];
sx q[2];
rz(-3.1040994) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7476269) q[1];
sx q[1];
rz(-2.6193301) q[1];
sx q[1];
rz(2.0678492) q[1];
rz(-pi) q[2];
rz(-2.3011175) q[3];
sx q[3];
rz(-2.2395036) q[3];
sx q[3];
rz(0.83276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.968154) q[2];
sx q[2];
rz(-2.3996694) q[2];
sx q[2];
rz(0.56419939) q[2];
rz(0.12600222) q[3];
sx q[3];
rz(-1.4583476) q[3];
sx q[3];
rz(-2.8469767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903704) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(2.4293161) q[0];
rz(2.6157216) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(0.66551048) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3367046) q[0];
sx q[0];
rz(-0.24222736) q[0];
sx q[0];
rz(-1.5750984) q[0];
rz(1.7475142) q[2];
sx q[2];
rz(-2.4364947) q[2];
sx q[2];
rz(1.2202028) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9741386) q[1];
sx q[1];
rz(-1.1145089) q[1];
sx q[1];
rz(-0.2445226) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3286367) q[3];
sx q[3];
rz(-1.4302963) q[3];
sx q[3];
rz(-1.5797918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9843288) q[2];
sx q[2];
rz(-2.4815346) q[2];
sx q[2];
rz(-1.7187913) q[2];
rz(3.026399) q[3];
sx q[3];
rz(-2.6264103) q[3];
sx q[3];
rz(0.19259024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049906235) q[0];
sx q[0];
rz(-1.1356857) q[0];
sx q[0];
rz(-0.19083047) q[0];
rz(2.514839) q[1];
sx q[1];
rz(-2.125506) q[1];
sx q[1];
rz(-0.33871067) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3074293) q[0];
sx q[0];
rz(-1.2273664) q[0];
sx q[0];
rz(-2.9817581) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4523925) q[2];
sx q[2];
rz(-2.6466742) q[2];
sx q[2];
rz(-2.6099043) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.415886) q[1];
sx q[1];
rz(-0.48735122) q[1];
sx q[1];
rz(-1.8094256) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8645105) q[3];
sx q[3];
rz(-1.2556561) q[3];
sx q[3];
rz(-1.806123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4954341) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(0.70043606) q[2];
rz(-0.8979848) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(-2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.609628) q[0];
sx q[0];
rz(-2.6599929) q[0];
sx q[0];
rz(-0.65834808) q[0];
rz(-0.61093962) q[1];
sx q[1];
rz(-1.2555723) q[1];
sx q[1];
rz(-3.0019965) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1698648) q[0];
sx q[0];
rz(-1.5507878) q[0];
sx q[0];
rz(-1.5097029) q[0];
rz(3.082824) q[2];
sx q[2];
rz(-2.3261056) q[2];
sx q[2];
rz(-1.2440484) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82397205) q[1];
sx q[1];
rz(-0.52597731) q[1];
sx q[1];
rz(-3.0082263) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0890907) q[3];
sx q[3];
rz(-2.2146261) q[3];
sx q[3];
rz(2.1569463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6167986) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(-2.810478) q[2];
rz(2.3838499) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(-0.087879114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.2963294) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(-2.9560126) q[0];
rz(-2.045385) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(-1.4846444) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.164924) q[0];
sx q[0];
rz(-2.6161368) q[0];
sx q[0];
rz(-2.1303961) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22557232) q[2];
sx q[2];
rz(-1.8403056) q[2];
sx q[2];
rz(3.0849506) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.46586043) q[1];
sx q[1];
rz(-1.670174) q[1];
sx q[1];
rz(1.194721) q[1];
x q[2];
rz(-3.1366523) q[3];
sx q[3];
rz(-0.86588174) q[3];
sx q[3];
rz(-2.4019394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2075656) q[2];
sx q[2];
rz(-0.88946122) q[2];
sx q[2];
rz(0.55220848) q[2];
rz(-2.3637555) q[3];
sx q[3];
rz(-2.2838897) q[3];
sx q[3];
rz(0.51789969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8778397) q[0];
sx q[0];
rz(-1.629566) q[0];
sx q[0];
rz(-1.4019479) q[0];
rz(-2.9539625) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(-0.4734584) q[2];
sx q[2];
rz(-2.8291694) q[2];
sx q[2];
rz(1.3419801) q[2];
rz(-0.68708146) q[3];
sx q[3];
rz(-1.0071181) q[3];
sx q[3];
rz(-1.3345171) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
