OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6501939) q[0];
sx q[0];
rz(-2.8770652) q[0];
sx q[0];
rz(-2.7471623) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(2.8013464) q[1];
sx q[1];
rz(10.624788) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7317176) q[0];
sx q[0];
rz(-0.048225064) q[0];
sx q[0];
rz(-1.4531141) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8067402) q[2];
sx q[2];
rz(-1.9960253) q[2];
sx q[2];
rz(-1.2644757) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4760121) q[1];
sx q[1];
rz(-1.952938) q[1];
sx q[1];
rz(2.0631454) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6148189) q[3];
sx q[3];
rz(-1.9485954) q[3];
sx q[3];
rz(1.6954741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8346943) q[2];
sx q[2];
rz(-1.6480185) q[2];
sx q[2];
rz(-0.28960323) q[2];
rz(-2.2662207) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(-3.0818821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5263379) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(-2.7976024) q[0];
rz(-0.084331766) q[1];
sx q[1];
rz(-0.66939676) q[1];
sx q[1];
rz(1.7864236) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25027572) q[0];
sx q[0];
rz(-2.180897) q[0];
sx q[0];
rz(-0.079205714) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9789574) q[2];
sx q[2];
rz(-1.4504823) q[2];
sx q[2];
rz(-1.6441117) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4210216) q[1];
sx q[1];
rz(-0.45383006) q[1];
sx q[1];
rz(-1.717091) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1739028) q[3];
sx q[3];
rz(-0.19860425) q[3];
sx q[3];
rz(-1.393115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29558674) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(-0.53768349) q[2];
rz(-0.50283557) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(-2.1311549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4780592) q[0];
sx q[0];
rz(-0.53776598) q[0];
sx q[0];
rz(0.64087254) q[0];
rz(2.3928941) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(-2.0764988) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.051620313) q[0];
sx q[0];
rz(-1.211477) q[0];
sx q[0];
rz(-0.006747147) q[0];
x q[1];
rz(-0.87419072) q[2];
sx q[2];
rz(-2.0421931) q[2];
sx q[2];
rz(3.0490321) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0703735) q[1];
sx q[1];
rz(-2.1871236) q[1];
sx q[1];
rz(0.044205772) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84824003) q[3];
sx q[3];
rz(-1.7602929) q[3];
sx q[3];
rz(3.0998067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37725267) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(2.4242145) q[2];
rz(-0.68850368) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(0.29754105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82729572) q[0];
sx q[0];
rz(-1.941444) q[0];
sx q[0];
rz(-2.8919019) q[0];
rz(-1.0149792) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(0.011118523) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3664368) q[0];
sx q[0];
rz(-1.1661134) q[0];
sx q[0];
rz(1.0767656) q[0];
rz(-pi) q[1];
rz(2.501802) q[2];
sx q[2];
rz(-2.1961229) q[2];
sx q[2];
rz(1.7196136) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3938155) q[1];
sx q[1];
rz(-1.897246) q[1];
sx q[1];
rz(1.0541037) q[1];
rz(1.7632145) q[3];
sx q[3];
rz(-1.8835861) q[3];
sx q[3];
rz(-0.49945143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.49542385) q[2];
sx q[2];
rz(-1.255722) q[2];
sx q[2];
rz(-2.8584976) q[2];
rz(2.4781573) q[3];
sx q[3];
rz(-0.59316558) q[3];
sx q[3];
rz(2.2535113) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11113142) q[0];
sx q[0];
rz(-2.9009394) q[0];
sx q[0];
rz(-2.8097613) q[0];
rz(0.49452531) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(-1.3269075) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3834284) q[0];
sx q[0];
rz(-0.82793068) q[0];
sx q[0];
rz(-1.7799671) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1660216) q[2];
sx q[2];
rz(-2.274548) q[2];
sx q[2];
rz(0.032035839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8718308) q[1];
sx q[1];
rz(-2.7360536) q[1];
sx q[1];
rz(-2.521442) q[1];
rz(1.1633478) q[3];
sx q[3];
rz(-1.8383998) q[3];
sx q[3];
rz(2.101055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.999324) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(-1.2456606) q[2];
rz(1.2549531) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(-2.3033223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44928837) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(0.6814878) q[0];
rz(2.9340414) q[1];
sx q[1];
rz(-0.46336585) q[1];
sx q[1];
rz(2.025827) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078243144) q[0];
sx q[0];
rz(-0.82394281) q[0];
sx q[0];
rz(-2.1372165) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0461147) q[2];
sx q[2];
rz(-0.74138734) q[2];
sx q[2];
rz(-1.807883) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7129732) q[1];
sx q[1];
rz(-1.5884807) q[1];
sx q[1];
rz(2.5531205) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3917771) q[3];
sx q[3];
rz(-2.2324623) q[3];
sx q[3];
rz(0.20565198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9289124) q[2];
sx q[2];
rz(-1.8447515) q[2];
sx q[2];
rz(2.7872655) q[2];
rz(2.8220693) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(-2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1324683) q[0];
sx q[0];
rz(-0.23129825) q[0];
sx q[0];
rz(-0.67434597) q[0];
rz(2.0293503) q[1];
sx q[1];
rz(-0.66450417) q[1];
sx q[1];
rz(-2.5792714) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.651197) q[0];
sx q[0];
rz(-2.2389452) q[0];
sx q[0];
rz(-0.018880318) q[0];
rz(-pi) q[1];
rz(2.855905) q[2];
sx q[2];
rz(-0.34491587) q[2];
sx q[2];
rz(3.0722741) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6866236) q[1];
sx q[1];
rz(-1.3775871) q[1];
sx q[1];
rz(1.4229694) q[1];
rz(-pi) q[2];
rz(-3.0543442) q[3];
sx q[3];
rz(-0.97594075) q[3];
sx q[3];
rz(-1.113254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.101863) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(-0.34004655) q[2];
rz(-2.9240821) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(2.8229009) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927004) q[0];
sx q[0];
rz(-3.0954439) q[0];
sx q[0];
rz(-0.39644077) q[0];
rz(-0.13892826) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(1.6202392) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8521261) q[0];
sx q[0];
rz(-0.20848256) q[0];
sx q[0];
rz(-0.33574386) q[0];
rz(0.87562008) q[2];
sx q[2];
rz(-2.5555829) q[2];
sx q[2];
rz(2.3919174) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.13008598) q[1];
sx q[1];
rz(-2.3475921) q[1];
sx q[1];
rz(-1.4873051) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89150724) q[3];
sx q[3];
rz(-1.659698) q[3];
sx q[3];
rz(0.31739435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56269318) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(-2.8472624) q[2];
rz(-1.1307905) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(2.14595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49333736) q[0];
sx q[0];
rz(-2.2264037) q[0];
sx q[0];
rz(-2.7822568) q[0];
rz(-0.94611478) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(-2.8709581) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04979241) q[0];
sx q[0];
rz(-1.6066178) q[0];
sx q[0];
rz(1.6019751) q[0];
rz(-pi) q[1];
rz(2.052202) q[2];
sx q[2];
rz(-2.319616) q[2];
sx q[2];
rz(-0.35441986) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.088085727) q[1];
sx q[1];
rz(-3.0778031) q[1];
sx q[1];
rz(-2.3699058) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98468303) q[3];
sx q[3];
rz(-0.95041785) q[3];
sx q[3];
rz(-1.8970722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3140807) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(-2.6861526) q[2];
rz(2.3296302) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6311326) q[0];
sx q[0];
rz(-1.6332508) q[0];
sx q[0];
rz(-0.73927885) q[0];
rz(-2.9108858) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(2.646692) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1586944) q[0];
sx q[0];
rz(-1.6086846) q[0];
sx q[0];
rz(-2.057103) q[0];
rz(2.8137389) q[2];
sx q[2];
rz(-1.5640537) q[2];
sx q[2];
rz(-1.1525796) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6634076) q[1];
sx q[1];
rz(-0.93485281) q[1];
sx q[1];
rz(-3.0031167) q[1];
x q[2];
rz(-3.0462618) q[3];
sx q[3];
rz(-2.2653326) q[3];
sx q[3];
rz(-1.3235843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0080002) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(0.32785329) q[2];
rz(0.19206364) q[3];
sx q[3];
rz(-2.8938507) q[3];
sx q[3];
rz(-1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(2.1223758) q[0];
sx q[0];
rz(-1.5126956) q[0];
sx q[0];
rz(1.4737286) q[0];
rz(-2.3090251) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(-2.9359948) q[2];
sx q[2];
rz(-1.0369221) q[2];
sx q[2];
rz(2.592929) q[2];
rz(-1.8264063) q[3];
sx q[3];
rz(-1.8288463) q[3];
sx q[3];
rz(2.0127206) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];