OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95076686) q[0];
sx q[0];
rz(-1.8523676) q[0];
sx q[0];
rz(-3.0486795) q[0];
rz(-0.095280401) q[1];
sx q[1];
rz(-0.73268259) q[1];
sx q[1];
rz(-1.2394152) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6042955) q[0];
sx q[0];
rz(-1.822246) q[0];
sx q[0];
rz(0.23668134) q[0];
rz(-1.3746512) q[2];
sx q[2];
rz(-1.8298379) q[2];
sx q[2];
rz(2.1262622) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.79151407) q[1];
sx q[1];
rz(-0.33153141) q[1];
sx q[1];
rz(3.0505807) q[1];
x q[2];
rz(2.0650571) q[3];
sx q[3];
rz(-2.0856033) q[3];
sx q[3];
rz(-0.47806397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6713082) q[2];
sx q[2];
rz(-1.6613864) q[2];
sx q[2];
rz(-1.7926463) q[2];
rz(0.74364439) q[3];
sx q[3];
rz(-2.8642004) q[3];
sx q[3];
rz(-2.9644137) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.141356) q[0];
sx q[0];
rz(-1.5445222) q[0];
sx q[0];
rz(2.8479688) q[0];
rz(2.1108421) q[1];
sx q[1];
rz(-1.7799957) q[1];
sx q[1];
rz(0.34293276) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4727508) q[0];
sx q[0];
rz(-2.0905417) q[0];
sx q[0];
rz(-0.23730554) q[0];
rz(-1.9428217) q[2];
sx q[2];
rz(-1.4497533) q[2];
sx q[2];
rz(-0.72265676) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3174177) q[1];
sx q[1];
rz(-1.7848099) q[1];
sx q[1];
rz(-0.70862464) q[1];
rz(-0.94128709) q[3];
sx q[3];
rz(-1.3664748) q[3];
sx q[3];
rz(-2.2245537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.12884101) q[2];
sx q[2];
rz(-1.7495456) q[2];
sx q[2];
rz(0.7592321) q[2];
rz(0.020708474) q[3];
sx q[3];
rz(-1.8755251) q[3];
sx q[3];
rz(0.43829632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6193806) q[0];
sx q[0];
rz(-2.5396357) q[0];
sx q[0];
rz(-0.0044599175) q[0];
rz(2.4941173) q[1];
sx q[1];
rz(-0.59440333) q[1];
sx q[1];
rz(2.3562145) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68092184) q[0];
sx q[0];
rz(-0.17853949) q[0];
sx q[0];
rz(2.6723249) q[0];
x q[1];
rz(2.0702281) q[2];
sx q[2];
rz(-2.7874261) q[2];
sx q[2];
rz(0.11812299) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0177893) q[1];
sx q[1];
rz(-0.54517704) q[1];
sx q[1];
rz(1.7032318) q[1];
rz(-2.4786199) q[3];
sx q[3];
rz(-2.1372652) q[3];
sx q[3];
rz(0.86323767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.37041) q[2];
sx q[2];
rz(-2.2213171) q[2];
sx q[2];
rz(1.3578337) q[2];
rz(-2.8010098) q[3];
sx q[3];
rz(-2.1850696) q[3];
sx q[3];
rz(2.3497439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1444645) q[0];
sx q[0];
rz(-2.0576532) q[0];
sx q[0];
rz(0.88515627) q[0];
rz(0.74686933) q[1];
sx q[1];
rz(-0.62364686) q[1];
sx q[1];
rz(-2.0451827) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089450739) q[0];
sx q[0];
rz(-0.47596395) q[0];
sx q[0];
rz(2.3030998) q[0];
rz(0.030641067) q[2];
sx q[2];
rz(-1.6302135) q[2];
sx q[2];
rz(-0.63842809) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1433574) q[1];
sx q[1];
rz(-1.1117742) q[1];
sx q[1];
rz(-1.281339) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73170264) q[3];
sx q[3];
rz(-1.9644794) q[3];
sx q[3];
rz(-1.3418152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5235644) q[2];
sx q[2];
rz(-2.9310493) q[2];
sx q[2];
rz(-3.1169685) q[2];
rz(-0.63557449) q[3];
sx q[3];
rz(-2.1533951) q[3];
sx q[3];
rz(-2.2583101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4516975) q[0];
sx q[0];
rz(-2.8391333) q[0];
sx q[0];
rz(0.46689335) q[0];
rz(-0.11058552) q[1];
sx q[1];
rz(-2.6778335) q[1];
sx q[1];
rz(-1.0708403) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12884451) q[0];
sx q[0];
rz(-1.8086047) q[0];
sx q[0];
rz(-2.0369916) q[0];
rz(-2.6314965) q[2];
sx q[2];
rz(-0.3613216) q[2];
sx q[2];
rz(0.23621836) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3339798) q[1];
sx q[1];
rz(-2.4845504) q[1];
sx q[1];
rz(0.64589898) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4399906) q[3];
sx q[3];
rz(-1.0402586) q[3];
sx q[3];
rz(2.4367133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0232627) q[2];
sx q[2];
rz(-1.7795965) q[2];
sx q[2];
rz(2.2892717) q[2];
rz(3.1039589) q[3];
sx q[3];
rz(-1.2331542) q[3];
sx q[3];
rz(0.5948624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1325876) q[0];
sx q[0];
rz(-1.1786893) q[0];
sx q[0];
rz(-1.8527385) q[0];
rz(-1.6291078) q[1];
sx q[1];
rz(-1.1237203) q[1];
sx q[1];
rz(1.9992794) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61296755) q[0];
sx q[0];
rz(-0.36699793) q[0];
sx q[0];
rz(-0.27530833) q[0];
rz(1.3630168) q[2];
sx q[2];
rz(-1.5099974) q[2];
sx q[2];
rz(-2.6039993) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.21396046) q[1];
sx q[1];
rz(-0.42812706) q[1];
sx q[1];
rz(0.63251782) q[1];
x q[2];
rz(1.7473502) q[3];
sx q[3];
rz(-1.8774967) q[3];
sx q[3];
rz(0.97489417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8196572) q[2];
sx q[2];
rz(-1.2677931) q[2];
sx q[2];
rz(0.12040559) q[2];
rz(-1.985792) q[3];
sx q[3];
rz(-1.5294231) q[3];
sx q[3];
rz(-2.9175478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8026546) q[0];
sx q[0];
rz(-1.8010362) q[0];
sx q[0];
rz(-1.4539723) q[0];
rz(1.5441719) q[1];
sx q[1];
rz(-1.495196) q[1];
sx q[1];
rz(2.6905751) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4225578) q[0];
sx q[0];
rz(-1.8717878) q[0];
sx q[0];
rz(1.5084356) q[0];
rz(-1.9372378) q[2];
sx q[2];
rz(-1.7882573) q[2];
sx q[2];
rz(-1.6161473) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6587131) q[1];
sx q[1];
rz(-0.45156839) q[1];
sx q[1];
rz(-2.8612479) q[1];
x q[2];
rz(2.3312206) q[3];
sx q[3];
rz(-1.5640508) q[3];
sx q[3];
rz(0.41616671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.94577998) q[2];
sx q[2];
rz(-1.4142298) q[2];
sx q[2];
rz(-1.3607402) q[2];
rz(-1.0055379) q[3];
sx q[3];
rz(-1.1755627) q[3];
sx q[3];
rz(-1.3895234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0176395) q[0];
sx q[0];
rz(-3.0160976) q[0];
sx q[0];
rz(-2.4355167) q[0];
rz(1.5438682) q[1];
sx q[1];
rz(-1.7192625) q[1];
sx q[1];
rz(2.2854038) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5498915) q[0];
sx q[0];
rz(-3.1208) q[0];
sx q[0];
rz(-1.2077622) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1107619) q[2];
sx q[2];
rz(-0.80321124) q[2];
sx q[2];
rz(0.75424657) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.44196196) q[1];
sx q[1];
rz(-2.5285427) q[1];
sx q[1];
rz(-2.7837201) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.007902) q[3];
sx q[3];
rz(-1.0858337) q[3];
sx q[3];
rz(-0.07459379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2989444) q[2];
sx q[2];
rz(-1.7405258) q[2];
sx q[2];
rz(-2.974158) q[2];
rz(1.5873448) q[3];
sx q[3];
rz(-0.7730248) q[3];
sx q[3];
rz(-2.1412444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3779959) q[0];
sx q[0];
rz(-1.6908228) q[0];
sx q[0];
rz(-2.7698621) q[0];
rz(-1.4338214) q[1];
sx q[1];
rz(-1.5377518) q[1];
sx q[1];
rz(-2.8415714) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8512631) q[0];
sx q[0];
rz(-1.4545049) q[0];
sx q[0];
rz(2.854748) q[0];
rz(1.6113564) q[2];
sx q[2];
rz(-0.99405655) q[2];
sx q[2];
rz(0.64740136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.91050657) q[1];
sx q[1];
rz(-2.4795737) q[1];
sx q[1];
rz(0.20259095) q[1];
x q[2];
rz(0.12436538) q[3];
sx q[3];
rz(-0.87166407) q[3];
sx q[3];
rz(-1.82064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6649449) q[2];
sx q[2];
rz(-2.6740394) q[2];
sx q[2];
rz(-0.42759582) q[2];
rz(-1.556373) q[3];
sx q[3];
rz(-2.0245602) q[3];
sx q[3];
rz(-0.75470406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4204191) q[0];
sx q[0];
rz(-0.54062802) q[0];
sx q[0];
rz(-0.46947259) q[0];
rz(-2.2162614) q[1];
sx q[1];
rz(-2.053849) q[1];
sx q[1];
rz(3.1184149) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1184922) q[0];
sx q[0];
rz(-1.4064624) q[0];
sx q[0];
rz(-0.53043764) q[0];
x q[1];
rz(2.2830354) q[2];
sx q[2];
rz(-0.80712592) q[2];
sx q[2];
rz(1.6940534) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.35649037) q[1];
sx q[1];
rz(-1.7447326) q[1];
sx q[1];
rz(2.9746303) q[1];
x q[2];
rz(-2.3007727) q[3];
sx q[3];
rz(-1.7125907) q[3];
sx q[3];
rz(0.26103448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2716486) q[2];
sx q[2];
rz(-1.2056489) q[2];
sx q[2];
rz(-2.6507157) q[2];
rz(-2.0513746) q[3];
sx q[3];
rz(-2.1722983) q[3];
sx q[3];
rz(-0.31392613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28521095) q[0];
sx q[0];
rz(-2.2638392) q[0];
sx q[0];
rz(1.0765156) q[0];
rz(-2.8648227) q[1];
sx q[1];
rz(-2.3634187) q[1];
sx q[1];
rz(1.707911) q[1];
rz(-1.8770915) q[2];
sx q[2];
rz(-2.2363792) q[2];
sx q[2];
rz(-2.1909942) q[2];
rz(2.6071045) q[3];
sx q[3];
rz(-1.5535003) q[3];
sx q[3];
rz(1.5737892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
