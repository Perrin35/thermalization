OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1420105) q[0];
sx q[0];
rz(-2.1394696) q[0];
sx q[0];
rz(0.89755091) q[0];
rz(-3.3759723) q[1];
sx q[1];
rz(3.4174089) q[1];
sx q[1];
rz(13.630907) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9165186) q[0];
sx q[0];
rz(-1.5009891) q[0];
sx q[0];
rz(0.9401456) q[0];
x q[1];
rz(0.93809442) q[2];
sx q[2];
rz(-1.2812867) q[2];
sx q[2];
rz(-3.0245568) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2805466) q[1];
sx q[1];
rz(-1.4006873) q[1];
sx q[1];
rz(0.64330805) q[1];
rz(-1.7552056) q[3];
sx q[3];
rz(-0.82114906) q[3];
sx q[3];
rz(2.8781995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.66427461) q[2];
sx q[2];
rz(-0.6330108) q[2];
sx q[2];
rz(-2.409639) q[2];
rz(0.96015635) q[3];
sx q[3];
rz(-0.82257661) q[3];
sx q[3];
rz(1.7094973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17523781) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(-0.91645855) q[0];
rz(-0.48049277) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(2.2629471) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4247596) q[0];
sx q[0];
rz(-0.75100198) q[0];
sx q[0];
rz(2.1483833) q[0];
rz(-1.8354561) q[2];
sx q[2];
rz(-0.96894962) q[2];
sx q[2];
rz(0.92483172) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0757054) q[1];
sx q[1];
rz(-1.9888708) q[1];
sx q[1];
rz(1.4904406) q[1];
rz(3.0733172) q[3];
sx q[3];
rz(-1.8797415) q[3];
sx q[3];
rz(-2.6451049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2800704) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(-2.4943165) q[2];
rz(0.17368008) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(2.98996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6140401) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(-2.8955984) q[0];
rz(1.7315158) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(-1.1211959) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58223984) q[0];
sx q[0];
rz(-1.9358578) q[0];
sx q[0];
rz(2.3854726) q[0];
rz(-pi) q[1];
rz(1.4175225) q[2];
sx q[2];
rz(-2.7445265) q[2];
sx q[2];
rz(3.0561471) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0908302) q[1];
sx q[1];
rz(-1.8659667) q[1];
sx q[1];
rz(-2.7961618) q[1];
rz(-pi) q[2];
rz(2.6477473) q[3];
sx q[3];
rz(-1.6742799) q[3];
sx q[3];
rz(-0.34085694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7401509) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(0.95139727) q[2];
rz(-0.65008632) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(-2.8459809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6275416) q[0];
sx q[0];
rz(-2.5294332) q[0];
sx q[0];
rz(-2.3535368) q[0];
rz(-2.0846941) q[1];
sx q[1];
rz(-1.1833271) q[1];
sx q[1];
rz(0.11638164) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7096036) q[0];
sx q[0];
rz(-0.63392144) q[0];
sx q[0];
rz(-3.0984127) q[0];
rz(-0.94159796) q[2];
sx q[2];
rz(-2.5430352) q[2];
sx q[2];
rz(1.7184005) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9434005) q[1];
sx q[1];
rz(-1.0775837) q[1];
sx q[1];
rz(2.0342159) q[1];
rz(-pi) q[2];
rz(2.651752) q[3];
sx q[3];
rz(-1.5252938) q[3];
sx q[3];
rz(-3.0385466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5300166) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(-2.8895565) q[2];
rz(-0.37825545) q[3];
sx q[3];
rz(-0.16246048) q[3];
sx q[3];
rz(0.3616412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5761121) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(-0.53260032) q[0];
rz(-1.700092) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(1.1486357) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.900433) q[0];
sx q[0];
rz(-2.1923089) q[0];
sx q[0];
rz(3.050699) q[0];
rz(0.45102851) q[2];
sx q[2];
rz(-0.94748679) q[2];
sx q[2];
rz(2.905068) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.83541218) q[1];
sx q[1];
rz(-2.7369943) q[1];
sx q[1];
rz(-0.49680357) q[1];
rz(-pi) q[2];
rz(0.43907459) q[3];
sx q[3];
rz(-0.77507654) q[3];
sx q[3];
rz(-0.72417688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0107515) q[2];
sx q[2];
rz(-2.4804219) q[2];
sx q[2];
rz(-1.032069) q[2];
rz(2.4268835) q[3];
sx q[3];
rz(-1.8604449) q[3];
sx q[3];
rz(-3.1177974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(-2.7456039) q[0];
rz(1.6962601) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(0.2063624) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2287184) q[0];
sx q[0];
rz(-1.657907) q[0];
sx q[0];
rz(-0.74761439) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80870734) q[2];
sx q[2];
rz(-0.84918298) q[2];
sx q[2];
rz(1.8546113) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5620835) q[1];
sx q[1];
rz(-1.675203) q[1];
sx q[1];
rz(-0.097192055) q[1];
x q[2];
rz(-2.0993125) q[3];
sx q[3];
rz(-2.7665666) q[3];
sx q[3];
rz(3.022775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6017194) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(0.27077857) q[2];
rz(0.74832908) q[3];
sx q[3];
rz(-1.7691408) q[3];
sx q[3];
rz(-1.07871) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8966184) q[0];
sx q[0];
rz(-2.0386319) q[0];
sx q[0];
rz(-0.57624972) q[0];
rz(-2.9684864) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(1.9304088) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9528708) q[0];
sx q[0];
rz(-1.0591918) q[0];
sx q[0];
rz(-2.1923724) q[0];
x q[1];
rz(1.7089825) q[2];
sx q[2];
rz(-2.0565363) q[2];
sx q[2];
rz(-1.1952343) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8163029) q[1];
sx q[1];
rz(-1.5029229) q[1];
sx q[1];
rz(-1.7273278) q[1];
rz(0.32254036) q[3];
sx q[3];
rz(-2.1493559) q[3];
sx q[3];
rz(1.6998147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8048191) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(3.1170735) q[2];
rz(2.426614) q[3];
sx q[3];
rz(-1.4638126) q[3];
sx q[3];
rz(1.582675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1361168) q[0];
sx q[0];
rz(-1.5908717) q[0];
sx q[0];
rz(2.1210282) q[0];
rz(-2.9868946) q[1];
sx q[1];
rz(-1.5203412) q[1];
sx q[1];
rz(1.221009) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97604254) q[0];
sx q[0];
rz(-1.8414458) q[0];
sx q[0];
rz(-2.3483777) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9646655) q[2];
sx q[2];
rz(-1.7311586) q[2];
sx q[2];
rz(-2.4307876) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.029433) q[1];
sx q[1];
rz(-0.1754079) q[1];
sx q[1];
rz(-0.68316858) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8052748) q[3];
sx q[3];
rz(-2.291403) q[3];
sx q[3];
rz(-1.99828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8528379) q[2];
sx q[2];
rz(-1.6483665) q[2];
sx q[2];
rz(0.70518804) q[2];
rz(2.5382036) q[3];
sx q[3];
rz(-2.2796977) q[3];
sx q[3];
rz(1.5195297) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0309546) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(3.0045793) q[0];
rz(-0.6048454) q[1];
sx q[1];
rz(-0.73692656) q[1];
sx q[1];
rz(0.12577122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9916358) q[0];
sx q[0];
rz(-1.7466674) q[0];
sx q[0];
rz(-1.7624904) q[0];
x q[1];
rz(0.33340402) q[2];
sx q[2];
rz(-1.3840904) q[2];
sx q[2];
rz(-0.29078996) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6313435) q[1];
sx q[1];
rz(-1.2878839) q[1];
sx q[1];
rz(0.81088539) q[1];
rz(-pi) q[2];
x q[2];
rz(2.845876) q[3];
sx q[3];
rz(-1.3774646) q[3];
sx q[3];
rz(-2.8545692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.003309) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(2.4323145) q[2];
rz(2.5214031) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(-2.9848849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.9897292) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(1.2003157) q[0];
rz(0.26750803) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(1.1402003) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88975554) q[0];
sx q[0];
rz(-1.6701084) q[0];
sx q[0];
rz(-1.3072144) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26009772) q[2];
sx q[2];
rz(-1.1220699) q[2];
sx q[2];
rz(1.5990433) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4789341) q[1];
sx q[1];
rz(-1.7591898) q[1];
sx q[1];
rz(1.6608095) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8285698) q[3];
sx q[3];
rz(-0.81837294) q[3];
sx q[3];
rz(2.5767874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.37529477) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(2.2429788) q[2];
rz(-3.1344154) q[3];
sx q[3];
rz(-2.3982748) q[3];
sx q[3];
rz(1.7858508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89467775) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(2.6208411) q[1];
sx q[1];
rz(-1.755935) q[1];
sx q[1];
rz(-1.4204949) q[1];
rz(2.3788135) q[2];
sx q[2];
rz(-0.63763466) q[2];
sx q[2];
rz(0.55479738) q[2];
rz(-1.981001) q[3];
sx q[3];
rz(-1.0241749) q[3];
sx q[3];
rz(0.54815626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
