OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.99958217) q[0];
sx q[0];
rz(-1.002123) q[0];
sx q[0];
rz(2.2440417) q[0];
rz(-3.3759723) q[1];
sx q[1];
rz(3.4174089) q[1];
sx q[1];
rz(13.630907) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8911154) q[0];
sx q[0];
rz(-2.5076137) q[0];
sx q[0];
rz(1.4527713) q[0];
rz(-pi) q[1];
rz(-1.1041553) q[2];
sx q[2];
rz(-2.4541514) q[2];
sx q[2];
rz(-1.825037) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0734288) q[1];
sx q[1];
rz(-2.4792719) q[1];
sx q[1];
rz(0.27889241) q[1];
rz(-pi) q[2];
rz(2.3834121) q[3];
sx q[3];
rz(-1.7054134) q[3];
sx q[3];
rz(-1.7077703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.477318) q[2];
sx q[2];
rz(-0.6330108) q[2];
sx q[2];
rz(0.73195362) q[2];
rz(0.96015635) q[3];
sx q[3];
rz(-0.82257661) q[3];
sx q[3];
rz(1.7094973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.5669211) q[1];
sx q[1];
rz(-2.2629471) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6966349) q[0];
sx q[0];
rz(-0.96224552) q[0];
sx q[0];
rz(-2.6702325) q[0];
rz(-pi) q[1];
rz(-0.36388134) q[2];
sx q[2];
rz(-0.65081396) q[2];
sx q[2];
rz(1.770307) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6138184) q[1];
sx q[1];
rz(-1.6442181) q[1];
sx q[1];
rz(0.41927494) q[1];
rz(-pi) q[2];
rz(1.8804178) q[3];
sx q[3];
rz(-1.6358346) q[3];
sx q[3];
rz(-1.0950973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2800704) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(2.4943165) q[2];
rz(0.17368008) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(2.98996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52755255) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(-0.24599427) q[0];
rz(1.4100769) q[1];
sx q[1];
rz(-1.9672111) q[1];
sx q[1];
rz(2.0203967) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58223984) q[0];
sx q[0];
rz(-1.2057349) q[0];
sx q[0];
rz(-0.75612005) q[0];
rz(-1.1779184) q[2];
sx q[2];
rz(-1.5117206) q[2];
sx q[2];
rz(1.3438366) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0908302) q[1];
sx q[1];
rz(-1.275626) q[1];
sx q[1];
rz(-0.34543085) q[1];
x q[2];
rz(-1.6882012) q[3];
sx q[3];
rz(-1.0798287) q[3];
sx q[3];
rz(1.1743869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7401509) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(-0.95139727) q[2];
rz(2.4915063) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51405108) q[0];
sx q[0];
rz(-2.5294332) q[0];
sx q[0];
rz(-2.3535368) q[0];
rz(1.0568985) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(3.025211) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1040092) q[0];
sx q[0];
rz(-1.596367) q[0];
sx q[0];
rz(0.63347647) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0747244) q[2];
sx q[2];
rz(-1.2328086) q[2];
sx q[2];
rz(2.4525814) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1310281) q[1];
sx q[1];
rz(-0.66337913) q[1];
sx q[1];
rz(0.69372155) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4898407) q[3];
sx q[3];
rz(-1.6162989) q[3];
sx q[3];
rz(0.10304606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6115761) q[2];
sx q[2];
rz(-1.896603) q[2];
sx q[2];
rz(2.8895565) q[2];
rz(-0.37825545) q[3];
sx q[3];
rz(-0.16246048) q[3];
sx q[3];
rz(-2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5761121) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(0.53260032) q[0];
rz(1.700092) q[1];
sx q[1];
rz(-0.36591995) q[1];
sx q[1];
rz(-1.9929569) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38265739) q[0];
sx q[0];
rz(-1.4969345) q[0];
sx q[0];
rz(2.1942684) q[0];
rz(-pi) q[1];
rz(-0.89678905) q[2];
sx q[2];
rz(-1.2090346) q[2];
sx q[2];
rz(-1.0587495) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3061805) q[1];
sx q[1];
rz(-2.7369943) q[1];
sx q[1];
rz(0.49680357) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43907459) q[3];
sx q[3];
rz(-2.3665161) q[3];
sx q[3];
rz(-2.4174158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0107515) q[2];
sx q[2];
rz(-2.4804219) q[2];
sx q[2];
rz(1.032069) q[2];
rz(-2.4268835) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(-3.1177974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5500568) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(-2.7456039) q[0];
rz(-1.6962601) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(2.9352303) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43564046) q[0];
sx q[0];
rz(-0.75169509) q[0];
sx q[0];
rz(-0.12775001) q[0];
x q[1];
rz(-0.80870734) q[2];
sx q[2];
rz(-2.2924097) q[2];
sx q[2];
rz(-1.8546113) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.57950912) q[1];
sx q[1];
rz(-1.675203) q[1];
sx q[1];
rz(3.0444006) q[1];
x q[2];
rz(-1.0422802) q[3];
sx q[3];
rz(-2.7665666) q[3];
sx q[3];
rz(0.11881766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6017194) q[2];
sx q[2];
rz(-2.0809934) q[2];
sx q[2];
rz(2.8708141) q[2];
rz(2.3932636) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(-1.07871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2449743) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(0.57624972) q[0];
rz(-0.17310625) q[1];
sx q[1];
rz(-2.3997967) q[1];
sx q[1];
rz(1.9304088) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21748397) q[0];
sx q[0];
rz(-2.3586914) q[0];
sx q[0];
rz(-2.3379675) q[0];
rz(-pi) q[1];
rz(-1.7089825) q[2];
sx q[2];
rz(-1.0850564) q[2];
sx q[2];
rz(1.9463584) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8163029) q[1];
sx q[1];
rz(-1.6386697) q[1];
sx q[1];
rz(-1.7273278) q[1];
rz(-pi) q[2];
rz(2.022642) q[3];
sx q[3];
rz(-0.65330905) q[3];
sx q[3];
rz(2.2484231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8048191) q[2];
sx q[2];
rz(-2.4833198) q[2];
sx q[2];
rz(0.024519196) q[2];
rz(-0.71497861) q[3];
sx q[3];
rz(-1.4638126) q[3];
sx q[3];
rz(-1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(3.1361168) q[0];
sx q[0];
rz(-1.5908717) q[0];
sx q[0];
rz(-2.1210282) q[0];
rz(-2.9868946) q[1];
sx q[1];
rz(-1.5203412) q[1];
sx q[1];
rz(-1.9205836) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85212612) q[0];
sx q[0];
rz(-2.3131436) q[0];
sx q[0];
rz(-0.37129398) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74322015) q[2];
sx q[2];
rz(-0.23822242) q[2];
sx q[2];
rz(1.5526349) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.029433) q[1];
sx q[1];
rz(-0.1754079) q[1];
sx q[1];
rz(-0.68316858) q[1];
rz(-pi) q[2];
rz(0.82151316) q[3];
sx q[3];
rz(-1.8213846) q[3];
sx q[3];
rz(0.65419765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8528379) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(-0.70518804) q[2];
rz(0.60338902) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(-1.6220629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.110638) q[0];
sx q[0];
rz(-1.5663261) q[0];
sx q[0];
rz(-3.0045793) q[0];
rz(-2.5367472) q[1];
sx q[1];
rz(-0.73692656) q[1];
sx q[1];
rz(-0.12577122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14995689) q[0];
sx q[0];
rz(-1.3949252) q[0];
sx q[0];
rz(1.7624904) q[0];
x q[1];
rz(-0.33340402) q[2];
sx q[2];
rz(-1.3840904) q[2];
sx q[2];
rz(-2.8508027) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3667664) q[1];
sx q[1];
rz(-0.80087304) q[1];
sx q[1];
rz(-1.9701387) q[1];
rz(-2.5500507) q[3];
sx q[3];
rz(-2.7898443) q[3];
sx q[3];
rz(1.2951617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13828364) q[2];
sx q[2];
rz(-0.97390276) q[2];
sx q[2];
rz(0.70927817) q[2];
rz(2.5214031) q[3];
sx q[3];
rz(-0.87738335) q[3];
sx q[3];
rz(-0.15670776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1518635) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(-1.2003157) q[0];
rz(2.8740846) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(-1.1402003) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2518371) q[0];
sx q[0];
rz(-1.6701084) q[0];
sx q[0];
rz(-1.3072144) q[0];
rz(-1.0802202) q[2];
sx q[2];
rz(-2.6274101) q[2];
sx q[2];
rz(-2.0928004) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4789341) q[1];
sx q[1];
rz(-1.3824029) q[1];
sx q[1];
rz(-1.6608095) q[1];
x q[2];
rz(-0.31302281) q[3];
sx q[3];
rz(-0.81837294) q[3];
sx q[3];
rz(-0.5648053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37529477) q[2];
sx q[2];
rz(-1.4582429) q[2];
sx q[2];
rz(-2.2429788) q[2];
rz(0.0071772655) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(-1.7858508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2469149) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(0.5207516) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(-2.3788135) q[2];
sx q[2];
rz(-2.503958) q[2];
sx q[2];
rz(-2.5867953) q[2];
rz(1.1605916) q[3];
sx q[3];
rz(-1.0241749) q[3];
sx q[3];
rz(0.54815626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
