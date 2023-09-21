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
rz(5.2810623) q[0];
sx q[0];
rz(5.3856344) q[0];
rz(-3.3759723) q[1];
sx q[1];
rz(3.4174089) q[1];
sx q[1];
rz(13.630907) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8911154) q[0];
sx q[0];
rz(-0.63397898) q[0];
sx q[0];
rz(-1.4527713) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93809442) q[2];
sx q[2];
rz(-1.860306) q[2];
sx q[2];
rz(0.11703581) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.068163888) q[1];
sx q[1];
rz(-0.66232077) q[1];
sx q[1];
rz(0.27889241) q[1];
rz(1.3863871) q[3];
sx q[3];
rz(-2.3204436) q[3];
sx q[3];
rz(0.26339312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.66427461) q[2];
sx q[2];
rz(-0.6330108) q[2];
sx q[2];
rz(-0.73195362) q[2];
rz(2.1814363) q[3];
sx q[3];
rz(-0.82257661) q[3];
sx q[3];
rz(1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17523781) q[0];
sx q[0];
rz(-2.2580999) q[0];
sx q[0];
rz(2.2251341) q[0];
rz(0.48049277) q[1];
sx q[1];
rz(-2.5669211) q[1];
sx q[1];
rz(-0.8786456) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7322313) q[0];
sx q[0];
rz(-1.9525813) q[0];
sx q[0];
rz(0.9071) q[0];
rz(2.7777113) q[2];
sx q[2];
rz(-2.4907787) q[2];
sx q[2];
rz(-1.770307) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52777427) q[1];
sx q[1];
rz(-1.6442181) q[1];
sx q[1];
rz(0.41927494) q[1];
rz(-pi) q[2];
rz(-0.068275498) q[3];
sx q[3];
rz(-1.8797415) q[3];
sx q[3];
rz(0.49648778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2800704) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(0.64727616) q[2];
rz(0.17368008) q[3];
sx q[3];
rz(-2.0208385) q[3];
sx q[3];
rz(-2.98996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52755255) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(2.8955984) q[0];
rz(1.7315158) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(2.0203967) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66378731) q[0];
sx q[0];
rz(-0.87513798) q[0];
sx q[0];
rz(2.0545161) q[0];
x q[1];
rz(3.0776575) q[2];
sx q[2];
rz(-1.1786412) q[2];
sx q[2];
rz(2.8901697) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.4157297) q[1];
sx q[1];
rz(-1.9007069) q[1];
sx q[1];
rz(1.8833453) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4533914) q[3];
sx q[3];
rz(-2.061764) q[3];
sx q[3];
rz(-1.1743869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40144172) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(0.95139727) q[2];
rz(-2.4915063) q[3];
sx q[3];
rz(-1.254046) q[3];
sx q[3];
rz(0.29561177) q[3];
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
rz(pi/2) q[0];
x q[1];
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
rz(2.0846941) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(-3.025211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1040092) q[0];
sx q[0];
rz(-1.5452256) q[0];
sx q[0];
rz(0.63347647) q[0];
rz(2.1999947) q[2];
sx q[2];
rz(-2.5430352) q[2];
sx q[2];
rz(1.7184005) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0013106) q[1];
sx q[1];
rz(-1.9754859) q[1];
sx q[1];
rz(-2.6005122) q[1];
rz(0.096480358) q[3];
sx q[3];
rz(-2.6498142) q[3];
sx q[3];
rz(1.5887367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5300166) q[2];
sx q[2];
rz(-1.896603) q[2];
sx q[2];
rz(2.8895565) q[2];
rz(2.7633372) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(2.7799515) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5761121) q[0];
sx q[0];
rz(-1.7385372) q[0];
sx q[0];
rz(-0.53260032) q[0];
rz(-1.4415007) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(1.9929569) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.900433) q[0];
sx q[0];
rz(-2.1923089) q[0];
sx q[0];
rz(0.09089367) q[0];
rz(-pi) q[1];
rz(-2.6905641) q[2];
sx q[2];
rz(-0.94748679) q[2];
sx q[2];
rz(-0.23652467) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7733113) q[1];
sx q[1];
rz(-1.2174264) q[1];
sx q[1];
rz(1.7721304) q[1];
x q[2];
rz(0.72539056) q[3];
sx q[3];
rz(-1.2687506) q[3];
sx q[3];
rz(0.52291742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1308412) q[2];
sx q[2];
rz(-2.4804219) q[2];
sx q[2];
rz(1.032069) q[2];
rz(2.4268835) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(3.1177974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(-2.7456039) q[0];
rz(-1.6962601) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(-0.2063624) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7059522) q[0];
sx q[0];
rz(-2.3898976) q[0];
sx q[0];
rz(0.12775001) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3328853) q[2];
sx q[2];
rz(-2.2924097) q[2];
sx q[2];
rz(1.2869814) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1604662) q[1];
sx q[1];
rz(-1.6674575) q[1];
sx q[1];
rz(-1.4658982) q[1];
rz(0.19595512) q[3];
sx q[3];
rz(-1.8926419) q[3];
sx q[3];
rz(-0.67924196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.53987327) q[2];
sx q[2];
rz(-2.0809934) q[2];
sx q[2];
rz(2.8708141) q[2];
rz(-2.3932636) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(-2.0628827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2449743) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(2.5653429) q[0];
rz(2.9684864) q[1];
sx q[1];
rz(-2.3997967) q[1];
sx q[1];
rz(1.9304088) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71935463) q[0];
sx q[0];
rz(-2.1033759) q[0];
sx q[0];
rz(-0.60441916) q[0];
x q[1];
rz(1.4326101) q[2];
sx q[2];
rz(-1.0850564) q[2];
sx q[2];
rz(-1.1952343) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.83963517) q[1];
sx q[1];
rz(-2.9710899) q[1];
sx q[1];
rz(1.9819928) q[1];
rz(2.1738449) q[3];
sx q[3];
rz(-1.3021886) q[3];
sx q[3];
rz(0.3097765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8048191) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(3.1170735) q[2];
rz(-0.71497861) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0054758469) q[0];
sx q[0];
rz(-1.5908717) q[0];
sx q[0];
rz(1.0205644) q[0];
rz(0.15469805) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(1.9205836) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97604254) q[0];
sx q[0];
rz(-1.8414458) q[0];
sx q[0];
rz(0.793215) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3983725) q[2];
sx q[2];
rz(-0.23822242) q[2];
sx q[2];
rz(1.5526349) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.029433) q[1];
sx q[1];
rz(-0.1754079) q[1];
sx q[1];
rz(-2.4584241) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.28875479) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(-0.70518804) q[2];
rz(-0.60338902) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(-1.5195297) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0309546) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(-0.13701339) q[0];
rz(0.6048454) q[1];
sx q[1];
rz(-0.73692656) q[1];
sx q[1];
rz(-0.12577122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14995689) q[0];
sx q[0];
rz(-1.3949252) q[0];
sx q[0];
rz(1.3791023) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7681098) q[2];
sx q[2];
rz(-1.8981877) q[2];
sx q[2];
rz(-1.9257853) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6313435) q[1];
sx q[1];
rz(-1.2878839) q[1];
sx q[1];
rz(2.3307073) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7726692) q[3];
sx q[3];
rz(-1.8608421) q[3];
sx q[3];
rz(1.2253075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13828364) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(0.70927817) q[2];
rz(-2.5214031) q[3];
sx q[3];
rz(-0.87738335) q[3];
sx q[3];
rz(0.15670776) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1518635) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(-1.2003157) q[0];
rz(0.26750803) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(1.1402003) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2518371) q[0];
sx q[0];
rz(-1.6701084) q[0];
sx q[0];
rz(1.3072144) q[0];
rz(-pi) q[1];
rz(-2.033038) q[2];
sx q[2];
rz(-1.8046364) q[2];
sx q[2];
rz(-3.0548981) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4789341) q[1];
sx q[1];
rz(-1.7591898) q[1];
sx q[1];
rz(-1.4807832) q[1];
rz(2.3481028) q[3];
sx q[3];
rz(-1.7975382) q[3];
sx q[3];
rz(-1.2236809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.37529477) q[2];
sx q[2];
rz(-1.4582429) q[2];
sx q[2];
rz(2.2429788) q[2];
rz(-0.0071772655) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(-1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2469149) q[0];
sx q[0];
rz(-1.4151731) q[0];
sx q[0];
rz(-1.7807501) q[0];
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
rz(0.58017147) q[3];
sx q[3];
rz(-0.67065722) q[3];
sx q[3];
rz(-1.8967659) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
