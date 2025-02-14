OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8283451) q[0];
sx q[0];
rz(-0.77594835) q[0];
sx q[0];
rz(2.0394072) q[0];
rz(1.6232396) q[1];
sx q[1];
rz(2.0145388) q[1];
sx q[1];
rz(9.1993499) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0505053) q[0];
sx q[0];
rz(-2.4417672) q[0];
sx q[0];
rz(-3.0449633) q[0];
x q[1];
rz(-0.20689865) q[2];
sx q[2];
rz(-0.62600905) q[2];
sx q[2];
rz(-1.8459143) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.082583383) q[1];
sx q[1];
rz(-1.2066926) q[1];
sx q[1];
rz(2.8541982) q[1];
x q[2];
rz(2.647764) q[3];
sx q[3];
rz(-1.2612311) q[3];
sx q[3];
rz(1.0012331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.48308358) q[2];
sx q[2];
rz(-2.4675214) q[2];
sx q[2];
rz(0.020641208) q[2];
rz(-2.9523197) q[3];
sx q[3];
rz(-2.7920189) q[3];
sx q[3];
rz(-1.4741723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3610158) q[0];
sx q[0];
rz(-0.78106946) q[0];
sx q[0];
rz(-2.1625157) q[0];
rz(3.0304404) q[1];
sx q[1];
rz(-2.4039098) q[1];
sx q[1];
rz(2.0806064) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8421335) q[0];
sx q[0];
rz(-0.7086904) q[0];
sx q[0];
rz(0.75229074) q[0];
x q[1];
rz(-1.5699638) q[2];
sx q[2];
rz(-1.4333409) q[2];
sx q[2];
rz(0.47333131) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0506405) q[1];
sx q[1];
rz(-2.0555858) q[1];
sx q[1];
rz(1.3371972) q[1];
rz(-0.33943601) q[3];
sx q[3];
rz(-2.7314679) q[3];
sx q[3];
rz(0.29537485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7079118) q[2];
sx q[2];
rz(-1.3206864) q[2];
sx q[2];
rz(0.39805463) q[2];
rz(-1.4237283) q[3];
sx q[3];
rz(-0.28702304) q[3];
sx q[3];
rz(0.48582336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61755359) q[0];
sx q[0];
rz(-2.6300639) q[0];
sx q[0];
rz(2.0892573) q[0];
rz(0.44318336) q[1];
sx q[1];
rz(-0.96142238) q[1];
sx q[1];
rz(2.4876432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67414647) q[0];
sx q[0];
rz(-2.3096414) q[0];
sx q[0];
rz(1.7427) q[0];
rz(-pi) q[1];
rz(-2.1939799) q[2];
sx q[2];
rz(-1.6296367) q[2];
sx q[2];
rz(-0.72966444) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3824275) q[1];
sx q[1];
rz(-2.2780721) q[1];
sx q[1];
rz(2.1334126) q[1];
rz(-1.3929358) q[3];
sx q[3];
rz(-1.3452969) q[3];
sx q[3];
rz(-2.5009843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7707278) q[2];
sx q[2];
rz(-2.3720522) q[2];
sx q[2];
rz(-0.23180836) q[2];
rz(1.2109463) q[3];
sx q[3];
rz(-1.1823267) q[3];
sx q[3];
rz(-0.32538357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66832191) q[0];
sx q[0];
rz(-0.40437651) q[0];
sx q[0];
rz(-2.7647198) q[0];
rz(-2.6301774) q[1];
sx q[1];
rz(-1.2579974) q[1];
sx q[1];
rz(-2.2976141) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91746861) q[0];
sx q[0];
rz(-3.030405) q[0];
sx q[0];
rz(1.412941) q[0];
x q[1];
rz(1.8966214) q[2];
sx q[2];
rz(-2.5375536) q[2];
sx q[2];
rz(-2.2789795) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6214024) q[1];
sx q[1];
rz(-2.108413) q[1];
sx q[1];
rz(-2.138184) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52167251) q[3];
sx q[3];
rz(-1.552201) q[3];
sx q[3];
rz(-2.2357495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1692928) q[2];
sx q[2];
rz(-1.0147213) q[2];
sx q[2];
rz(-0.69804066) q[2];
rz(1.9418779) q[3];
sx q[3];
rz(-2.2359087) q[3];
sx q[3];
rz(-0.5996632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19071628) q[0];
sx q[0];
rz(-2.8671725) q[0];
sx q[0];
rz(3.0158667) q[0];
rz(2.114864) q[1];
sx q[1];
rz(-1.7544361) q[1];
sx q[1];
rz(0.52621192) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3644117) q[0];
sx q[0];
rz(-2.6573349) q[0];
sx q[0];
rz(-0.15987349) q[0];
rz(-pi) q[1];
rz(1.7650928) q[2];
sx q[2];
rz(-0.70617968) q[2];
sx q[2];
rz(-2.6841109) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.80568589) q[1];
sx q[1];
rz(-1.9503924) q[1];
sx q[1];
rz(-1.4604881) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86889699) q[3];
sx q[3];
rz(-0.90889895) q[3];
sx q[3];
rz(-1.6493662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0478263) q[2];
sx q[2];
rz(-0.83790773) q[2];
sx q[2];
rz(2.8313336) q[2];
rz(-3.1066762) q[3];
sx q[3];
rz(-1.3860476) q[3];
sx q[3];
rz(-2.3518899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0763615) q[0];
sx q[0];
rz(-1.4820453) q[0];
sx q[0];
rz(0.37102997) q[0];
rz(-3.0767483) q[1];
sx q[1];
rz(-0.82227451) q[1];
sx q[1];
rz(-1.8745905) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2081535) q[0];
sx q[0];
rz(-1.9133718) q[0];
sx q[0];
rz(-0.88985635) q[0];
rz(-pi) q[1];
rz(0.94159884) q[2];
sx q[2];
rz(-1.8591188) q[2];
sx q[2];
rz(1.1278111) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1888926) q[1];
sx q[1];
rz(-1.6052393) q[1];
sx q[1];
rz(-0.39826213) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46601136) q[3];
sx q[3];
rz(-1.8309085) q[3];
sx q[3];
rz(-2.8046908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.12396699) q[2];
sx q[2];
rz(-1.4767246) q[2];
sx q[2];
rz(-1.5383447) q[2];
rz(-2.7247143) q[3];
sx q[3];
rz(-2.6346801) q[3];
sx q[3];
rz(-0.1703593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040319547) q[0];
sx q[0];
rz(-0.23389255) q[0];
sx q[0];
rz(-0.42863578) q[0];
rz(-2.0173232) q[1];
sx q[1];
rz(-1.602403) q[1];
sx q[1];
rz(0.7695778) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3807629) q[0];
sx q[0];
rz(-1.5676982) q[0];
sx q[0];
rz(0.0018381434) q[0];
x q[1];
rz(0.52346241) q[2];
sx q[2];
rz(-2.0085196) q[2];
sx q[2];
rz(1.4283534) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.10248549) q[1];
sx q[1];
rz(-1.277248) q[1];
sx q[1];
rz(2.5204646) q[1];
rz(-1.0060723) q[3];
sx q[3];
rz(-0.58593633) q[3];
sx q[3];
rz(1.0387109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.6577242) q[2];
sx q[2];
rz(-1.8572073) q[2];
sx q[2];
rz(-0.095495187) q[2];
rz(2.5996082) q[3];
sx q[3];
rz(-0.46613765) q[3];
sx q[3];
rz(3.0123762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99590456) q[0];
sx q[0];
rz(-0.93038428) q[0];
sx q[0];
rz(-0.11181871) q[0];
rz(-1.3189141) q[1];
sx q[1];
rz(-1.8967862) q[1];
sx q[1];
rz(-1.8359312) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3506136) q[0];
sx q[0];
rz(-1.0270938) q[0];
sx q[0];
rz(0.14474317) q[0];
rz(2.6839031) q[2];
sx q[2];
rz(-0.31186843) q[2];
sx q[2];
rz(2.092474) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9010734) q[1];
sx q[1];
rz(-1.8522634) q[1];
sx q[1];
rz(2.752384) q[1];
rz(-pi) q[2];
rz(0.54763973) q[3];
sx q[3];
rz(-1.9353011) q[3];
sx q[3];
rz(0.49307666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.32885113) q[2];
sx q[2];
rz(-1.1527656) q[2];
sx q[2];
rz(1.1197155) q[2];
rz(2.8730734) q[3];
sx q[3];
rz(-1.7374141) q[3];
sx q[3];
rz(-1.9560811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8609404) q[0];
sx q[0];
rz(-2.379731) q[0];
sx q[0];
rz(-2.5326488) q[0];
rz(-2.2641585) q[1];
sx q[1];
rz(-2.139822) q[1];
sx q[1];
rz(-2.529349) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5800047) q[0];
sx q[0];
rz(-1.9733493) q[0];
sx q[0];
rz(-0.7640362) q[0];
x q[1];
rz(-1.2127146) q[2];
sx q[2];
rz(-2.1109258) q[2];
sx q[2];
rz(-3.0703406) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.60499268) q[1];
sx q[1];
rz(-0.92291622) q[1];
sx q[1];
rz(0.98439321) q[1];
rz(2.6004535) q[3];
sx q[3];
rz(-1.637646) q[3];
sx q[3];
rz(-2.5802286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0683384) q[2];
sx q[2];
rz(-1.5404258) q[2];
sx q[2];
rz(0.18745984) q[2];
rz(-2.0295664) q[3];
sx q[3];
rz(-0.91809648) q[3];
sx q[3];
rz(2.658127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40792313) q[0];
sx q[0];
rz(-2.3275571) q[0];
sx q[0];
rz(0.2440051) q[0];
rz(-1.4834652) q[1];
sx q[1];
rz(-1.6226059) q[1];
sx q[1];
rz(-0.77267486) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84757728) q[0];
sx q[0];
rz(-0.71888598) q[0];
sx q[0];
rz(-0.80475828) q[0];
rz(-1.7079321) q[2];
sx q[2];
rz(-1.8923645) q[2];
sx q[2];
rz(-1.1346045) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.15068842) q[1];
sx q[1];
rz(-2.7953618) q[1];
sx q[1];
rz(-2.8105286) q[1];
x q[2];
rz(2.3540007) q[3];
sx q[3];
rz(-2.1887472) q[3];
sx q[3];
rz(1.7389115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3385758) q[2];
sx q[2];
rz(-0.96119857) q[2];
sx q[2];
rz(2.4147066) q[2];
rz(1.1072985) q[3];
sx q[3];
rz(-0.50873435) q[3];
sx q[3];
rz(1.0072964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1399287) q[0];
sx q[0];
rz(-1.7353084) q[0];
sx q[0];
rz(-1.1575862) q[0];
rz(-2.6955556) q[1];
sx q[1];
rz(-2.0560494) q[1];
sx q[1];
rz(2.5037419) q[1];
rz(-1.695651) q[2];
sx q[2];
rz(-1.1706252) q[2];
sx q[2];
rz(-1.3782383) q[2];
rz(1.4230396) q[3];
sx q[3];
rz(-1.9922517) q[3];
sx q[3];
rz(-1.7094517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
