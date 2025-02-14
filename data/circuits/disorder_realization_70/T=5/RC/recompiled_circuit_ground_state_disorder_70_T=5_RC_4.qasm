OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.31324759) q[0];
sx q[0];
rz(-2.3656443) q[0];
sx q[0];
rz(1.1021855) q[0];
rz(1.6232396) q[1];
sx q[1];
rz(2.0145388) q[1];
sx q[1];
rz(9.1993499) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091087393) q[0];
sx q[0];
rz(-0.69982547) q[0];
sx q[0];
rz(0.096629337) q[0];
x q[1];
rz(2.5257843) q[2];
sx q[2];
rz(-1.6914509) q[2];
sx q[2];
rz(0.10663154) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0590093) q[1];
sx q[1];
rz(-1.9349001) q[1];
sx q[1];
rz(2.8541982) q[1];
x q[2];
rz(-1.2223689) q[3];
sx q[3];
rz(-2.0392079) q[3];
sx q[3];
rz(0.40701501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6585091) q[2];
sx q[2];
rz(-2.4675214) q[2];
sx q[2];
rz(0.020641208) q[2];
rz(-0.189273) q[3];
sx q[3];
rz(-2.7920189) q[3];
sx q[3];
rz(1.4741723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3610158) q[0];
sx q[0];
rz(-2.3605232) q[0];
sx q[0];
rz(0.97907698) q[0];
rz(0.11115221) q[1];
sx q[1];
rz(-0.73768288) q[1];
sx q[1];
rz(2.0806064) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7951668) q[0];
sx q[0];
rz(-2.0316664) q[0];
sx q[0];
rz(2.5823442) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0060180863) q[2];
sx q[2];
rz(-3.0041347) q[2];
sx q[2];
rz(2.6743367) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0909522) q[1];
sx q[1];
rz(-2.0555858) q[1];
sx q[1];
rz(1.8043955) q[1];
x q[2];
rz(-2.7525192) q[3];
sx q[3];
rz(-1.437646) q[3];
sx q[3];
rz(-1.5885858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4336808) q[2];
sx q[2];
rz(-1.8209063) q[2];
sx q[2];
rz(2.743538) q[2];
rz(1.7178644) q[3];
sx q[3];
rz(-0.28702304) q[3];
sx q[3];
rz(0.48582336) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5240391) q[0];
sx q[0];
rz(-2.6300639) q[0];
sx q[0];
rz(2.0892573) q[0];
rz(2.6984093) q[1];
sx q[1];
rz(-0.96142238) q[1];
sx q[1];
rz(0.65394941) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67414647) q[0];
sx q[0];
rz(-0.83195126) q[0];
sx q[0];
rz(-1.3988926) q[0];
rz(3.0691745) q[2];
sx q[2];
rz(-0.94885599) q[2];
sx q[2];
rz(0.79889132) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.941085) q[1];
sx q[1];
rz(-1.9882747) q[1];
sx q[1];
rz(-2.3509376) q[1];
x q[2];
rz(0.65695936) q[3];
sx q[3];
rz(-2.8553319) q[3];
sx q[3];
rz(3.1050499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7707278) q[2];
sx q[2];
rz(-2.3720522) q[2];
sx q[2];
rz(0.23180836) q[2];
rz(1.2109463) q[3];
sx q[3];
rz(-1.959266) q[3];
sx q[3];
rz(0.32538357) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66832191) q[0];
sx q[0];
rz(-0.40437651) q[0];
sx q[0];
rz(-0.37687287) q[0];
rz(-2.6301774) q[1];
sx q[1];
rz(-1.8835953) q[1];
sx q[1];
rz(-0.84397856) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.224124) q[0];
sx q[0];
rz(-0.11118764) q[0];
sx q[0];
rz(1.7286517) q[0];
rz(-pi) q[1];
x q[1];
rz(2.149821) q[2];
sx q[2];
rz(-1.7536153) q[2];
sx q[2];
rz(2.7046375) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.36605814) q[1];
sx q[1];
rz(-2.0506128) q[1];
sx q[1];
rz(2.5262031) q[1];
x q[2];
rz(-0.52167251) q[3];
sx q[3];
rz(-1.552201) q[3];
sx q[3];
rz(0.90584318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9722998) q[2];
sx q[2];
rz(-2.1268714) q[2];
sx q[2];
rz(0.69804066) q[2];
rz(-1.1997148) q[3];
sx q[3];
rz(-2.2359087) q[3];
sx q[3];
rz(-0.5996632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9508764) q[0];
sx q[0];
rz(-0.27442014) q[0];
sx q[0];
rz(-0.12572591) q[0];
rz(-1.0267286) q[1];
sx q[1];
rz(-1.3871565) q[1];
sx q[1];
rz(-0.52621192) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77718098) q[0];
sx q[0];
rz(-2.6573349) q[0];
sx q[0];
rz(0.15987349) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16321142) q[2];
sx q[2];
rz(-2.2610352) q[2];
sx q[2];
rz(-2.9371967) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6263585) q[1];
sx q[1];
rz(-2.7470416) q[1];
sx q[1];
rz(-0.26923967) q[1];
x q[2];
rz(0.86889699) q[3];
sx q[3];
rz(-2.2326937) q[3];
sx q[3];
rz(1.4922265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0937664) q[2];
sx q[2];
rz(-0.83790773) q[2];
sx q[2];
rz(-0.31025904) q[2];
rz(0.034916498) q[3];
sx q[3];
rz(-1.3860476) q[3];
sx q[3];
rz(-2.3518899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0763615) q[0];
sx q[0];
rz(-1.4820453) q[0];
sx q[0];
rz(-0.37102997) q[0];
rz(3.0767483) q[1];
sx q[1];
rz(-2.3193181) q[1];
sx q[1];
rz(1.2670021) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3716301) q[0];
sx q[0];
rz(-2.2055059) q[0];
sx q[0];
rz(0.43031613) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7899988) q[2];
sx q[2];
rz(-2.1702655) q[2];
sx q[2];
rz(-2.4945187) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.95270009) q[1];
sx q[1];
rz(-1.6052393) q[1];
sx q[1];
rz(-0.39826213) q[1];
rz(-pi) q[2];
rz(-1.2812595) q[3];
sx q[3];
rz(-1.1216333) q[3];
sx q[3];
rz(-1.7790573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0176257) q[2];
sx q[2];
rz(-1.664868) q[2];
sx q[2];
rz(-1.5383447) q[2];
rz(2.7247143) q[3];
sx q[3];
rz(-2.6346801) q[3];
sx q[3];
rz(-2.9712334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040319547) q[0];
sx q[0];
rz(-0.23389255) q[0];
sx q[0];
rz(2.7129569) q[0];
rz(-2.0173232) q[1];
sx q[1];
rz(-1.5391896) q[1];
sx q[1];
rz(-0.7695778) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80996087) q[0];
sx q[0];
rz(-1.5689582) q[0];
sx q[0];
rz(1.5738945) q[0];
x q[1];
rz(-2.0662082) q[2];
sx q[2];
rz(-2.0406463) q[2];
sx q[2];
rz(0.097493492) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0575233) q[1];
sx q[1];
rz(-0.67862679) q[1];
sx q[1];
rz(-2.6625257) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0598211) q[3];
sx q[3];
rz(-1.2703514) q[3];
sx q[3];
rz(2.1238126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.6577242) q[2];
sx q[2];
rz(-1.8572073) q[2];
sx q[2];
rz(0.095495187) q[2];
rz(-0.5419845) q[3];
sx q[3];
rz(-2.675455) q[3];
sx q[3];
rz(-3.0123762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1456881) q[0];
sx q[0];
rz(-2.2112084) q[0];
sx q[0];
rz(-0.11181871) q[0];
rz(1.8226786) q[1];
sx q[1];
rz(-1.8967862) q[1];
sx q[1];
rz(-1.8359312) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14492098) q[0];
sx q[0];
rz(-1.4470417) q[0];
sx q[0];
rz(-2.1191639) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4576896) q[2];
sx q[2];
rz(-0.31186843) q[2];
sx q[2];
rz(-1.0491187) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9255786) q[1];
sx q[1];
rz(-2.6655201) q[1];
sx q[1];
rz(-2.4904597) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9910802) q[3];
sx q[3];
rz(-2.0788361) q[3];
sx q[3];
rz(-2.2779494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8127415) q[2];
sx q[2];
rz(-1.1527656) q[2];
sx q[2];
rz(-2.0218772) q[2];
rz(0.26851922) q[3];
sx q[3];
rz(-1.7374141) q[3];
sx q[3];
rz(-1.1855116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28065228) q[0];
sx q[0];
rz(-2.379731) q[0];
sx q[0];
rz(-2.5326488) q[0];
rz(-0.87743419) q[1];
sx q[1];
rz(-1.0017706) q[1];
sx q[1];
rz(-2.529349) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36831066) q[0];
sx q[0];
rz(-0.88079494) q[0];
sx q[0];
rz(-2.1036214) q[0];
rz(-pi) q[1];
rz(-2.5721278) q[2];
sx q[2];
rz(-1.8761022) q[2];
sx q[2];
rz(-1.832163) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3471191) q[1];
sx q[1];
rz(-1.11382) q[1];
sx q[1];
rz(0.73758482) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12924592) q[3];
sx q[3];
rz(-0.5448444) q[3];
sx q[3];
rz(0.89873492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0683384) q[2];
sx q[2];
rz(-1.5404258) q[2];
sx q[2];
rz(-0.18745984) q[2];
rz(2.0295664) q[3];
sx q[3];
rz(-0.91809648) q[3];
sx q[3];
rz(-2.658127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7336695) q[0];
sx q[0];
rz(-2.3275571) q[0];
sx q[0];
rz(-2.8975876) q[0];
rz(-1.6581274) q[1];
sx q[1];
rz(-1.5189867) q[1];
sx q[1];
rz(2.3689178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84757728) q[0];
sx q[0];
rz(-2.4227067) q[0];
sx q[0];
rz(-2.3368344) q[0];
rz(-pi) q[1];
rz(1.4336606) q[2];
sx q[2];
rz(-1.8923645) q[2];
sx q[2];
rz(-1.1346045) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.15068842) q[1];
sx q[1];
rz(-2.7953618) q[1];
sx q[1];
rz(0.33106403) q[1];
rz(2.3540007) q[3];
sx q[3];
rz(-0.95284546) q[3];
sx q[3];
rz(1.4026812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3385758) q[2];
sx q[2];
rz(-0.96119857) q[2];
sx q[2];
rz(-0.72688603) q[2];
rz(-1.1072985) q[3];
sx q[3];
rz(-0.50873435) q[3];
sx q[3];
rz(-1.0072964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0016639391) q[0];
sx q[0];
rz(-1.7353084) q[0];
sx q[0];
rz(-1.1575862) q[0];
rz(-2.6955556) q[1];
sx q[1];
rz(-2.0560494) q[1];
sx q[1];
rz(2.5037419) q[1];
rz(2.8552796) q[2];
sx q[2];
rz(-2.7234017) q[2];
sx q[2];
rz(-1.0665733) q[2];
rz(1.7185531) q[3];
sx q[3];
rz(-1.1493409) q[3];
sx q[3];
rz(1.4321409) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
