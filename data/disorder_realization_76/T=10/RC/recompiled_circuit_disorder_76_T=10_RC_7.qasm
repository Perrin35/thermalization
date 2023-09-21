OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27622142) q[0];
sx q[0];
rz(5.4260317) q[0];
sx q[0];
rz(9.5572588) q[0];
rz(0.26710701) q[1];
sx q[1];
rz(5.6981882) q[1];
sx q[1];
rz(11.873801) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4953295) q[0];
sx q[0];
rz(-1.4709934) q[0];
sx q[0];
rz(-2.0908337) q[0];
rz(1.1039256) q[2];
sx q[2];
rz(-0.97969998) q[2];
sx q[2];
rz(2.1819654) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8323648) q[1];
sx q[1];
rz(-1.698306) q[1];
sx q[1];
rz(-0.072673701) q[1];
rz(-pi) q[2];
rz(-2.0444144) q[3];
sx q[3];
rz(-1.3485104) q[3];
sx q[3];
rz(-1.8412631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0341558) q[2];
sx q[2];
rz(-0.52705708) q[2];
sx q[2];
rz(-1.6050603) q[2];
rz(-1.5213373) q[3];
sx q[3];
rz(-1.6531569) q[3];
sx q[3];
rz(-3.1055514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3261616) q[0];
sx q[0];
rz(-2.8679929) q[0];
sx q[0];
rz(-1.2492299) q[0];
rz(-2.5800887) q[1];
sx q[1];
rz(-0.77604547) q[1];
sx q[1];
rz(0.5805648) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1626494) q[0];
sx q[0];
rz(-1.8793086) q[0];
sx q[0];
rz(-1.4573775) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63802982) q[2];
sx q[2];
rz(-1.7472162) q[2];
sx q[2];
rz(-1.0002913) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21743449) q[1];
sx q[1];
rz(-2.1719451) q[1];
sx q[1];
rz(2.2648328) q[1];
rz(-pi) q[2];
rz(2.0248236) q[3];
sx q[3];
rz(-1.8899711) q[3];
sx q[3];
rz(-1.9516731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4124174) q[2];
sx q[2];
rz(-2.9340332) q[2];
sx q[2];
rz(2.2632329) q[2];
rz(0.39204028) q[3];
sx q[3];
rz(-1.4441898) q[3];
sx q[3];
rz(-2.5382606) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6648401) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(-2.3679249) q[0];
rz(-0.0013008612) q[1];
sx q[1];
rz(-1.6157849) q[1];
sx q[1];
rz(0.032827854) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5531909) q[0];
sx q[0];
rz(-1.2598039) q[0];
sx q[0];
rz(-1.3283967) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1187395) q[2];
sx q[2];
rz(-2.5411122) q[2];
sx q[2];
rz(-0.26417363) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.58363885) q[1];
sx q[1];
rz(-1.4065521) q[1];
sx q[1];
rz(0.79973952) q[1];
x q[2];
rz(-2.7955416) q[3];
sx q[3];
rz(-2.0882574) q[3];
sx q[3];
rz(-0.82511653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.34439987) q[2];
sx q[2];
rz(-2.1160647) q[2];
sx q[2];
rz(-0.27734217) q[2];
rz(0.39595655) q[3];
sx q[3];
rz(-1.5405416) q[3];
sx q[3];
rz(2.4424281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4375777) q[0];
sx q[0];
rz(-0.33575785) q[0];
sx q[0];
rz(0.26279703) q[0];
rz(-2.908005) q[1];
sx q[1];
rz(-2.3065152) q[1];
sx q[1];
rz(0.77082005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9642826) q[0];
sx q[0];
rz(-1.5543803) q[0];
sx q[0];
rz(1.55127) q[0];
rz(-pi) q[1];
rz(0.1044354) q[2];
sx q[2];
rz(-1.0829139) q[2];
sx q[2];
rz(-1.7761531) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.82476393) q[1];
sx q[1];
rz(-2.0889805) q[1];
sx q[1];
rz(-1.9903241) q[1];
x q[2];
rz(-0.1905258) q[3];
sx q[3];
rz(-1.7912205) q[3];
sx q[3];
rz(2.4998553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.16584855) q[2];
sx q[2];
rz(-1.5074915) q[2];
sx q[2];
rz(2.4285994) q[2];
rz(2.1285848) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(-0.95389429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35904303) q[0];
sx q[0];
rz(-2.0694216) q[0];
sx q[0];
rz(-1.7011401) q[0];
rz(0.094093181) q[1];
sx q[1];
rz(-0.73939878) q[1];
sx q[1];
rz(-2.9715911) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3225587) q[0];
sx q[0];
rz(-1.8043307) q[0];
sx q[0];
rz(0.76604953) q[0];
rz(-pi) q[1];
rz(-1.0988118) q[2];
sx q[2];
rz(-1.7760279) q[2];
sx q[2];
rz(2.1489378) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1483037) q[1];
sx q[1];
rz(-1.171604) q[1];
sx q[1];
rz(2.1972375) q[1];
rz(2.2682297) q[3];
sx q[3];
rz(-1.0358827) q[3];
sx q[3];
rz(1.3267335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99469441) q[2];
sx q[2];
rz(-2.6165104) q[2];
sx q[2];
rz(1.4040995) q[2];
rz(-1.5935625) q[3];
sx q[3];
rz(-0.78909767) q[3];
sx q[3];
rz(1.7061957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65574044) q[0];
sx q[0];
rz(-1.1463373) q[0];
sx q[0];
rz(0.45853841) q[0];
rz(-0.25587747) q[1];
sx q[1];
rz(-1.8829388) q[1];
sx q[1];
rz(2.4564254) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79294369) q[0];
sx q[0];
rz(-0.92175882) q[0];
sx q[0];
rz(-0.087260212) q[0];
rz(-pi) q[1];
rz(-2.9951727) q[2];
sx q[2];
rz(-0.78260566) q[2];
sx q[2];
rz(0.27461068) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.398715) q[1];
sx q[1];
rz(-2.0600973) q[1];
sx q[1];
rz(0.87280416) q[1];
rz(-1.54322) q[3];
sx q[3];
rz(-0.58610361) q[3];
sx q[3];
rz(-2.5930149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.3926065) q[2];
sx q[2];
rz(-0.41967732) q[2];
sx q[2];
rz(-1.7123327) q[2];
rz(-2.0424992) q[3];
sx q[3];
rz(-0.50656879) q[3];
sx q[3];
rz(-2.9523622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-3.1290865) q[0];
sx q[0];
rz(-1.5959473) q[0];
sx q[0];
rz(-2.4122453) q[0];
rz(-0.29306456) q[1];
sx q[1];
rz(-2.9022419) q[1];
sx q[1];
rz(-1.9940631) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1345394) q[0];
sx q[0];
rz(-0.42047406) q[0];
sx q[0];
rz(0.18376952) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2044264) q[2];
sx q[2];
rz(-1.8735777) q[2];
sx q[2];
rz(-2.5123951) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3299678) q[1];
sx q[1];
rz(-0.48626541) q[1];
sx q[1];
rz(-0.087841308) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64671867) q[3];
sx q[3];
rz(-0.59520703) q[3];
sx q[3];
rz(-1.3884384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.78836936) q[2];
sx q[2];
rz(-2.1187783) q[2];
sx q[2];
rz(0.74907556) q[2];
rz(2.4979112) q[3];
sx q[3];
rz(-2.1285074) q[3];
sx q[3];
rz(2.2275887) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1483243) q[0];
sx q[0];
rz(-1.1060306) q[0];
sx q[0];
rz(2.3102982) q[0];
rz(-1.3759026) q[1];
sx q[1];
rz(-2.3283236) q[1];
sx q[1];
rz(-2.7430699) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4115899) q[0];
sx q[0];
rz(-1.411479) q[0];
sx q[0];
rz(1.0986064) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5404846) q[2];
sx q[2];
rz(-1.483344) q[2];
sx q[2];
rz(-0.89278883) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7290013) q[1];
sx q[1];
rz(-1.593643) q[1];
sx q[1];
rz(0.97357915) q[1];
rz(-pi) q[2];
rz(-1.0141482) q[3];
sx q[3];
rz(-1.1245407) q[3];
sx q[3];
rz(-0.92170148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9514256) q[2];
sx q[2];
rz(-1.3484893) q[2];
sx q[2];
rz(1.8438967) q[2];
rz(-2.0166345) q[3];
sx q[3];
rz(-1.2599895) q[3];
sx q[3];
rz(-2.4979533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1289537) q[0];
sx q[0];
rz(-2.5057827) q[0];
sx q[0];
rz(-1.7522316) q[0];
rz(-1.6268436) q[1];
sx q[1];
rz(-1.4667958) q[1];
sx q[1];
rz(-2.0432037) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2929045) q[0];
sx q[0];
rz(-2.0540037) q[0];
sx q[0];
rz(-0.26341652) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82010014) q[2];
sx q[2];
rz(-1.2250049) q[2];
sx q[2];
rz(-1.0320013) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2201982) q[1];
sx q[1];
rz(-0.99153334) q[1];
sx q[1];
rz(1.5701576) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25402756) q[3];
sx q[3];
rz(-1.7367749) q[3];
sx q[3];
rz(-2.0167054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2417458) q[2];
sx q[2];
rz(-1.8490303) q[2];
sx q[2];
rz(0.51952726) q[2];
rz(1.3351006) q[3];
sx q[3];
rz(-0.83659187) q[3];
sx q[3];
rz(1.3174723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7463503) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(2.675132) q[0];
rz(0.17164224) q[1];
sx q[1];
rz(-1.9263093) q[1];
sx q[1];
rz(0.62896532) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0590203) q[0];
sx q[0];
rz(-1.1963545) q[0];
sx q[0];
rz(-3.0890205) q[0];
rz(-pi) q[1];
rz(-1.2538818) q[2];
sx q[2];
rz(-2.4744518) q[2];
sx q[2];
rz(0.72310477) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.65044636) q[1];
sx q[1];
rz(-1.568734) q[1];
sx q[1];
rz(-0.44209977) q[1];
x q[2];
rz(-1.915515) q[3];
sx q[3];
rz(-0.88092062) q[3];
sx q[3];
rz(2.7310839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8978867) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(1.0591327) q[2];
rz(0.6774261) q[3];
sx q[3];
rz(-2.1493561) q[3];
sx q[3];
rz(-0.89390755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28329904) q[0];
sx q[0];
rz(-2.0170006) q[0];
sx q[0];
rz(2.2977805) q[0];
rz(2.8876866) q[1];
sx q[1];
rz(-2.0575247) q[1];
sx q[1];
rz(2.5622096) q[1];
rz(3.0918448) q[2];
sx q[2];
rz(-2.4174953) q[2];
sx q[2];
rz(0.13725666) q[2];
rz(-1.4678636) q[3];
sx q[3];
rz(-2.5522305) q[3];
sx q[3];
rz(1.1141368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
