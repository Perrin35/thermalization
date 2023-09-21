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
rz(-2.8744856) q[1];
sx q[1];
rz(-2.5565956) q[1];
sx q[1];
rz(-2.4490228) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4953295) q[0];
sx q[0];
rz(-1.6705992) q[0];
sx q[0];
rz(-2.0908337) q[0];
rz(-0.59074596) q[2];
sx q[2];
rz(-2.4060537) q[2];
sx q[2];
rz(2.9172446) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.270826) q[1];
sx q[1];
rz(-1.4987136) q[1];
sx q[1];
rz(1.6986398) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0308706) q[3];
sx q[3];
rz(-0.51957909) q[3];
sx q[3];
rz(-0.67668623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1074368) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(1.5365323) q[2];
rz(1.6202554) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(3.1055514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.3261616) q[0];
sx q[0];
rz(-2.8679929) q[0];
sx q[0];
rz(-1.8923627) q[0];
rz(0.56150395) q[1];
sx q[1];
rz(-0.77604547) q[1];
sx q[1];
rz(-2.5610279) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55727977) q[0];
sx q[0];
rz(-1.4627539) q[0];
sx q[0];
rz(2.8312107) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63802982) q[2];
sx q[2];
rz(-1.3943765) q[2];
sx q[2];
rz(2.1413013) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7932574) q[1];
sx q[1];
rz(-1.0151334) q[1];
sx q[1];
rz(0.72850119) q[1];
rz(-2.0248236) q[3];
sx q[3];
rz(-1.2516216) q[3];
sx q[3];
rz(-1.9516731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4124174) q[2];
sx q[2];
rz(-2.9340332) q[2];
sx q[2];
rz(2.2632329) q[2];
rz(-2.7495524) q[3];
sx q[3];
rz(-1.4441898) q[3];
sx q[3];
rz(0.6033321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6648401) q[0];
sx q[0];
rz(-2.219253) q[0];
sx q[0];
rz(-0.77366775) q[0];
rz(3.1402918) q[1];
sx q[1];
rz(-1.6157849) q[1];
sx q[1];
rz(-3.1087648) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5531909) q[0];
sx q[0];
rz(-1.2598039) q[0];
sx q[0];
rz(1.3283967) q[0];
rz(-pi) q[1];
rz(-0.29067729) q[2];
sx q[2];
rz(-1.0376087) q[2];
sx q[2];
rz(2.3454587) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58363885) q[1];
sx q[1];
rz(-1.7350405) q[1];
sx q[1];
rz(-2.3418531) q[1];
rz(-pi) q[2];
rz(0.34605108) q[3];
sx q[3];
rz(-1.0533353) q[3];
sx q[3];
rz(-2.3164761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7971928) q[2];
sx q[2];
rz(-1.025528) q[2];
sx q[2];
rz(0.27734217) q[2];
rz(2.7456361) q[3];
sx q[3];
rz(-1.6010511) q[3];
sx q[3];
rz(-0.69916454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
rz(-0.83507744) q[1];
sx q[1];
rz(2.3707726) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9642826) q[0];
sx q[0];
rz(-1.5872123) q[0];
sx q[0];
rz(-1.5903227) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0609444) q[2];
sx q[2];
rz(-1.66301) q[2];
sx q[2];
rz(-0.25445081) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.091688823) q[1];
sx q[1];
rz(-2.4871475) q[1];
sx q[1];
rz(-2.5212538) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9510668) q[3];
sx q[3];
rz(-1.7912205) q[3];
sx q[3];
rz(-0.64173736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.16584855) q[2];
sx q[2];
rz(-1.5074915) q[2];
sx q[2];
rz(-0.7129933) q[2];
rz(2.1285848) q[3];
sx q[3];
rz(-0.37390798) q[3];
sx q[3];
rz(-2.1876984) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35904303) q[0];
sx q[0];
rz(-1.0721711) q[0];
sx q[0];
rz(1.7011401) q[0];
rz(-0.094093181) q[1];
sx q[1];
rz(-0.73939878) q[1];
sx q[1];
rz(2.9715911) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48430303) q[0];
sx q[0];
rz(-2.347725) q[0];
sx q[0];
rz(-2.8110709) q[0];
rz(2.0427809) q[2];
sx q[2];
rz(-1.7760279) q[2];
sx q[2];
rz(2.1489378) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4449094) q[1];
sx q[1];
rz(-1.0001567) q[1];
sx q[1];
rz(0.48008227) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4835028) q[3];
sx q[3];
rz(-2.1562025) q[3];
sx q[3];
rz(-0.64774367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.99469441) q[2];
sx q[2];
rz(-2.6165104) q[2];
sx q[2];
rz(-1.4040995) q[2];
rz(1.5480301) q[3];
sx q[3];
rz(-2.352495) q[3];
sx q[3];
rz(-1.7061957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.4858522) q[0];
sx q[0];
rz(-1.9952554) q[0];
sx q[0];
rz(-0.45853841) q[0];
rz(2.8857152) q[1];
sx q[1];
rz(-1.2586539) q[1];
sx q[1];
rz(0.68516723) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79294369) q[0];
sx q[0];
rz(-0.92175882) q[0];
sx q[0];
rz(3.0543324) q[0];
rz(0.77722705) q[2];
sx q[2];
rz(-1.6738552) q[2];
sx q[2];
rz(-1.7412141) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31720686) q[1];
sx q[1];
rz(-2.3134391) q[1];
sx q[1];
rz(0.87888996) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1567247) q[3];
sx q[3];
rz(-1.5555447) q[3];
sx q[3];
rz(-1.0451942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7489862) q[2];
sx q[2];
rz(-2.7219153) q[2];
sx q[2];
rz(1.7123327) q[2];
rz(1.0990934) q[3];
sx q[3];
rz(-0.50656879) q[3];
sx q[3];
rz(-2.9523622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-3.1290865) q[0];
sx q[0];
rz(-1.5959473) q[0];
sx q[0];
rz(-2.4122453) q[0];
rz(0.29306456) q[1];
sx q[1];
rz(-0.23935071) q[1];
sx q[1];
rz(1.1475295) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1345394) q[0];
sx q[0];
rz(-0.42047406) q[0];
sx q[0];
rz(0.18376952) q[0];
rz(-pi) q[1];
rz(2.7717934) q[2];
sx q[2];
rz(-2.1714006) q[2];
sx q[2];
rz(-0.72593867) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.16312576) q[1];
sx q[1];
rz(-1.5297869) q[1];
sx q[1];
rz(-2.6569215) q[1];
rz(-pi) q[2];
rz(1.958193) q[3];
sx q[3];
rz(-2.0347188) q[3];
sx q[3];
rz(1.0138318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.78836936) q[2];
sx q[2];
rz(-1.0228144) q[2];
sx q[2];
rz(0.74907556) q[2];
rz(-0.64368147) q[3];
sx q[3];
rz(-2.1285074) q[3];
sx q[3];
rz(-0.914004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1483243) q[0];
sx q[0];
rz(-2.0355621) q[0];
sx q[0];
rz(-0.83129445) q[0];
rz(1.7656901) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(2.7430699) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9996223) q[0];
sx q[0];
rz(-2.6451783) q[0];
sx q[0];
rz(-1.9103785) q[0];
rz(0.60110809) q[2];
sx q[2];
rz(-1.6582487) q[2];
sx q[2];
rz(-2.2488038) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1426705) q[1];
sx q[1];
rz(-2.167836) q[1];
sx q[1];
rz(3.1139657) q[1];
rz(-pi) q[2];
rz(2.1274444) q[3];
sx q[3];
rz(-1.1245407) q[3];
sx q[3];
rz(2.2198912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9514256) q[2];
sx q[2];
rz(-1.7931033) q[2];
sx q[2];
rz(1.8438967) q[2];
rz(1.1249582) q[3];
sx q[3];
rz(-1.2599895) q[3];
sx q[3];
rz(0.64363939) q[3];
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
rz(0.012638906) q[0];
sx q[0];
rz(-2.5057827) q[0];
sx q[0];
rz(-1.3893611) q[0];
rz(1.6268436) q[1];
sx q[1];
rz(-1.4667958) q[1];
sx q[1];
rz(2.0432037) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5441355) q[0];
sx q[0];
rz(-1.3381334) q[0];
sx q[0];
rz(1.0730037) q[0];
x q[1];
rz(0.45778747) q[2];
sx q[2];
rz(-0.87399235) q[2];
sx q[2];
rz(0.84504499) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.64975148) q[1];
sx q[1];
rz(-1.5702618) q[1];
sx q[1];
rz(-0.57926308) q[1];
rz(-pi) q[2];
rz(1.3994201) q[3];
sx q[3];
rz(-1.3203353) q[3];
sx q[3];
rz(-2.7385538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2417458) q[2];
sx q[2];
rz(-1.8490303) q[2];
sx q[2];
rz(0.51952726) q[2];
rz(1.8064921) q[3];
sx q[3];
rz(-2.3050008) q[3];
sx q[3];
rz(-1.8241204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39524233) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(-0.46646068) q[0];
rz(0.17164224) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(-0.62896532) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2254612) q[0];
sx q[0];
rz(-0.3779419) q[0];
sx q[0];
rz(1.7037237) q[0];
rz(-pi) q[1];
rz(-1.2538818) q[2];
sx q[2];
rz(-2.4744518) q[2];
sx q[2];
rz(0.72310477) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.65044636) q[1];
sx q[1];
rz(-1.568734) q[1];
sx q[1];
rz(-2.6994929) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4217989) q[3];
sx q[3];
rz(-1.3070953) q[3];
sx q[3];
rz(-1.7566453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.24370596) q[2];
sx q[2];
rz(-1.0906929) q[2];
sx q[2];
rz(2.0824599) q[2];
rz(0.6774261) q[3];
sx q[3];
rz(-2.1493561) q[3];
sx q[3];
rz(-0.89390755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28329904) q[0];
sx q[0];
rz(-1.1245921) q[0];
sx q[0];
rz(-0.84381214) q[0];
rz(2.8876866) q[1];
sx q[1];
rz(-2.0575247) q[1];
sx q[1];
rz(2.5622096) q[1];
rz(-1.614744) q[2];
sx q[2];
rz(-0.84779253) q[2];
sx q[2];
rz(0.070889125) q[2];
rz(3.0729978) q[3];
sx q[3];
rz(-2.1566236) q[3];
sx q[3];
rz(-1.9038283) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
