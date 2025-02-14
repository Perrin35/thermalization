OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1287542) q[0];
sx q[0];
rz(-1.1829809) q[0];
sx q[0];
rz(0.24721375) q[0];
rz(-1.5850868) q[1];
sx q[1];
rz(4.1559846) q[1];
sx q[1];
rz(13.520887) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5069141) q[0];
sx q[0];
rz(-1.5999113) q[0];
sx q[0];
rz(1.1965189) q[0];
rz(-pi) q[1];
rz(0.45329161) q[2];
sx q[2];
rz(-1.5201836) q[2];
sx q[2];
rz(2.5091189) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40084154) q[1];
sx q[1];
rz(-0.84226006) q[1];
sx q[1];
rz(0.42899112) q[1];
rz(-pi) q[2];
rz(-2.3684816) q[3];
sx q[3];
rz(-1.793981) q[3];
sx q[3];
rz(3.0497568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9556094) q[2];
sx q[2];
rz(-2.535203) q[2];
sx q[2];
rz(0.28140226) q[2];
rz(-1.227281) q[3];
sx q[3];
rz(-1.3325007) q[3];
sx q[3];
rz(1.9277771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1948552) q[0];
sx q[0];
rz(-0.40732107) q[0];
sx q[0];
rz(2.5545004) q[0];
rz(0.52014822) q[1];
sx q[1];
rz(-1.2112434) q[1];
sx q[1];
rz(-0.21336666) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0555455) q[0];
sx q[0];
rz(-0.16250691) q[0];
sx q[0];
rz(-1.9616198) q[0];
rz(0.51220973) q[2];
sx q[2];
rz(-0.98645112) q[2];
sx q[2];
rz(-2.0456631) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50700106) q[1];
sx q[1];
rz(-1.4774972) q[1];
sx q[1];
rz(-0.044528151) q[1];
rz(-pi) q[2];
rz(-2.0961833) q[3];
sx q[3];
rz(-1.8118708) q[3];
sx q[3];
rz(1.8808533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1063891) q[2];
sx q[2];
rz(-0.28855244) q[2];
sx q[2];
rz(-2.8698548) q[2];
rz(2.4954097) q[3];
sx q[3];
rz(-1.313442) q[3];
sx q[3];
rz(1.9338231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81994098) q[0];
sx q[0];
rz(-0.8701179) q[0];
sx q[0];
rz(2.8735549) q[0];
rz(2.9053814) q[1];
sx q[1];
rz(-1.7832489) q[1];
sx q[1];
rz(1.7265629) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38430308) q[0];
sx q[0];
rz(-1.4688604) q[0];
sx q[0];
rz(1.6166286) q[0];
rz(-0.36249749) q[2];
sx q[2];
rz(-1.4115184) q[2];
sx q[2];
rz(2.8822219) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.336795) q[1];
sx q[1];
rz(-2.6141254) q[1];
sx q[1];
rz(2.5673366) q[1];
rz(0.3424267) q[3];
sx q[3];
rz(-0.70243109) q[3];
sx q[3];
rz(-0.46063761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95431027) q[2];
sx q[2];
rz(-2.8128746) q[2];
sx q[2];
rz(-1.7598565) q[2];
rz(0.36661822) q[3];
sx q[3];
rz(-1.6691875) q[3];
sx q[3];
rz(3.0439175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5919375) q[0];
sx q[0];
rz(-2.9357935) q[0];
sx q[0];
rz(0.016121443) q[0];
rz(-1.6751809) q[1];
sx q[1];
rz(-1.5212395) q[1];
sx q[1];
rz(0.88502562) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21986905) q[0];
sx q[0];
rz(-1.9737287) q[0];
sx q[0];
rz(1.7582488) q[0];
rz(-pi) q[1];
rz(2.1145203) q[2];
sx q[2];
rz(-2.7197065) q[2];
sx q[2];
rz(-2.2510901) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.414847) q[1];
sx q[1];
rz(-1.1903084) q[1];
sx q[1];
rz(0.49729113) q[1];
rz(-pi) q[2];
rz(-2.0871986) q[3];
sx q[3];
rz(-1.2149016) q[3];
sx q[3];
rz(3.0404224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.57899388) q[2];
sx q[2];
rz(-2.6851974) q[2];
sx q[2];
rz(-1.558051) q[2];
rz(1.7700178) q[3];
sx q[3];
rz(-1.2746425) q[3];
sx q[3];
rz(2.9461327) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9850995) q[0];
sx q[0];
rz(-1.766196) q[0];
sx q[0];
rz(-2.9920355) q[0];
rz(2.0984446) q[1];
sx q[1];
rz(-0.97277343) q[1];
sx q[1];
rz(-2.3057888) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5468109) q[0];
sx q[0];
rz(-1.321209) q[0];
sx q[0];
rz(-1.9052192) q[0];
x q[1];
rz(-1.582242) q[2];
sx q[2];
rz(-0.69037333) q[2];
sx q[2];
rz(2.9857218) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0479314) q[1];
sx q[1];
rz(-1.7788299) q[1];
sx q[1];
rz(0.00079493513) q[1];
x q[2];
rz(-1.3204783) q[3];
sx q[3];
rz(-1.3930818) q[3];
sx q[3];
rz(1.632393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.122637) q[2];
sx q[2];
rz(-1.476172) q[2];
sx q[2];
rz(-2.7254851) q[2];
rz(3.05919) q[3];
sx q[3];
rz(-2.1722138) q[3];
sx q[3];
rz(-2.940322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0788954) q[0];
sx q[0];
rz(-2.1128928) q[0];
sx q[0];
rz(1.0908352) q[0];
rz(2.003032) q[1];
sx q[1];
rz(-1.9416315) q[1];
sx q[1];
rz(-2.535215) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2045468) q[0];
sx q[0];
rz(-1.318745) q[0];
sx q[0];
rz(-2.4840647) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44754812) q[2];
sx q[2];
rz(-0.85361631) q[2];
sx q[2];
rz(0.56499519) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.55058897) q[1];
sx q[1];
rz(-0.85500756) q[1];
sx q[1];
rz(2.4229933) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71833937) q[3];
sx q[3];
rz(-1.4786309) q[3];
sx q[3];
rz(-0.80764333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.828317) q[2];
sx q[2];
rz(-0.54855359) q[2];
sx q[2];
rz(0.11087785) q[2];
rz(1.5634792) q[3];
sx q[3];
rz(-0.21136798) q[3];
sx q[3];
rz(-1.3255239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9368847) q[0];
sx q[0];
rz(-0.81975833) q[0];
sx q[0];
rz(2.4844266) q[0];
rz(1.8183297) q[1];
sx q[1];
rz(-0.75138775) q[1];
sx q[1];
rz(2.0164067) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.781625) q[0];
sx q[0];
rz(-3.1198745) q[0];
sx q[0];
rz(2.3853073) q[0];
x q[1];
rz(0.31382896) q[2];
sx q[2];
rz(-1.4782566) q[2];
sx q[2];
rz(1.0864663) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4328496) q[1];
sx q[1];
rz(-1.1869935) q[1];
sx q[1];
rz(-0.52977458) q[1];
rz(1.8784889) q[3];
sx q[3];
rz(-2.2224226) q[3];
sx q[3];
rz(-0.055881413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7323759) q[2];
sx q[2];
rz(-1.9595307) q[2];
sx q[2];
rz(2.1803975) q[2];
rz(-3.1327278) q[3];
sx q[3];
rz(-2.2564087) q[3];
sx q[3];
rz(-1.2113021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0601785) q[0];
sx q[0];
rz(-0.1135122) q[0];
sx q[0];
rz(-0.42763448) q[0];
rz(-2.0291406) q[1];
sx q[1];
rz(-1.8266269) q[1];
sx q[1];
rz(-2.1449259) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0286197) q[0];
sx q[0];
rz(-2.1899208) q[0];
sx q[0];
rz(2.1834247) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2995424) q[2];
sx q[2];
rz(-2.315186) q[2];
sx q[2];
rz(-2.8959993) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.61178798) q[1];
sx q[1];
rz(-0.63358297) q[1];
sx q[1];
rz(1.1128913) q[1];
x q[2];
rz(0.57076591) q[3];
sx q[3];
rz(-1.7460572) q[3];
sx q[3];
rz(1.4890763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4619649) q[2];
sx q[2];
rz(-2.6042016) q[2];
sx q[2];
rz(-2.4930387) q[2];
rz(-2.8280761) q[3];
sx q[3];
rz(-0.88839141) q[3];
sx q[3];
rz(3.0891109) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40470966) q[0];
sx q[0];
rz(-0.22060224) q[0];
sx q[0];
rz(2.1746461) q[0];
rz(2.1838358) q[1];
sx q[1];
rz(-2.0499332) q[1];
sx q[1];
rz(-1.0831833) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5388018) q[0];
sx q[0];
rz(-1.9968918) q[0];
sx q[0];
rz(-1.2998796) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71523209) q[2];
sx q[2];
rz(-0.87022793) q[2];
sx q[2];
rz(2.689555) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5965278) q[1];
sx q[1];
rz(-1.178243) q[1];
sx q[1];
rz(0.82253455) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6622436) q[3];
sx q[3];
rz(-1.5707301) q[3];
sx q[3];
rz(0.48653761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0208685) q[2];
sx q[2];
rz(-0.82896295) q[2];
sx q[2];
rz(-1.6797569) q[2];
rz(-0.36137897) q[3];
sx q[3];
rz(-1.8166108) q[3];
sx q[3];
rz(-2.1593275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7892889) q[0];
sx q[0];
rz(-2.4109349) q[0];
sx q[0];
rz(-0.58255449) q[0];
rz(1.7763058) q[1];
sx q[1];
rz(-1.0368232) q[1];
sx q[1];
rz(2.9127311) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.796237) q[0];
sx q[0];
rz(-1.7526049) q[0];
sx q[0];
rz(2.7391091) q[0];
x q[1];
rz(2.1624915) q[2];
sx q[2];
rz(-2.8276463) q[2];
sx q[2];
rz(0.57127956) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4881674) q[1];
sx q[1];
rz(-1.5060802) q[1];
sx q[1];
rz(1.4117227) q[1];
rz(-1.1353605) q[3];
sx q[3];
rz(-1.2139911) q[3];
sx q[3];
rz(0.26998587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2710496) q[2];
sx q[2];
rz(-2.5898263) q[2];
sx q[2];
rz(-1.7331227) q[2];
rz(-2.4208141) q[3];
sx q[3];
rz(-1.9920789) q[3];
sx q[3];
rz(-1.6225947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085717399) q[0];
sx q[0];
rz(-1.4376823) q[0];
sx q[0];
rz(1.8327376) q[0];
rz(1.5695288) q[1];
sx q[1];
rz(-2.5253898) q[1];
sx q[1];
rz(0.22402221) q[1];
rz(-2.1355676) q[2];
sx q[2];
rz(-1.7823162) q[2];
sx q[2];
rz(1.6020365) q[2];
rz(1.4103945) q[3];
sx q[3];
rz(-1.4199724) q[3];
sx q[3];
rz(1.7589105) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
