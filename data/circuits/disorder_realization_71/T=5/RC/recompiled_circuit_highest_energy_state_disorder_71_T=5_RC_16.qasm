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
rz(-2.0105536) q[0];
sx q[0];
rz(3.717489) q[0];
sx q[0];
rz(5.1562638) q[0];
rz(-0.95070401) q[1];
sx q[1];
rz(-0.60125142) q[1];
sx q[1];
rz(2.8079005) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93984825) q[0];
sx q[0];
rz(-2.6932095) q[0];
sx q[0];
rz(-1.9223619) q[0];
rz(-pi) q[1];
rz(0.26716455) q[2];
sx q[2];
rz(-1.1295302) q[2];
sx q[2];
rz(1.6187606) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.04449439) q[1];
sx q[1];
rz(-1.3584029) q[1];
sx q[1];
rz(-0.35081671) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2980311) q[3];
sx q[3];
rz(-0.78744167) q[3];
sx q[3];
rz(1.3389325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2600962) q[2];
sx q[2];
rz(-2.0659955) q[2];
sx q[2];
rz(-1.0502226) q[2];
rz(2.8460734) q[3];
sx q[3];
rz(-2.3727356) q[3];
sx q[3];
rz(1.7723005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9722026) q[0];
sx q[0];
rz(-2.1222332) q[0];
sx q[0];
rz(-2.4743359) q[0];
rz(-0.15070209) q[1];
sx q[1];
rz(-1.0548016) q[1];
sx q[1];
rz(-2.8299455) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.920645) q[0];
sx q[0];
rz(-1.049982) q[0];
sx q[0];
rz(-3.060311) q[0];
x q[1];
rz(-2.0097334) q[2];
sx q[2];
rz(-1.8510185) q[2];
sx q[2];
rz(-2.9228022) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2978486) q[1];
sx q[1];
rz(-2.9103855) q[1];
sx q[1];
rz(2.48965) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4979579) q[3];
sx q[3];
rz(-1.967481) q[3];
sx q[3];
rz(1.14183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6665337) q[2];
sx q[2];
rz(-0.89343137) q[2];
sx q[2];
rz(2.3739572) q[2];
rz(0.4808937) q[3];
sx q[3];
rz(-0.77156639) q[3];
sx q[3];
rz(-1.5546999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84592205) q[0];
sx q[0];
rz(-1.5376115) q[0];
sx q[0];
rz(2.3620102) q[0];
rz(-2.7698611) q[1];
sx q[1];
rz(-2.1340243) q[1];
sx q[1];
rz(-2.380611) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8801988) q[0];
sx q[0];
rz(-0.43873271) q[0];
sx q[0];
rz(-2.3068271) q[0];
rz(-pi) q[1];
rz(0.8849319) q[2];
sx q[2];
rz(-1.2239309) q[2];
sx q[2];
rz(0.82277966) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5058712) q[1];
sx q[1];
rz(-0.85800115) q[1];
sx q[1];
rz(0.21943032) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1366421) q[3];
sx q[3];
rz(-2.288886) q[3];
sx q[3];
rz(-0.45482054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.80428213) q[2];
sx q[2];
rz(-0.90039841) q[2];
sx q[2];
rz(0.4551746) q[2];
rz(-2.4496487) q[3];
sx q[3];
rz(-2.4271836) q[3];
sx q[3];
rz(0.80879319) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1728425) q[0];
sx q[0];
rz(-1.3469561) q[0];
sx q[0];
rz(0.2562879) q[0];
rz(1.434727) q[1];
sx q[1];
rz(-1.4204357) q[1];
sx q[1];
rz(0.11875471) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8855806) q[0];
sx q[0];
rz(-0.47614723) q[0];
sx q[0];
rz(0.34389596) q[0];
rz(-pi) q[1];
rz(-3.0135113) q[2];
sx q[2];
rz(-1.2799147) q[2];
sx q[2];
rz(-2.0258935) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.465047) q[1];
sx q[1];
rz(-0.80261723) q[1];
sx q[1];
rz(1.2458744) q[1];
rz(1.1944653) q[3];
sx q[3];
rz(-2.2897165) q[3];
sx q[3];
rz(0.014499078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8541096) q[2];
sx q[2];
rz(-2.867925) q[2];
sx q[2];
rz(1.0556833) q[2];
rz(1.4500729) q[3];
sx q[3];
rz(-1.6943211) q[3];
sx q[3];
rz(2.292574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5797822) q[0];
sx q[0];
rz(-2.9636443) q[0];
sx q[0];
rz(-1.5280888) q[0];
rz(-1.1771419) q[1];
sx q[1];
rz(-1.8714995) q[1];
sx q[1];
rz(-1.5783763) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5185243) q[0];
sx q[0];
rz(-0.94618624) q[0];
sx q[0];
rz(-1.6562069) q[0];
rz(-pi) q[1];
rz(-1.8991532) q[2];
sx q[2];
rz(-1.9617726) q[2];
sx q[2];
rz(-2.7308488) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8635378) q[1];
sx q[1];
rz(-1.7124933) q[1];
sx q[1];
rz(-0.63512953) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1354229) q[3];
sx q[3];
rz(-0.4094165) q[3];
sx q[3];
rz(0.96640429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13581181) q[2];
sx q[2];
rz(-1.9848738) q[2];
sx q[2];
rz(1.8394252) q[2];
rz(2.8016413) q[3];
sx q[3];
rz(-0.8067185) q[3];
sx q[3];
rz(-0.66238856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0078916773) q[0];
sx q[0];
rz(-1.5971203) q[0];
sx q[0];
rz(1.9418465) q[0];
rz(-1.3613191) q[1];
sx q[1];
rz(-0.97313762) q[1];
sx q[1];
rz(2.5637085) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91980714) q[0];
sx q[0];
rz(-3.1263906) q[0];
sx q[0];
rz(-1.1342681) q[0];
rz(-pi) q[1];
rz(0.43283845) q[2];
sx q[2];
rz(-1.4269526) q[2];
sx q[2];
rz(1.8707616) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.51398811) q[1];
sx q[1];
rz(-2.0236694) q[1];
sx q[1];
rz(-1.7220366) q[1];
rz(3.1158524) q[3];
sx q[3];
rz(-1.2719858) q[3];
sx q[3];
rz(-0.65626496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63152385) q[2];
sx q[2];
rz(-0.57454595) q[2];
sx q[2];
rz(-2.6681382) q[2];
rz(1.8691285) q[3];
sx q[3];
rz(-1.748964) q[3];
sx q[3];
rz(2.9191391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1943787) q[0];
sx q[0];
rz(-2.4838303) q[0];
sx q[0];
rz(-0.33313242) q[0];
rz(-2.9130452) q[1];
sx q[1];
rz(-1.8014149) q[1];
sx q[1];
rz(-2.5659335) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7299907) q[0];
sx q[0];
rz(-1.6766929) q[0];
sx q[0];
rz(-2.3292755) q[0];
rz(-pi) q[1];
rz(0.47708738) q[2];
sx q[2];
rz(-1.2735575) q[2];
sx q[2];
rz(1.5971668) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.595069) q[1];
sx q[1];
rz(-1.8966339) q[1];
sx q[1];
rz(-0.21104668) q[1];
rz(-pi) q[2];
rz(1.3153399) q[3];
sx q[3];
rz(-1.1029579) q[3];
sx q[3];
rz(-1.5327191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.26611844) q[2];
sx q[2];
rz(-0.55142752) q[2];
sx q[2];
rz(1.4761338) q[2];
rz(-1.1009781) q[3];
sx q[3];
rz(-1.063238) q[3];
sx q[3];
rz(2.6235918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(1.5440893) q[0];
sx q[0];
rz(-2.2552555) q[0];
sx q[0];
rz(-0.30712095) q[0];
rz(-2.8111474) q[1];
sx q[1];
rz(-1.5437061) q[1];
sx q[1];
rz(1.365136) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86548818) q[0];
sx q[0];
rz(-1.5886602) q[0];
sx q[0];
rz(-3.0602656) q[0];
rz(-pi) q[1];
rz(2.6583985) q[2];
sx q[2];
rz(-0.38887923) q[2];
sx q[2];
rz(0.70975297) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0338034) q[1];
sx q[1];
rz(-2.6328936) q[1];
sx q[1];
rz(-2.1230039) q[1];
rz(-pi) q[2];
rz(0.14162986) q[3];
sx q[3];
rz(-2.3034769) q[3];
sx q[3];
rz(-1.4328626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.27552989) q[2];
sx q[2];
rz(-0.58774647) q[2];
sx q[2];
rz(0.24869871) q[2];
rz(1.2134877) q[3];
sx q[3];
rz(-1.4444193) q[3];
sx q[3];
rz(0.061633751) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98951775) q[0];
sx q[0];
rz(-0.9062506) q[0];
sx q[0];
rz(0.686598) q[0];
rz(2.0286512) q[1];
sx q[1];
rz(-2.3390892) q[1];
sx q[1];
rz(-2.8259605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.21232) q[0];
sx q[0];
rz(-1.2480191) q[0];
sx q[0];
rz(-3.0930711) q[0];
rz(-pi) q[1];
rz(1.8191843) q[2];
sx q[2];
rz(-2.3150716) q[2];
sx q[2];
rz(0.404111) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0084969) q[1];
sx q[1];
rz(-1.4766472) q[1];
sx q[1];
rz(-2.2341683) q[1];
rz(0.078048869) q[3];
sx q[3];
rz(-2.4025318) q[3];
sx q[3];
rz(-0.8821677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.45372066) q[2];
sx q[2];
rz(-2.5747955) q[2];
sx q[2];
rz(0.68457121) q[2];
rz(1.941393) q[3];
sx q[3];
rz(-1.129351) q[3];
sx q[3];
rz(2.9620192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0612563) q[0];
sx q[0];
rz(-0.48114023) q[0];
sx q[0];
rz(3.1367593) q[0];
rz(-1.5176516) q[1];
sx q[1];
rz(-2.3976517) q[1];
sx q[1];
rz(2.8878816) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1954945) q[0];
sx q[0];
rz(-1.4854447) q[0];
sx q[0];
rz(-1.781989) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4058622) q[2];
sx q[2];
rz(-1.2085739) q[2];
sx q[2];
rz(-0.74184201) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1349134) q[1];
sx q[1];
rz(-2.2354534) q[1];
sx q[1];
rz(2.8921739) q[1];
x q[2];
rz(-1.3718448) q[3];
sx q[3];
rz(-1.0291489) q[3];
sx q[3];
rz(-1.0197848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1634875) q[2];
sx q[2];
rz(-1.9065964) q[2];
sx q[2];
rz(1.9592436) q[2];
rz(-0.67665082) q[3];
sx q[3];
rz(-2.7550321) q[3];
sx q[3];
rz(2.1277908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5273298) q[0];
sx q[0];
rz(-2.3557721) q[0];
sx q[0];
rz(-2.8489805) q[0];
rz(-1.0489427) q[1];
sx q[1];
rz(-1.7221778) q[1];
sx q[1];
rz(-2.4377951) q[1];
rz(0.47164698) q[2];
sx q[2];
rz(-1.2670965) q[2];
sx q[2];
rz(-0.22128174) q[2];
rz(1.9378035) q[3];
sx q[3];
rz(-1.6631775) q[3];
sx q[3];
rz(0.1803995) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
