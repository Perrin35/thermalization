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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.010095) q[0];
sx q[0];
rz(-0.37535497) q[0];
sx q[0];
rz(-1.4913056) q[0];
rz(-pi) q[1];
rz(-1.5145095) q[2];
sx q[2];
rz(-2.0234642) q[2];
sx q[2];
rz(-2.2279091) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1387922) q[1];
sx q[1];
rz(-0.82502581) q[1];
sx q[1];
rz(1.1345798) q[1];
x q[2];
rz(0.31422796) q[3];
sx q[3];
rz(-2.3433821) q[3];
sx q[3];
rz(1.2558921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1859833) q[2];
sx q[2];
rz(-0.60638967) q[2];
sx q[2];
rz(-2.8601904) q[2];
rz(1.9143117) q[3];
sx q[3];
rz(-1.3325007) q[3];
sx q[3];
rz(-1.2138155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94673741) q[0];
sx q[0];
rz(-2.7342716) q[0];
sx q[0];
rz(0.58709225) q[0];
rz(-0.52014822) q[1];
sx q[1];
rz(-1.2112434) q[1];
sx q[1];
rz(0.21336666) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086047202) q[0];
sx q[0];
rz(-0.16250691) q[0];
sx q[0];
rz(1.1799728) q[0];
rz(-2.6293829) q[2];
sx q[2];
rz(-2.1551415) q[2];
sx q[2];
rz(2.0456631) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0679464) q[1];
sx q[1];
rz(-1.526462) q[1];
sx q[1];
rz(1.6641875) q[1];
rz(1.1150241) q[3];
sx q[3];
rz(-2.5682862) q[3];
sx q[3];
rz(0.70070964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0352036) q[2];
sx q[2];
rz(-2.8530402) q[2];
sx q[2];
rz(-2.8698548) q[2];
rz(-0.64618293) q[3];
sx q[3];
rz(-1.8281507) q[3];
sx q[3];
rz(1.2077695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.81994098) q[0];
sx q[0];
rz(-0.8701179) q[0];
sx q[0];
rz(0.2680378) q[0];
rz(0.23621121) q[1];
sx q[1];
rz(-1.3583438) q[1];
sx q[1];
rz(1.7265629) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1911604) q[0];
sx q[0];
rz(-1.6163905) q[0];
sx q[0];
rz(3.0395503) q[0];
x q[1];
rz(2.7790952) q[2];
sx q[2];
rz(-1.4115184) q[2];
sx q[2];
rz(2.8822219) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.336795) q[1];
sx q[1];
rz(-2.6141254) q[1];
sx q[1];
rz(0.57425604) q[1];
x q[2];
rz(0.3424267) q[3];
sx q[3];
rz(-0.70243109) q[3];
sx q[3];
rz(2.680955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1872824) q[2];
sx q[2];
rz(-0.3287181) q[2];
sx q[2];
rz(1.3817361) q[2];
rz(0.36661822) q[3];
sx q[3];
rz(-1.6691875) q[3];
sx q[3];
rz(-0.097675145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5919375) q[0];
sx q[0];
rz(-2.9357935) q[0];
sx q[0];
rz(-0.016121443) q[0];
rz(1.4664117) q[1];
sx q[1];
rz(-1.6203531) q[1];
sx q[1];
rz(-0.88502562) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2766874) q[0];
sx q[0];
rz(-1.7430796) q[0];
sx q[0];
rz(2.7322463) q[0];
rz(-pi) q[1];
rz(0.22815223) q[2];
sx q[2];
rz(-1.9288262) q[2];
sx q[2];
rz(1.6659074) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4440892) q[1];
sx q[1];
rz(-0.61629811) q[1];
sx q[1];
rz(0.69774903) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0543941) q[3];
sx q[3];
rz(-1.2149016) q[3];
sx q[3];
rz(0.10117029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5625988) q[2];
sx q[2];
rz(-0.4563953) q[2];
sx q[2];
rz(-1.558051) q[2];
rz(1.3715749) q[3];
sx q[3];
rz(-1.8669502) q[3];
sx q[3];
rz(2.9461327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9850995) q[0];
sx q[0];
rz(-1.3753966) q[0];
sx q[0];
rz(0.14955713) q[0];
rz(1.043148) q[1];
sx q[1];
rz(-2.1688192) q[1];
sx q[1];
rz(0.8358039) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1096032) q[0];
sx q[0];
rz(-1.2471203) q[0];
sx q[0];
rz(0.26357414) q[0];
rz(-0.0094532254) q[2];
sx q[2];
rz(-0.88047709) q[2];
sx q[2];
rz(-3.0005665) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0975102) q[1];
sx q[1];
rz(-2.9335576) q[1];
sx q[1];
rz(1.5670304) q[1];
x q[2];
rz(2.9582865) q[3];
sx q[3];
rz(-1.8170895) q[3];
sx q[3];
rz(3.1251647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0189556) q[2];
sx q[2];
rz(-1.6654207) q[2];
sx q[2];
rz(2.7254851) q[2];
rz(-3.05919) q[3];
sx q[3];
rz(-2.1722138) q[3];
sx q[3];
rz(2.940322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06269726) q[0];
sx q[0];
rz(-1.0286999) q[0];
sx q[0];
rz(-2.0507574) q[0];
rz(1.1385607) q[1];
sx q[1];
rz(-1.9416315) q[1];
sx q[1];
rz(2.535215) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67865506) q[0];
sx q[0];
rz(-2.444164) q[0];
sx q[0];
rz(0.3987958) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1101704) q[2];
sx q[2];
rz(-2.3177882) q[2];
sx q[2];
rz(0.065814171) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5910037) q[1];
sx q[1];
rz(-0.85500756) q[1];
sx q[1];
rz(2.4229933) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0020671) q[3];
sx q[3];
rz(-0.72318132) q[3];
sx q[3];
rz(-0.65825247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.31327569) q[2];
sx q[2];
rz(-2.5930391) q[2];
sx q[2];
rz(-0.11087785) q[2];
rz(-1.5634792) q[3];
sx q[3];
rz(-2.9302247) q[3];
sx q[3];
rz(-1.3255239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.9368847) q[0];
sx q[0];
rz(-2.3218343) q[0];
sx q[0];
rz(2.4844266) q[0];
rz(1.8183297) q[1];
sx q[1];
rz(-0.75138775) q[1];
sx q[1];
rz(2.0164067) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.781625) q[0];
sx q[0];
rz(-0.021718135) q[0];
sx q[0];
rz(0.75628539) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29203307) q[2];
sx q[2];
rz(-2.8148373) q[2];
sx q[2];
rz(-2.9346443) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4309498) q[1];
sx q[1];
rz(-0.64326566) q[1];
sx q[1];
rz(-2.4673854) q[1];
rz(0.67496211) q[3];
sx q[3];
rz(-1.3276023) q[3];
sx q[3];
rz(1.7053108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40921673) q[2];
sx q[2];
rz(-1.9595307) q[2];
sx q[2];
rz(-2.1803975) q[2];
rz(-3.1327278) q[3];
sx q[3];
rz(-0.88518393) q[3];
sx q[3];
rz(1.2113021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0601785) q[0];
sx q[0];
rz(-0.1135122) q[0];
sx q[0];
rz(0.42763448) q[0];
rz(-2.0291406) q[1];
sx q[1];
rz(-1.8266269) q[1];
sx q[1];
rz(-2.1449259) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8450981) q[0];
sx q[0];
rz(-1.0834435) q[0];
sx q[0];
rz(-0.71655207) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3786167) q[2];
sx q[2];
rz(-1.7691649) q[2];
sx q[2];
rz(1.5114443) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5298047) q[1];
sx q[1];
rz(-2.5080097) q[1];
sx q[1];
rz(1.1128913) q[1];
rz(-pi) q[2];
rz(1.7782042) q[3];
sx q[3];
rz(-1.0098352) q[3];
sx q[3];
rz(0.029766729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4619649) q[2];
sx q[2];
rz(-2.6042016) q[2];
sx q[2];
rz(0.648554) q[2];
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
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40470966) q[0];
sx q[0];
rz(-2.9209904) q[0];
sx q[0];
rz(-0.96694651) q[0];
rz(-0.95775682) q[1];
sx q[1];
rz(-1.0916595) q[1];
sx q[1];
rz(-2.0584094) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082302467) q[0];
sx q[0];
rz(-1.8169615) q[0];
sx q[0];
rz(0.44025006) q[0];
rz(-pi) q[1];
rz(-0.71523209) q[2];
sx q[2];
rz(-2.2713647) q[2];
sx q[2];
rz(2.689555) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8260767) q[1];
sx q[1];
rz(-0.89102645) q[1];
sx q[1];
rz(-2.6273374) q[1];
rz(-pi) q[2];
rz(-1.5708709) q[3];
sx q[3];
rz(-1.0914472) q[3];
sx q[3];
rz(-1.0842931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0208685) q[2];
sx q[2];
rz(-2.3126297) q[2];
sx q[2];
rz(-1.6797569) q[2];
rz(2.7802137) q[3];
sx q[3];
rz(-1.3249818) q[3];
sx q[3];
rz(2.1593275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7892889) q[0];
sx q[0];
rz(-2.4109349) q[0];
sx q[0];
rz(2.5590382) q[0];
rz(1.3652868) q[1];
sx q[1];
rz(-1.0368232) q[1];
sx q[1];
rz(0.22886151) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34535566) q[0];
sx q[0];
rz(-1.3889878) q[0];
sx q[0];
rz(-0.40248351) q[0];
rz(-2.1624915) q[2];
sx q[2];
rz(-2.8276463) q[2];
sx q[2];
rz(-0.57127956) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4881674) q[1];
sx q[1];
rz(-1.5060802) q[1];
sx q[1];
rz(1.72987) q[1];
rz(-pi) q[2];
rz(1.1353605) q[3];
sx q[3];
rz(-1.9276016) q[3];
sx q[3];
rz(0.26998587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.87054306) q[2];
sx q[2];
rz(-0.55176631) q[2];
sx q[2];
rz(1.40847) q[2];
rz(-2.4208141) q[3];
sx q[3];
rz(-1.9920789) q[3];
sx q[3];
rz(1.518998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.0558753) q[0];
sx q[0];
rz(-1.4376823) q[0];
sx q[0];
rz(1.8327376) q[0];
rz(1.5695288) q[1];
sx q[1];
rz(-2.5253898) q[1];
sx q[1];
rz(0.22402221) q[1];
rz(-0.24893473) q[2];
sx q[2];
rz(-1.0200844) q[2];
sx q[2];
rz(3.040584) q[2];
rz(-0.15275501) q[3];
sx q[3];
rz(-1.7293617) q[3];
sx q[3];
rz(0.21241906) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
