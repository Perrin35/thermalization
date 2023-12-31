OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8653712) q[0];
sx q[0];
rz(-2.2844391) q[0];
sx q[0];
rz(-0.13248086) q[0];
rz(-2.8744856) q[1];
sx q[1];
rz(-2.5565956) q[1];
sx q[1];
rz(-2.4490228) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8938457) q[0];
sx q[0];
rz(-0.5286628) q[0];
sx q[0];
rz(-1.7696487) q[0];
rz(2.0376671) q[2];
sx q[2];
rz(-0.97969998) q[2];
sx q[2];
rz(0.95962722) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.2107271) q[1];
sx q[1];
rz(-2.9949246) q[1];
sx q[1];
rz(2.0861097) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0308706) q[3];
sx q[3];
rz(-2.6220136) q[3];
sx q[3];
rz(-0.67668623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1074368) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(-1.5365323) q[2];
rz(-1.5213373) q[3];
sx q[3];
rz(-1.6531569) q[3];
sx q[3];
rz(0.03604123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.81543106) q[0];
sx q[0];
rz(-0.27359971) q[0];
sx q[0];
rz(1.8923627) q[0];
rz(-0.56150395) q[1];
sx q[1];
rz(-2.3655472) q[1];
sx q[1];
rz(0.5805648) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8037572) q[0];
sx q[0];
rz(-0.32807402) q[0];
sx q[0];
rz(-2.8003545) q[0];
rz(-1.3524019) q[2];
sx q[2];
rz(-2.1973655) q[2];
sx q[2];
rz(0.44109694) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7932574) q[1];
sx q[1];
rz(-1.0151334) q[1];
sx q[1];
rz(-0.72850119) q[1];
x q[2];
rz(2.21653) q[3];
sx q[3];
rz(-2.5930773) q[3];
sx q[3];
rz(-0.95228449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4124174) q[2];
sx q[2];
rz(-0.20755945) q[2];
sx q[2];
rz(-0.87835971) q[2];
rz(0.39204028) q[3];
sx q[3];
rz(-1.6974028) q[3];
sx q[3];
rz(-0.6033321) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47675258) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(-0.77366775) q[0];
rz(0.0013008612) q[1];
sx q[1];
rz(-1.6157849) q[1];
sx q[1];
rz(3.1087648) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87356991) q[0];
sx q[0];
rz(-2.7497254) q[0];
sx q[0];
rz(-2.5001532) q[0];
rz(-2.122934) q[2];
sx q[2];
rz(-1.3214006) q[2];
sx q[2];
rz(2.2160335) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3120011) q[1];
sx q[1];
rz(-0.81273505) q[1];
sx q[1];
rz(2.9144822) q[1];
rz(-pi) q[2];
rz(-2.7955416) q[3];
sx q[3];
rz(-1.0533353) q[3];
sx q[3];
rz(0.82511653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34439987) q[2];
sx q[2];
rz(-1.025528) q[2];
sx q[2];
rz(2.8642505) q[2];
rz(2.7456361) q[3];
sx q[3];
rz(-1.5405416) q[3];
sx q[3];
rz(-2.4424281) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70401496) q[0];
sx q[0];
rz(-0.33575785) q[0];
sx q[0];
rz(-2.8787956) q[0];
rz(0.2335877) q[1];
sx q[1];
rz(-0.83507744) q[1];
sx q[1];
rz(-0.77082005) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0924776) q[0];
sx q[0];
rz(-3.1160833) q[0];
sx q[0];
rz(-2.269948) q[0];
x q[1];
rz(1.0806482) q[2];
sx q[2];
rz(-1.4785826) q[2];
sx q[2];
rz(-0.25445081) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1781436) q[1];
sx q[1];
rz(-1.2091067) q[1];
sx q[1];
rz(2.5835035) q[1];
x q[2];
rz(-0.869107) q[3];
sx q[3];
rz(-2.8512555) q[3];
sx q[3];
rz(1.7770191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16584855) q[2];
sx q[2];
rz(-1.6341012) q[2];
sx q[2];
rz(-2.4285994) q[2];
rz(-1.0130079) q[3];
sx q[3];
rz(-0.37390798) q[3];
sx q[3];
rz(0.95389429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7825496) q[0];
sx q[0];
rz(-1.0721711) q[0];
sx q[0];
rz(-1.7011401) q[0];
rz(-3.0474995) q[1];
sx q[1];
rz(-0.73939878) q[1];
sx q[1];
rz(-2.9715911) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6572896) q[0];
sx q[0];
rz(-2.347725) q[0];
sx q[0];
rz(0.33052175) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0988118) q[2];
sx q[2];
rz(-1.3655647) q[2];
sx q[2];
rz(0.99265487) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1483037) q[1];
sx q[1];
rz(-1.9699886) q[1];
sx q[1];
rz(0.94435512) q[1];
rz(-pi) q[2];
rz(2.2682297) q[3];
sx q[3];
rz(-2.1057099) q[3];
sx q[3];
rz(1.8148592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1468982) q[2];
sx q[2];
rz(-2.6165104) q[2];
sx q[2];
rz(-1.7374932) q[2];
rz(-1.5480301) q[3];
sx q[3];
rz(-0.78909767) q[3];
sx q[3];
rz(1.4353969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65574044) q[0];
sx q[0];
rz(-1.1463373) q[0];
sx q[0];
rz(-2.6830542) q[0];
rz(-2.8857152) q[1];
sx q[1];
rz(-1.8829388) q[1];
sx q[1];
rz(-2.4564254) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72502575) q[0];
sx q[0];
rz(-1.6402813) q[0];
sx q[0];
rz(2.2216703) q[0];
rz(2.3643656) q[2];
sx q[2];
rz(-1.4677375) q[2];
sx q[2];
rz(1.4003786) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.203478) q[1];
sx q[1];
rz(-2.1739829) q[1];
sx q[1];
rz(0.60738648) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1232883) q[3];
sx q[3];
rz(-0.98494512) q[3];
sx q[3];
rz(2.626112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.3926065) q[2];
sx q[2];
rz(-0.41967732) q[2];
sx q[2];
rz(-1.4292599) q[2];
rz(-1.0990934) q[3];
sx q[3];
rz(-2.6350239) q[3];
sx q[3];
rz(-2.9523622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1290865) q[0];
sx q[0];
rz(-1.5456454) q[0];
sx q[0];
rz(-2.4122453) q[0];
rz(-2.8485281) q[1];
sx q[1];
rz(-2.9022419) q[1];
sx q[1];
rz(-1.1475295) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4097737) q[0];
sx q[0];
rz(-1.496135) q[0];
sx q[0];
rz(2.7274107) q[0];
x q[1];
rz(2.7717934) q[2];
sx q[2];
rz(-2.1714006) q[2];
sx q[2];
rz(2.415654) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7123375) q[1];
sx q[1];
rz(-2.0550248) q[1];
sx q[1];
rz(-1.5244563) q[1];
rz(-pi) q[2];
rz(0.64671867) q[3];
sx q[3];
rz(-2.5463856) q[3];
sx q[3];
rz(-1.7531542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3532233) q[2];
sx q[2];
rz(-2.1187783) q[2];
sx q[2];
rz(-2.3925171) q[2];
rz(-2.4979112) q[3];
sx q[3];
rz(-1.0130853) q[3];
sx q[3];
rz(2.2275887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-0.81326905) q[1];
sx q[1];
rz(-0.39852279) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9015394) q[0];
sx q[0];
rz(-1.1050637) q[0];
sx q[0];
rz(2.9630911) q[0];
rz(1.6767098) q[2];
sx q[2];
rz(-2.1692861) q[2];
sx q[2];
rz(0.61818365) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1426705) q[1];
sx q[1];
rz(-2.167836) q[1];
sx q[1];
rz(0.02762694) q[1];
rz(-1.0141482) q[3];
sx q[3];
rz(-1.1245407) q[3];
sx q[3];
rz(-0.92170148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1901671) q[2];
sx q[2];
rz(-1.3484893) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012638906) q[0];
sx q[0];
rz(-2.5057827) q[0];
sx q[0];
rz(1.7522316) q[0];
rz(1.6268436) q[1];
sx q[1];
rz(-1.4667958) q[1];
sx q[1];
rz(-1.0983889) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2929045) q[0];
sx q[0];
rz(-1.0875889) q[0];
sx q[0];
rz(0.26341652) q[0];
rz(-pi) q[1];
rz(-0.45778747) q[2];
sx q[2];
rz(-0.87399235) q[2];
sx q[2];
rz(-0.84504499) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2201982) q[1];
sx q[1];
rz(-0.99153334) q[1];
sx q[1];
rz(1.571435) q[1];
rz(-pi) q[2];
rz(0.58795712) q[3];
sx q[3];
rz(-0.30246099) q[3];
sx q[3];
rz(1.0126225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2417458) q[2];
sx q[2];
rz(-1.2925623) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.39524233) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(-0.46646068) q[0];
rz(2.9699504) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(-2.5126273) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2254612) q[0];
sx q[0];
rz(-2.7636508) q[0];
sx q[0];
rz(-1.4378689) q[0];
rz(-pi) q[1];
rz(0.92832698) q[2];
sx q[2];
rz(-1.3767585) q[2];
sx q[2];
rz(2.5460668) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.65044636) q[1];
sx q[1];
rz(-1.5728587) q[1];
sx q[1];
rz(0.44209977) q[1];
rz(-1.2260776) q[3];
sx q[3];
rz(-2.260672) q[3];
sx q[3];
rz(-0.41050875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24370596) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(2.0824599) q[2];
rz(-2.4641666) q[3];
sx q[3];
rz(-0.99223653) q[3];
sx q[3];
rz(0.89390755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8582936) q[0];
sx q[0];
rz(-1.1245921) q[0];
sx q[0];
rz(-0.84381214) q[0];
rz(-2.8876866) q[1];
sx q[1];
rz(-1.084068) q[1];
sx q[1];
rz(-0.57938309) q[1];
rz(2.4181096) q[2];
sx q[2];
rz(-1.6037446) q[2];
sx q[2];
rz(-1.4708191) q[2];
rz(-0.068594882) q[3];
sx q[3];
rz(-2.1566236) q[3];
sx q[3];
rz(-1.9038283) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
