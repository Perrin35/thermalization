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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8938457) q[0];
sx q[0];
rz(-2.6129299) q[0];
sx q[0];
rz(1.7696487) q[0];
rz(0.59074596) q[2];
sx q[2];
rz(-2.4060537) q[2];
sx q[2];
rz(0.22434805) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9308656) q[1];
sx q[1];
rz(-2.9949246) q[1];
sx q[1];
rz(1.055483) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2487189) q[3];
sx q[3];
rz(-1.1097483) q[3];
sx q[3];
rz(-2.98364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0341558) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(1.5365323) q[2];
rz(1.5213373) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(0.03604123) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81543106) q[0];
sx q[0];
rz(-2.8679929) q[0];
sx q[0];
rz(-1.8923627) q[0];
rz(-0.56150395) q[1];
sx q[1];
rz(-0.77604547) q[1];
sx q[1];
rz(2.5610279) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97894325) q[0];
sx q[0];
rz(-1.2622841) q[0];
sx q[0];
rz(-1.6842151) q[0];
rz(-pi) q[1];
rz(2.5035628) q[2];
sx q[2];
rz(-1.3943765) q[2];
sx q[2];
rz(1.0002913) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7932574) q[1];
sx q[1];
rz(-1.0151334) q[1];
sx q[1];
rz(-0.72850119) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7892116) q[3];
sx q[3];
rz(-1.1412732) q[3];
sx q[3];
rz(-2.9126715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7291752) q[2];
sx q[2];
rz(-2.9340332) q[2];
sx q[2];
rz(0.87835971) q[2];
rz(2.7495524) q[3];
sx q[3];
rz(-1.6974028) q[3];
sx q[3];
rz(0.6033321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47675258) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(2.3679249) q[0];
rz(3.1402918) q[1];
sx q[1];
rz(-1.5258077) q[1];
sx q[1];
rz(3.1087648) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5531909) q[0];
sx q[0];
rz(-1.8817888) q[0];
sx q[0];
rz(-1.8131959) q[0];
x q[1];
rz(0.29067729) q[2];
sx q[2];
rz(-2.103984) q[2];
sx q[2];
rz(2.3454587) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58363885) q[1];
sx q[1];
rz(-1.4065521) q[1];
sx q[1];
rz(-0.79973952) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.114931) q[3];
sx q[3];
rz(-1.2715724) q[3];
sx q[3];
rz(-0.92218375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34439987) q[2];
sx q[2];
rz(-1.025528) q[2];
sx q[2];
rz(2.8642505) q[2];
rz(-0.39595655) q[3];
sx q[3];
rz(-1.6010511) q[3];
sx q[3];
rz(2.4424281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4375777) q[0];
sx q[0];
rz(-0.33575785) q[0];
sx q[0];
rz(-0.26279703) q[0];
rz(2.908005) q[1];
sx q[1];
rz(-2.3065152) q[1];
sx q[1];
rz(2.3707726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.049115) q[0];
sx q[0];
rz(-3.1160833) q[0];
sx q[0];
rz(2.269948) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3768251) q[2];
sx q[2];
rz(-0.49805222) q[2];
sx q[2];
rz(-1.1454524) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1781436) q[1];
sx q[1];
rz(-1.2091067) q[1];
sx q[1];
rz(-2.5835035) q[1];
x q[2];
rz(-1.3464438) q[3];
sx q[3];
rz(-1.3849349) q[3];
sx q[3];
rz(-0.88691521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.16584855) q[2];
sx q[2];
rz(-1.5074915) q[2];
sx q[2];
rz(-0.7129933) q[2];
rz(1.0130079) q[3];
sx q[3];
rz(-2.7676847) q[3];
sx q[3];
rz(0.95389429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35904303) q[0];
sx q[0];
rz(-2.0694216) q[0];
sx q[0];
rz(-1.4404526) q[0];
rz(3.0474995) q[1];
sx q[1];
rz(-0.73939878) q[1];
sx q[1];
rz(2.9715911) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6572896) q[0];
sx q[0];
rz(-0.7938677) q[0];
sx q[0];
rz(-0.33052175) q[0];
x q[1];
rz(-2.0427809) q[2];
sx q[2];
rz(-1.3655647) q[2];
sx q[2];
rz(2.1489378) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0712142) q[1];
sx q[1];
rz(-0.7281248) q[1];
sx q[1];
rz(-0.94707625) q[1];
rz(-pi) q[2];
rz(2.3159536) q[3];
sx q[3];
rz(-0.85074433) q[3];
sx q[3];
rz(-2.839098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1468982) q[2];
sx q[2];
rz(-0.52508223) q[2];
sx q[2];
rz(1.7374932) q[2];
rz(1.5480301) q[3];
sx q[3];
rz(-2.352495) q[3];
sx q[3];
rz(1.4353969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65574044) q[0];
sx q[0];
rz(-1.1463373) q[0];
sx q[0];
rz(2.6830542) q[0];
rz(-2.8857152) q[1];
sx q[1];
rz(-1.2586539) q[1];
sx q[1];
rz(2.4564254) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4165669) q[0];
sx q[0];
rz(-1.6402813) q[0];
sx q[0];
rz(2.2216703) q[0];
x q[1];
rz(2.9951727) q[2];
sx q[2];
rz(-2.358987) q[2];
sx q[2];
rz(-2.866982) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8243858) q[1];
sx q[1];
rz(-2.3134391) q[1];
sx q[1];
rz(-0.87888996) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98486793) q[3];
sx q[3];
rz(-1.586048) q[3];
sx q[3];
rz(1.0451942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7489862) q[2];
sx q[2];
rz(-0.41967732) q[2];
sx q[2];
rz(1.4292599) q[2];
rz(2.0424992) q[3];
sx q[3];
rz(-2.6350239) q[3];
sx q[3];
rz(0.18923047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1290865) q[0];
sx q[0];
rz(-1.5959473) q[0];
sx q[0];
rz(2.4122453) q[0];
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
rz(0.731819) q[0];
sx q[0];
rz(-1.6454576) q[0];
sx q[0];
rz(-2.7274107) q[0];
rz(-pi) q[1];
rz(-2.2044264) q[2];
sx q[2];
rz(-1.8735777) q[2];
sx q[2];
rz(-0.62919754) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.16312576) q[1];
sx q[1];
rz(-1.5297869) q[1];
sx q[1];
rz(-0.48467111) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.494874) q[3];
sx q[3];
rz(-0.59520703) q[3];
sx q[3];
rz(1.7531542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.78836936) q[2];
sx q[2];
rz(-2.1187783) q[2];
sx q[2];
rz(2.3925171) q[2];
rz(-0.64368147) q[3];
sx q[3];
rz(-1.0130853) q[3];
sx q[3];
rz(-2.2275887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1483243) q[0];
sx q[0];
rz(-2.0355621) q[0];
sx q[0];
rz(2.3102982) q[0];
rz(1.7656901) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(2.7430699) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7300028) q[0];
sx q[0];
rz(-1.7301136) q[0];
sx q[0];
rz(-2.0429862) q[0];
rz(-pi) q[1];
rz(1.6767098) q[2];
sx q[2];
rz(-2.1692861) q[2];
sx q[2];
rz(-2.523409) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.41259137) q[1];
sx q[1];
rz(-1.593643) q[1];
sx q[1];
rz(0.97357915) q[1];
rz(-1.0141482) q[3];
sx q[3];
rz(-1.1245407) q[3];
sx q[3];
rz(-0.92170148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1901671) q[2];
sx q[2];
rz(-1.3484893) q[2];
sx q[2];
rz(1.8438967) q[2];
rz(-2.0166345) q[3];
sx q[3];
rz(-1.8816032) q[3];
sx q[3];
rz(2.4979533) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1289537) q[0];
sx q[0];
rz(-0.63580996) q[0];
sx q[0];
rz(1.3893611) q[0];
rz(1.6268436) q[1];
sx q[1];
rz(-1.6747968) q[1];
sx q[1];
rz(-2.0432037) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59745715) q[0];
sx q[0];
rz(-1.8034593) q[0];
sx q[0];
rz(2.068589) q[0];
rz(-2.3214925) q[2];
sx q[2];
rz(-1.9165877) q[2];
sx q[2];
rz(-2.1095914) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2201982) q[1];
sx q[1];
rz(-0.99153334) q[1];
sx q[1];
rz(-1.571435) q[1];
x q[2];
rz(2.5536355) q[3];
sx q[3];
rz(-2.8391317) q[3];
sx q[3];
rz(1.0126225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2417458) q[2];
sx q[2];
rz(-1.8490303) q[2];
sx q[2];
rz(0.51952726) q[2];
rz(-1.3351006) q[3];
sx q[3];
rz(-2.3050008) q[3];
sx q[3];
rz(-1.8241204) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39524233) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(-2.675132) q[0];
rz(2.9699504) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(-2.5126273) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2254612) q[0];
sx q[0];
rz(-0.3779419) q[0];
sx q[0];
rz(1.7037237) q[0];
rz(-pi) q[1];
rz(-1.8877108) q[2];
sx q[2];
rz(-2.4744518) q[2];
sx q[2];
rz(2.4184879) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.92470691) q[1];
sx q[1];
rz(-0.44210426) q[1];
sx q[1];
rz(0.0048203992) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71979379) q[3];
sx q[3];
rz(-1.8344973) q[3];
sx q[3];
rz(1.3849474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8978867) q[2];
sx q[2];
rz(-1.0906929) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28329904) q[0];
sx q[0];
rz(-2.0170006) q[0];
sx q[0];
rz(2.2977805) q[0];
rz(-2.8876866) q[1];
sx q[1];
rz(-1.084068) q[1];
sx q[1];
rz(-0.57938309) q[1];
rz(-0.72348307) q[2];
sx q[2];
rz(-1.6037446) q[2];
sx q[2];
rz(-1.4708191) q[2];
rz(0.068594882) q[3];
sx q[3];
rz(-0.98496901) q[3];
sx q[3];
rz(1.2377644) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
