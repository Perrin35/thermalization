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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1231175) q[0];
sx q[0];
rz(-1.0536061) q[0];
sx q[0];
rz(0.11488199) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0376671) q[2];
sx q[2];
rz(-0.97969998) q[2];
sx q[2];
rz(0.95962722) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9308656) q[1];
sx q[1];
rz(-2.9949246) q[1];
sx q[1];
rz(-1.055483) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2487189) q[3];
sx q[3];
rz(-1.1097483) q[3];
sx q[3];
rz(0.15795262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1074368) q[2];
sx q[2];
rz(-2.6145356) q[2];
sx q[2];
rz(-1.6050603) q[2];
rz(1.5213373) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(-3.1055514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81543106) q[0];
sx q[0];
rz(-2.8679929) q[0];
sx q[0];
rz(-1.8923627) q[0];
rz(2.5800887) q[1];
sx q[1];
rz(-0.77604547) q[1];
sx q[1];
rz(2.5610279) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3378355) q[0];
sx q[0];
rz(-2.8135186) q[0];
sx q[0];
rz(2.8003545) q[0];
rz(-pi) q[1];
rz(2.5035628) q[2];
sx q[2];
rz(-1.7472162) q[2];
sx q[2];
rz(2.1413013) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3851871) q[1];
sx q[1];
rz(-2.2573973) q[1];
sx q[1];
rz(2.3910206) q[1];
rz(2.21653) q[3];
sx q[3];
rz(-2.5930773) q[3];
sx q[3];
rz(2.1893082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4124174) q[2];
sx q[2];
rz(-2.9340332) q[2];
sx q[2];
rz(2.2632329) q[2];
rz(-2.7495524) q[3];
sx q[3];
rz(-1.6974028) q[3];
sx q[3];
rz(-0.6033321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6648401) q[0];
sx q[0];
rz(-0.92233962) q[0];
sx q[0];
rz(-2.3679249) q[0];
rz(0.0013008612) q[1];
sx q[1];
rz(-1.6157849) q[1];
sx q[1];
rz(-0.032827854) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5884018) q[0];
sx q[0];
rz(-1.2598039) q[0];
sx q[0];
rz(-1.8131959) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1187395) q[2];
sx q[2];
rz(-0.60048044) q[2];
sx q[2];
rz(2.877419) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9877316) q[1];
sx q[1];
rz(-2.3567833) q[1];
sx q[1];
rz(-1.8042817) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34605108) q[3];
sx q[3];
rz(-2.0882574) q[3];
sx q[3];
rz(0.82511653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7971928) q[2];
sx q[2];
rz(-2.1160647) q[2];
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
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4375777) q[0];
sx q[0];
rz(-2.8058348) q[0];
sx q[0];
rz(-2.8787956) q[0];
rz(2.908005) q[1];
sx q[1];
rz(-0.83507744) q[1];
sx q[1];
rz(-2.3707726) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3931657) q[0];
sx q[0];
rz(-1.5512726) q[0];
sx q[0];
rz(-3.1251735) q[0];
rz(1.7647676) q[2];
sx q[2];
rz(-0.49805222) q[2];
sx q[2];
rz(1.9961403) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.091688823) q[1];
sx q[1];
rz(-0.65444512) q[1];
sx q[1];
rz(-2.5212538) q[1];
rz(0.1905258) q[3];
sx q[3];
rz(-1.3503721) q[3];
sx q[3];
rz(2.4998553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16584855) q[2];
sx q[2];
rz(-1.6341012) q[2];
sx q[2];
rz(-2.4285994) q[2];
rz(-1.0130079) q[3];
sx q[3];
rz(-0.37390798) q[3];
sx q[3];
rz(-2.1876984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(-0.17000155) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6572896) q[0];
sx q[0];
rz(-0.7938677) q[0];
sx q[0];
rz(-0.33052175) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0427809) q[2];
sx q[2];
rz(-1.3655647) q[2];
sx q[2];
rz(0.99265487) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4449094) q[1];
sx q[1];
rz(-1.0001567) q[1];
sx q[1];
rz(0.48008227) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4835028) q[3];
sx q[3];
rz(-2.1562025) q[3];
sx q[3];
rz(0.64774367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99469441) q[2];
sx q[2];
rz(-0.52508223) q[2];
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
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65574044) q[0];
sx q[0];
rz(-1.9952554) q[0];
sx q[0];
rz(-2.6830542) q[0];
rz(-0.25587747) q[1];
sx q[1];
rz(-1.8829388) q[1];
sx q[1];
rz(-0.68516723) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.348649) q[0];
sx q[0];
rz(-2.2198338) q[0];
sx q[0];
rz(-3.0543324) q[0];
x q[1];
rz(2.9951727) q[2];
sx q[2];
rz(-2.358987) q[2];
sx q[2];
rz(-2.866982) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.7428776) q[1];
sx q[1];
rz(-2.0600973) q[1];
sx q[1];
rz(2.2687885) q[1];
rz(-pi) q[2];
rz(2.1567247) q[3];
sx q[3];
rz(-1.5555447) q[3];
sx q[3];
rz(2.0963984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7489862) q[2];
sx q[2];
rz(-0.41967732) q[2];
sx q[2];
rz(1.4292599) q[2];
rz(-2.0424992) q[3];
sx q[3];
rz(-2.6350239) q[3];
sx q[3];
rz(2.9523622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(0.012506164) q[0];
sx q[0];
rz(-1.5959473) q[0];
sx q[0];
rz(-2.4122453) q[0];
rz(-2.8485281) q[1];
sx q[1];
rz(-2.9022419) q[1];
sx q[1];
rz(1.9940631) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4097737) q[0];
sx q[0];
rz(-1.6454576) q[0];
sx q[0];
rz(0.41418196) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36979923) q[2];
sx q[2];
rz(-0.97019201) q[2];
sx q[2];
rz(2.415654) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4292552) q[1];
sx q[1];
rz(-2.0550248) q[1];
sx q[1];
rz(1.6171364) q[1];
rz(-pi) q[2];
rz(2.494874) q[3];
sx q[3];
rz(-0.59520703) q[3];
sx q[3];
rz(1.3884384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.78836936) q[2];
sx q[2];
rz(-1.0228144) q[2];
sx q[2];
rz(2.3925171) q[2];
rz(0.64368147) q[3];
sx q[3];
rz(-2.1285074) q[3];
sx q[3];
rz(0.914004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1483243) q[0];
sx q[0];
rz(-2.0355621) q[0];
sx q[0];
rz(0.83129445) q[0];
rz(1.7656901) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(2.7430699) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9015394) q[0];
sx q[0];
rz(-1.1050637) q[0];
sx q[0];
rz(-2.9630911) q[0];
x q[1];
rz(-1.4648828) q[2];
sx q[2];
rz(-2.1692861) q[2];
sx q[2];
rz(0.61818365) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9989222) q[1];
sx q[1];
rz(-2.167836) q[1];
sx q[1];
rz(-0.02762694) q[1];
x q[2];
rz(-1.0141482) q[3];
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
rz(1.1901671) q[2];
sx q[2];
rz(-1.7931033) q[2];
sx q[2];
rz(-1.297696) q[2];
rz(-2.0166345) q[3];
sx q[3];
rz(-1.2599895) q[3];
sx q[3];
rz(0.64363939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1289537) q[0];
sx q[0];
rz(-0.63580996) q[0];
sx q[0];
rz(1.7522316) q[0];
rz(-1.6268436) q[1];
sx q[1];
rz(-1.6747968) q[1];
sx q[1];
rz(-1.0983889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59745715) q[0];
sx q[0];
rz(-1.8034593) q[0];
sx q[0];
rz(2.068589) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82010014) q[2];
sx q[2];
rz(-1.9165877) q[2];
sx q[2];
rz(2.1095914) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.92022773) q[1];
sx q[1];
rz(-0.5792633) q[1];
sx q[1];
rz(-0.00097640493) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3994201) q[3];
sx q[3];
rz(-1.3203353) q[3];
sx q[3];
rz(-2.7385538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8998469) q[2];
sx q[2];
rz(-1.2925623) q[2];
sx q[2];
rz(0.51952726) q[2];
rz(-1.3351006) q[3];
sx q[3];
rz(-2.3050008) q[3];
sx q[3];
rz(1.3174723) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7463503) q[0];
sx q[0];
rz(-1.1205751) q[0];
sx q[0];
rz(0.46646068) q[0];
rz(-2.9699504) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(2.5126273) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2254612) q[0];
sx q[0];
rz(-0.3779419) q[0];
sx q[0];
rz(1.4378689) q[0];
rz(-pi) q[1];
rz(-0.24068995) q[2];
sx q[2];
rz(-0.94229892) q[2];
sx q[2];
rz(2.0230053) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2222189) q[1];
sx q[1];
rz(-2.0128951) q[1];
sx q[1];
rz(1.5685146) q[1];
rz(-pi) q[2];
rz(-1.2260776) q[3];
sx q[3];
rz(-2.260672) q[3];
sx q[3];
rz(2.7310839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8978867) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(1.0591327) q[2];
rz(-2.4641666) q[3];
sx q[3];
rz(-2.1493561) q[3];
sx q[3];
rz(2.2476851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
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
rz(0.25390608) q[1];
sx q[1];
rz(-1.084068) q[1];
sx q[1];
rz(-0.57938309) q[1];
rz(-1.5268486) q[2];
sx q[2];
rz(-2.2938001) q[2];
sx q[2];
rz(-3.0707035) q[2];
rz(0.98388381) q[3];
sx q[3];
rz(-1.513653) q[3];
sx q[3];
rz(-0.37099864) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
