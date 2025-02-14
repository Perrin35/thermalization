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
rz(1.1310391) q[0];
sx q[0];
rz(-0.57589632) q[0];
sx q[0];
rz(-2.0146712) q[0];
rz(-0.95070401) q[1];
sx q[1];
rz(-0.60125142) q[1];
sx q[1];
rz(2.8079005) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5883049) q[0];
sx q[0];
rz(-1.9899564) q[0];
sx q[0];
rz(2.9774211) q[0];
x q[1];
rz(2.026171) q[2];
sx q[2];
rz(-1.3297594) q[2];
sx q[2];
rz(0.068403989) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6032927) q[1];
sx q[1];
rz(-1.2281907) q[1];
sx q[1];
rz(-1.3450772) q[1];
rz(-pi) q[2];
rz(0.26417664) q[3];
sx q[3];
rz(-0.81981711) q[3];
sx q[3];
rz(1.7163484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8814964) q[2];
sx q[2];
rz(-2.0659955) q[2];
sx q[2];
rz(2.0913701) q[2];
rz(-0.29551926) q[3];
sx q[3];
rz(-0.76885709) q[3];
sx q[3];
rz(1.3692921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1693901) q[0];
sx q[0];
rz(-1.0193595) q[0];
sx q[0];
rz(-0.66725677) q[0];
rz(-2.9908906) q[1];
sx q[1];
rz(-1.0548016) q[1];
sx q[1];
rz(-0.31164718) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2209476) q[0];
sx q[0];
rz(-1.049982) q[0];
sx q[0];
rz(0.08128165) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0097334) q[2];
sx q[2];
rz(-1.8510185) q[2];
sx q[2];
rz(0.2187905) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2296621) q[1];
sx q[1];
rz(-1.7102825) q[1];
sx q[1];
rz(0.18499495) q[1];
x q[2];
rz(2.7439609) q[3];
sx q[3];
rz(-1.5036229) q[3];
sx q[3];
rz(2.6844418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.47505891) q[2];
sx q[2];
rz(-2.2481613) q[2];
sx q[2];
rz(0.7676355) q[2];
rz(-2.660699) q[3];
sx q[3];
rz(-2.3700263) q[3];
sx q[3];
rz(-1.5868928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2956706) q[0];
sx q[0];
rz(-1.5376115) q[0];
sx q[0];
rz(-0.77958244) q[0];
rz(-0.37173158) q[1];
sx q[1];
rz(-1.0075684) q[1];
sx q[1];
rz(-2.380611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8801988) q[0];
sx q[0];
rz(-2.7028599) q[0];
sx q[0];
rz(2.3068271) q[0];
rz(1.0521616) q[2];
sx q[2];
rz(-2.385879) q[2];
sx q[2];
rz(-1.99988) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63572143) q[1];
sx q[1];
rz(-2.2835915) q[1];
sx q[1];
rz(-2.9221623) q[1];
rz(-2.2888921) q[3];
sx q[3];
rz(-1.5745244) q[3];
sx q[3];
rz(-2.0223597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.80428213) q[2];
sx q[2];
rz(-2.2411942) q[2];
sx q[2];
rz(2.6864181) q[2];
rz(0.69194397) q[3];
sx q[3];
rz(-0.71440905) q[3];
sx q[3];
rz(2.3327995) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96875018) q[0];
sx q[0];
rz(-1.3469561) q[0];
sx q[0];
rz(-0.2562879) q[0];
rz(-1.434727) q[1];
sx q[1];
rz(-1.4204357) q[1];
sx q[1];
rz(-0.11875471) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8855806) q[0];
sx q[0];
rz(-2.6654454) q[0];
sx q[0];
rz(-0.34389596) q[0];
x q[1];
rz(0.1280814) q[2];
sx q[2];
rz(-1.8616779) q[2];
sx q[2];
rz(-1.1156991) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9164711) q[1];
sx q[1];
rz(-2.3206704) q[1];
sx q[1];
rz(0.31912843) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3975735) q[3];
sx q[3];
rz(-2.3459917) q[3];
sx q[3];
rz(-2.6156486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8541096) q[2];
sx q[2];
rz(-2.867925) q[2];
sx q[2];
rz(-1.0556833) q[2];
rz(-1.4500729) q[3];
sx q[3];
rz(-1.4472716) q[3];
sx q[3];
rz(-0.84901866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.56181041) q[0];
sx q[0];
rz(-2.9636443) q[0];
sx q[0];
rz(1.6135038) q[0];
rz(1.9644507) q[1];
sx q[1];
rz(-1.2700932) q[1];
sx q[1];
rz(1.5783763) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10229853) q[0];
sx q[0];
rz(-1.5015409) q[0];
sx q[0];
rz(-0.62634344) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2424395) q[2];
sx q[2];
rz(-1.17982) q[2];
sx q[2];
rz(-0.41074387) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.745143) q[1];
sx q[1];
rz(-2.1985594) q[1];
sx q[1];
rz(1.3954074) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0061697841) q[3];
sx q[3];
rz(-2.7321762) q[3];
sx q[3];
rz(-2.1751884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.13581181) q[2];
sx q[2];
rz(-1.1567189) q[2];
sx q[2];
rz(1.8394252) q[2];
rz(-0.33995134) q[3];
sx q[3];
rz(-2.3348742) q[3];
sx q[3];
rz(-2.4792041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.133701) q[0];
sx q[0];
rz(-1.5971203) q[0];
sx q[0];
rz(1.9418465) q[0];
rz(1.3613191) q[1];
sx q[1];
rz(-2.168455) q[1];
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
rz(2.0073246) q[0];
rz(-0.33249929) q[2];
sx q[2];
rz(-2.6869171) q[2];
sx q[2];
rz(-2.5408059) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6276045) q[1];
sx q[1];
rz(-1.1179233) q[1];
sx q[1];
rz(-1.4195561) q[1];
rz(-pi) q[2];
rz(1.6541566) q[3];
sx q[3];
rz(-2.8417086) q[3];
sx q[3];
rz(0.56902992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5100688) q[2];
sx q[2];
rz(-0.57454595) q[2];
sx q[2];
rz(2.6681382) q[2];
rz(1.2724642) q[3];
sx q[3];
rz(-1.3926287) q[3];
sx q[3];
rz(2.9191391) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1943787) q[0];
sx q[0];
rz(-0.65776238) q[0];
sx q[0];
rz(-0.33313242) q[0];
rz(-0.22854742) q[1];
sx q[1];
rz(-1.8014149) q[1];
sx q[1];
rz(-0.57565912) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27028337) q[0];
sx q[0];
rz(-0.76437274) q[0];
sx q[0];
rz(1.4174711) q[0];
rz(-pi) q[1];
rz(-2.5533122) q[2];
sx q[2];
rz(-2.5856087) q[2];
sx q[2];
rz(2.5997272) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0927432) q[1];
sx q[1];
rz(-1.7705838) q[1];
sx q[1];
rz(1.9034991) q[1];
x q[2];
rz(0.48124619) q[3];
sx q[3];
rz(-1.3433045) q[3];
sx q[3];
rz(0.15531048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26611844) q[2];
sx q[2];
rz(-0.55142752) q[2];
sx q[2];
rz(1.4761338) q[2];
rz(-1.1009781) q[3];
sx q[3];
rz(-1.063238) q[3];
sx q[3];
rz(-0.5180009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70676409) q[0];
sx q[0];
rz(-1.6521104) q[0];
sx q[0];
rz(-1.5887194) q[0];
x q[1];
rz(2.7935289) q[2];
sx q[2];
rz(-1.3937147) q[2];
sx q[2];
rz(2.7325316) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0982837) q[1];
sx q[1];
rz(-1.829147) q[1];
sx q[1];
rz(-1.1275379) q[1];
rz(1.7264113) q[3];
sx q[3];
rz(-2.3978516) q[3];
sx q[3];
rz(-1.9187601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8660628) q[2];
sx q[2];
rz(-0.58774647) q[2];
sx q[2];
rz(2.8928939) q[2];
rz(1.9281049) q[3];
sx q[3];
rz(-1.6971734) q[3];
sx q[3];
rz(-3.0799589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1520749) q[0];
sx q[0];
rz(-2.2353421) q[0];
sx q[0];
rz(2.4549947) q[0];
rz(1.1129414) q[1];
sx q[1];
rz(-2.3390892) q[1];
sx q[1];
rz(-0.31563219) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7677143) q[0];
sx q[0];
rz(-1.5247824) q[0];
sx q[0];
rz(1.2476646) q[0];
rz(0.75982538) q[2];
sx q[2];
rz(-1.3889599) q[2];
sx q[2];
rz(2.1450617) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5110398) q[1];
sx q[1];
rz(-0.91088089) q[1];
sx q[1];
rz(3.0223165) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0635438) q[3];
sx q[3];
rz(-2.4025318) q[3];
sx q[3];
rz(-2.259425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.687872) q[2];
sx q[2];
rz(-0.56679711) q[2];
sx q[2];
rz(-0.68457121) q[2];
rz(-1.941393) q[3];
sx q[3];
rz(-2.0122416) q[3];
sx q[3];
rz(-0.17957345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0803364) q[0];
sx q[0];
rz(-2.6604524) q[0];
sx q[0];
rz(-0.0048333724) q[0];
rz(-1.5176516) q[1];
sx q[1];
rz(-2.3976517) q[1];
sx q[1];
rz(2.8878816) q[1];
x q[2];
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
rz(-2.6275835) q[2];
sx q[2];
rz(-2.3367662) q[2];
sx q[2];
rz(-2.6857306) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.0066792329) q[1];
sx q[1];
rz(-0.90613922) q[1];
sx q[1];
rz(0.24941872) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5911035) q[3];
sx q[3];
rz(-1.4006249) q[3];
sx q[3];
rz(-0.44743928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1634875) q[2];
sx q[2];
rz(-1.2349962) q[2];
sx q[2];
rz(1.9592436) q[2];
rz(-2.4649418) q[3];
sx q[3];
rz(-2.7550321) q[3];
sx q[3];
rz(-2.1277908) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5273298) q[0];
sx q[0];
rz(-0.78582055) q[0];
sx q[0];
rz(0.2926122) q[0];
rz(1.0489427) q[1];
sx q[1];
rz(-1.4194149) q[1];
sx q[1];
rz(0.70379757) q[1];
rz(1.232515) q[2];
sx q[2];
rz(-1.1223553) q[2];
sx q[2];
rz(-1.9434402) q[2];
rz(3.0426619) q[3];
sx q[3];
rz(-1.9361648) q[3];
sx q[3];
rz(1.786644) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
