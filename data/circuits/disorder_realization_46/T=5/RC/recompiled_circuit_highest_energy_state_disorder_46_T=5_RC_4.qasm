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
rz(1.8019851) q[0];
sx q[0];
rz(-2.7924502) q[0];
sx q[0];
rz(-2.1139297) q[0];
rz(3.1066306) q[1];
sx q[1];
rz(-2.1244013) q[1];
sx q[1];
rz(-1.6962167) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6562614) q[0];
sx q[0];
rz(-2.4527271) q[0];
sx q[0];
rz(0.553225) q[0];
x q[1];
rz(0.17757444) q[2];
sx q[2];
rz(-1.53731) q[2];
sx q[2];
rz(0.78508961) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.45856341) q[1];
sx q[1];
rz(-1.7818319) q[1];
sx q[1];
rz(-2.8082735) q[1];
x q[2];
rz(-0.39528109) q[3];
sx q[3];
rz(-0.78215212) q[3];
sx q[3];
rz(0.30667337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8029636) q[2];
sx q[2];
rz(-2.6821319) q[2];
sx q[2];
rz(-2.6098693) q[2];
rz(-1.7470597) q[3];
sx q[3];
rz(-1.2732384) q[3];
sx q[3];
rz(1.1751706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26674536) q[0];
sx q[0];
rz(-1.5044455) q[0];
sx q[0];
rz(-0.7919842) q[0];
rz(2.2199471) q[1];
sx q[1];
rz(-1.6519203) q[1];
sx q[1];
rz(2.0335061) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9757248) q[0];
sx q[0];
rz(-2.168403) q[0];
sx q[0];
rz(-3.0063666) q[0];
rz(-pi) q[1];
rz(-0.048594995) q[2];
sx q[2];
rz(-0.94779769) q[2];
sx q[2];
rz(-2.5041265) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.82906065) q[1];
sx q[1];
rz(-0.8735114) q[1];
sx q[1];
rz(-2.6379449) q[1];
rz(-pi) q[2];
rz(1.2543801) q[3];
sx q[3];
rz(-1.5174447) q[3];
sx q[3];
rz(-0.90094588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7532928) q[2];
sx q[2];
rz(-1.029195) q[2];
sx q[2];
rz(1.252906) q[2];
rz(-2.5168354) q[3];
sx q[3];
rz(-1.2478991) q[3];
sx q[3];
rz(-2.6784082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6817634) q[0];
sx q[0];
rz(-0.1777996) q[0];
sx q[0];
rz(-2.8681927) q[0];
rz(-0.071648486) q[1];
sx q[1];
rz(-2.007808) q[1];
sx q[1];
rz(-1.515548) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0833383) q[0];
sx q[0];
rz(-1.951955) q[0];
sx q[0];
rz(-0.55321557) q[0];
rz(0.46892445) q[2];
sx q[2];
rz(-1.8919626) q[2];
sx q[2];
rz(1.657148) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7715986) q[1];
sx q[1];
rz(-2.150476) q[1];
sx q[1];
rz(-1.7607776) q[1];
rz(0.97730009) q[3];
sx q[3];
rz(-1.7827099) q[3];
sx q[3];
rz(0.43928543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.36149725) q[2];
sx q[2];
rz(-0.7146892) q[2];
sx q[2];
rz(1.3698461) q[2];
rz(1.0684377) q[3];
sx q[3];
rz(-1.3911824) q[3];
sx q[3];
rz(2.1865602) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84900981) q[0];
sx q[0];
rz(-0.42790258) q[0];
sx q[0];
rz(-0.12246116) q[0];
rz(0.017008688) q[1];
sx q[1];
rz(-1.6288792) q[1];
sx q[1];
rz(-0.88517991) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95851446) q[0];
sx q[0];
rz(-1.4349271) q[0];
sx q[0];
rz(-2.0625173) q[0];
rz(-pi) q[1];
rz(-2.3058168) q[2];
sx q[2];
rz(-1.3124497) q[2];
sx q[2];
rz(-1.6222184) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.0065476848) q[1];
sx q[1];
rz(-1.1741877) q[1];
sx q[1];
rz(1.0746677) q[1];
rz(-2.9131575) q[3];
sx q[3];
rz(-0.71492787) q[3];
sx q[3];
rz(-2.950516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.433832) q[2];
sx q[2];
rz(-1.8261352) q[2];
sx q[2];
rz(3.0641277) q[2];
rz(-0.42099434) q[3];
sx q[3];
rz(-1.1869895) q[3];
sx q[3];
rz(2.8903956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0780599) q[0];
sx q[0];
rz(-0.01883004) q[0];
sx q[0];
rz(-0.68191648) q[0];
rz(0.49184999) q[1];
sx q[1];
rz(-2.1147155) q[1];
sx q[1];
rz(-1.9815725) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9710142) q[0];
sx q[0];
rz(-1.7587629) q[0];
sx q[0];
rz(2.2925799) q[0];
rz(-0.57038681) q[2];
sx q[2];
rz(-1.5231992) q[2];
sx q[2];
rz(0.80709761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1067137) q[1];
sx q[1];
rz(-0.95297652) q[1];
sx q[1];
rz(2.930083) q[1];
rz(1.9914636) q[3];
sx q[3];
rz(-1.3004817) q[3];
sx q[3];
rz(-0.079416954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.68449768) q[2];
sx q[2];
rz(-0.86898494) q[2];
sx q[2];
rz(2.1333466) q[2];
rz(1.5021987) q[3];
sx q[3];
rz(-1.1049756) q[3];
sx q[3];
rz(-0.26386279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0428001) q[0];
sx q[0];
rz(-1.5473939) q[0];
sx q[0];
rz(2.7368326) q[0];
rz(-0.81819355) q[1];
sx q[1];
rz(-1.8408366) q[1];
sx q[1];
rz(0.71664804) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7966864) q[0];
sx q[0];
rz(-0.27936882) q[0];
sx q[0];
rz(3.1059103) q[0];
rz(0.76397447) q[2];
sx q[2];
rz(-1.171955) q[2];
sx q[2];
rz(0.13834902) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4617052) q[1];
sx q[1];
rz(-0.53820921) q[1];
sx q[1];
rz(-2.6160282) q[1];
rz(-0.53933177) q[3];
sx q[3];
rz(-1.8908653) q[3];
sx q[3];
rz(1.1160929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8019668) q[2];
sx q[2];
rz(-0.64029396) q[2];
sx q[2];
rz(2.4862945) q[2];
rz(-2.7845434) q[3];
sx q[3];
rz(-1.3909631) q[3];
sx q[3];
rz(-2.6220139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94012564) q[0];
sx q[0];
rz(-0.31620142) q[0];
sx q[0];
rz(2.1215718) q[0];
rz(2.5992498) q[1];
sx q[1];
rz(-2.0261363) q[1];
sx q[1];
rz(1.382359) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4931204) q[0];
sx q[0];
rz(-2.1628404) q[0];
sx q[0];
rz(0.73541321) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5711201) q[2];
sx q[2];
rz(-2.8920724) q[2];
sx q[2];
rz(-0.0970627) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89703694) q[1];
sx q[1];
rz(-0.4677597) q[1];
sx q[1];
rz(1.3964584) q[1];
rz(-pi) q[2];
x q[2];
rz(0.027023836) q[3];
sx q[3];
rz(-0.083649717) q[3];
sx q[3];
rz(0.25318709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0688087) q[2];
sx q[2];
rz(-1.5280318) q[2];
sx q[2];
rz(2.7911348) q[2];
rz(0.46640486) q[3];
sx q[3];
rz(-0.91752183) q[3];
sx q[3];
rz(0.94356999) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7826409) q[0];
sx q[0];
rz(-2.751001) q[0];
sx q[0];
rz(0.95640957) q[0];
rz(1.6920754) q[1];
sx q[1];
rz(-0.78069514) q[1];
sx q[1];
rz(2.4704959) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9195477) q[0];
sx q[0];
rz(-1.4853108) q[0];
sx q[0];
rz(-1.7128904) q[0];
x q[1];
rz(0.568957) q[2];
sx q[2];
rz(-1.3929741) q[2];
sx q[2];
rz(-2.5513955) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2626604) q[1];
sx q[1];
rz(-1.2178246) q[1];
sx q[1];
rz(-0.067990677) q[1];
x q[2];
rz(-2.4375979) q[3];
sx q[3];
rz(-1.3823095) q[3];
sx q[3];
rz(2.2437167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.87390071) q[2];
sx q[2];
rz(-0.42517069) q[2];
sx q[2];
rz(-1.1262061) q[2];
rz(0.33024427) q[3];
sx q[3];
rz(-1.4010022) q[3];
sx q[3];
rz(-0.11708524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0432334) q[0];
sx q[0];
rz(-2.5625304) q[0];
sx q[0];
rz(-0.057057127) q[0];
rz(0.13380274) q[1];
sx q[1];
rz(-2.0963142) q[1];
sx q[1];
rz(0.71279508) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063449115) q[0];
sx q[0];
rz(-0.3219372) q[0];
sx q[0];
rz(0.19970317) q[0];
x q[1];
rz(0.78929094) q[2];
sx q[2];
rz(-2.5954491) q[2];
sx q[2];
rz(-0.30420732) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.61049938) q[1];
sx q[1];
rz(-1.1720017) q[1];
sx q[1];
rz(-2.628875) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2650293) q[3];
sx q[3];
rz(-1.5086155) q[3];
sx q[3];
rz(-2.2060195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5130676) q[2];
sx q[2];
rz(-1.2863938) q[2];
sx q[2];
rz(0.081550278) q[2];
rz(-1.2201803) q[3];
sx q[3];
rz(-2.8704075) q[3];
sx q[3];
rz(2.0512106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.14209014) q[0];
sx q[0];
rz(-1.6195848) q[0];
sx q[0];
rz(1.5600486) q[0];
rz(1.1355431) q[1];
sx q[1];
rz(-1.138843) q[1];
sx q[1];
rz(1.9948237) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97871507) q[0];
sx q[0];
rz(-2.204127) q[0];
sx q[0];
rz(-2.1970554) q[0];
rz(-0.76374526) q[2];
sx q[2];
rz(-0.95852214) q[2];
sx q[2];
rz(-2.6291922) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1022328) q[1];
sx q[1];
rz(-0.22645849) q[1];
sx q[1];
rz(0.38180085) q[1];
rz(-pi) q[2];
rz(-1.7079748) q[3];
sx q[3];
rz(-0.56932031) q[3];
sx q[3];
rz(-0.63976054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5494988) q[2];
sx q[2];
rz(-1.8213976) q[2];
sx q[2];
rz(-2.7756694) q[2];
rz(-2.7538815) q[3];
sx q[3];
rz(-2.0230899) q[3];
sx q[3];
rz(-1.4661192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53030071) q[0];
sx q[0];
rz(-2.3983751) q[0];
sx q[0];
rz(-0.30839738) q[0];
rz(-0.74956924) q[1];
sx q[1];
rz(-1.1309962) q[1];
sx q[1];
rz(1.8684734) q[1];
rz(-1.1403492) q[2];
sx q[2];
rz(-2.4348209) q[2];
sx q[2];
rz(-1.5904279) q[2];
rz(1.9182792) q[3];
sx q[3];
rz(-0.96746222) q[3];
sx q[3];
rz(-0.37012561) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
