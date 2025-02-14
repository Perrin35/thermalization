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
rz(-0.55627745) q[0];
sx q[0];
rz(3.102432) q[0];
sx q[0];
rz(10.089212) q[0];
rz(-0.18222624) q[1];
sx q[1];
rz(4.8290904) q[1];
sx q[1];
rz(11.157142) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0213703) q[0];
sx q[0];
rz(-0.5628559) q[0];
sx q[0];
rz(-0.31546202) q[0];
x q[1];
rz(1.9033236) q[2];
sx q[2];
rz(-1.1477787) q[2];
sx q[2];
rz(0.11655434) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4709027) q[1];
sx q[1];
rz(-1.8876612) q[1];
sx q[1];
rz(2.5554727) q[1];
rz(-1.9326747) q[3];
sx q[3];
rz(-1.3986949) q[3];
sx q[3];
rz(0.73735305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0613784) q[2];
sx q[2];
rz(-2.5900216) q[2];
sx q[2];
rz(-2.3870094) q[2];
rz(1.6825698) q[3];
sx q[3];
rz(-0.85586923) q[3];
sx q[3];
rz(-1.7350908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.215312) q[0];
sx q[0];
rz(-0.54406852) q[0];
sx q[0];
rz(-1.9790443) q[0];
rz(-0.7112208) q[1];
sx q[1];
rz(-2.4237207) q[1];
sx q[1];
rz(-0.1444764) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1294589) q[0];
sx q[0];
rz(-2.8367865) q[0];
sx q[0];
rz(0.27979677) q[0];
x q[1];
rz(0.87717339) q[2];
sx q[2];
rz(-2.2899004) q[2];
sx q[2];
rz(2.4404877) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4481758) q[1];
sx q[1];
rz(-0.87478335) q[1];
sx q[1];
rz(-2.059518) q[1];
x q[2];
rz(-1.9599592) q[3];
sx q[3];
rz(-2.4414821) q[3];
sx q[3];
rz(-0.35193974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43313906) q[2];
sx q[2];
rz(-1.3044367) q[2];
sx q[2];
rz(-0.30678314) q[2];
rz(-0.77543801) q[3];
sx q[3];
rz(-2.756835) q[3];
sx q[3];
rz(-2.9966089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48258346) q[0];
sx q[0];
rz(-0.44335303) q[0];
sx q[0];
rz(0.79202598) q[0];
rz(-1.5576942) q[1];
sx q[1];
rz(-1.3521103) q[1];
sx q[1];
rz(-2.1477594) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3344468) q[0];
sx q[0];
rz(-0.69753555) q[0];
sx q[0];
rz(-0.56484449) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5517919) q[2];
sx q[2];
rz(-1.4695393) q[2];
sx q[2];
rz(2.4112005) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5372588) q[1];
sx q[1];
rz(-0.52642614) q[1];
sx q[1];
rz(-2.1597305) q[1];
rz(-pi) q[2];
rz(-1.1937856) q[3];
sx q[3];
rz(-2.6600231) q[3];
sx q[3];
rz(1.6653614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1899679) q[2];
sx q[2];
rz(-0.8364532) q[2];
sx q[2];
rz(-2.9659029) q[2];
rz(-0.22819337) q[3];
sx q[3];
rz(-0.56932813) q[3];
sx q[3];
rz(-0.5051676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2100385) q[0];
sx q[0];
rz(-2.2966972) q[0];
sx q[0];
rz(-1.5799874) q[0];
rz(-2.2700894) q[1];
sx q[1];
rz(-1.5305287) q[1];
sx q[1];
rz(-0.16564381) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79194389) q[0];
sx q[0];
rz(-1.569618) q[0];
sx q[0];
rz(-3.1399957) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5524988) q[2];
sx q[2];
rz(-2.6958272) q[2];
sx q[2];
rz(-0.95761567) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9479182) q[1];
sx q[1];
rz(-2.3478824) q[1];
sx q[1];
rz(3.0315184) q[1];
x q[2];
rz(-0.58067643) q[3];
sx q[3];
rz(-1.3370974) q[3];
sx q[3];
rz(-3.0729938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6988301) q[2];
sx q[2];
rz(-0.58219588) q[2];
sx q[2];
rz(-0.77787918) q[2];
rz(-0.88464087) q[3];
sx q[3];
rz(-1.9886465) q[3];
sx q[3];
rz(-0.9790498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1409461) q[0];
sx q[0];
rz(-1.0362754) q[0];
sx q[0];
rz(2.272814) q[0];
rz(-2.2390305) q[1];
sx q[1];
rz(-1.9156009) q[1];
sx q[1];
rz(-0.57194078) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7392699) q[0];
sx q[0];
rz(-2.2796541) q[0];
sx q[0];
rz(0.92415442) q[0];
rz(-pi) q[1];
rz(-0.50359132) q[2];
sx q[2];
rz(-1.9707754) q[2];
sx q[2];
rz(-1.3445889) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1871698) q[1];
sx q[1];
rz(-2.5989418) q[1];
sx q[1];
rz(-0.95466787) q[1];
x q[2];
rz(-2.8945893) q[3];
sx q[3];
rz(-1.9387783) q[3];
sx q[3];
rz(2.0968269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9945485) q[2];
sx q[2];
rz(-0.93783718) q[2];
sx q[2];
rz(0.63329548) q[2];
rz(2.5607732) q[3];
sx q[3];
rz(-1.7827026) q[3];
sx q[3];
rz(2.4766428) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13034114) q[0];
sx q[0];
rz(-2.7524152) q[0];
sx q[0];
rz(-1.8036386) q[0];
rz(-1.7300216) q[1];
sx q[1];
rz(-2.7346225) q[1];
sx q[1];
rz(0.79708797) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4914843) q[0];
sx q[0];
rz(-2.8239125) q[0];
sx q[0];
rz(-1.0450715) q[0];
x q[1];
rz(2.1363968) q[2];
sx q[2];
rz(-1.505064) q[2];
sx q[2];
rz(-2.1594723) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2458249) q[1];
sx q[1];
rz(-1.3029769) q[1];
sx q[1];
rz(2.9755249) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8800297) q[3];
sx q[3];
rz(-0.54094523) q[3];
sx q[3];
rz(0.23829392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.36593124) q[2];
sx q[2];
rz(-2.4727827) q[2];
sx q[2];
rz(-2.456341) q[2];
rz(0.25964409) q[3];
sx q[3];
rz(-0.97340596) q[3];
sx q[3];
rz(1.9297622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24945666) q[0];
sx q[0];
rz(-0.15661713) q[0];
sx q[0];
rz(-2.6020965) q[0];
rz(-0.19142137) q[1];
sx q[1];
rz(-1.4861264) q[1];
sx q[1];
rz(-2.6991381) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0809796) q[0];
sx q[0];
rz(-0.0091305841) q[0];
sx q[0];
rz(1.5332444) q[0];
rz(-pi) q[1];
rz(2.7436951) q[2];
sx q[2];
rz(-0.98491633) q[2];
sx q[2];
rz(2.2316124) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6574508) q[1];
sx q[1];
rz(-1.9062348) q[1];
sx q[1];
rz(2.2141333) q[1];
rz(2.3549162) q[3];
sx q[3];
rz(-1.0287675) q[3];
sx q[3];
rz(1.4876798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29360867) q[2];
sx q[2];
rz(-0.30985761) q[2];
sx q[2];
rz(2.4388745) q[2];
rz(2.8637049) q[3];
sx q[3];
rz(-1.0669471) q[3];
sx q[3];
rz(-1.3387298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24669312) q[0];
sx q[0];
rz(-0.86240697) q[0];
sx q[0];
rz(2.605751) q[0];
rz(2.4193173) q[1];
sx q[1];
rz(-1.4403789) q[1];
sx q[1];
rz(2.3586418) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0412126) q[0];
sx q[0];
rz(-1.4836652) q[0];
sx q[0];
rz(-1.9308596) q[0];
x q[1];
rz(-2.7442619) q[2];
sx q[2];
rz(-1.6226007) q[2];
sx q[2];
rz(-1.0578294) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6440575) q[1];
sx q[1];
rz(-2.3573207) q[1];
sx q[1];
rz(-2.3375744) q[1];
rz(2.1184741) q[3];
sx q[3];
rz(-1.191282) q[3];
sx q[3];
rz(1.4364157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4202412) q[2];
sx q[2];
rz(-0.46478096) q[2];
sx q[2];
rz(-2.6381524) q[2];
rz(0.33468801) q[3];
sx q[3];
rz(-1.3379593) q[3];
sx q[3];
rz(2.9065342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6675785) q[0];
sx q[0];
rz(-2.773371) q[0];
sx q[0];
rz(2.0817122) q[0];
rz(-0.31479752) q[1];
sx q[1];
rz(-1.7960725) q[1];
sx q[1];
rz(1.559929) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81875694) q[0];
sx q[0];
rz(-2.4222932) q[0];
sx q[0];
rz(-0.34157217) q[0];
rz(-pi) q[1];
rz(-1.306172) q[2];
sx q[2];
rz(-1.5381863) q[2];
sx q[2];
rz(0.6219686) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.25645721) q[1];
sx q[1];
rz(-2.4607435) q[1];
sx q[1];
rz(1.9155986) q[1];
rz(-0.022059343) q[3];
sx q[3];
rz(-1.1556871) q[3];
sx q[3];
rz(-0.65655113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9554837) q[2];
sx q[2];
rz(-1.1649705) q[2];
sx q[2];
rz(-1.0355518) q[2];
rz(-1.1459076) q[3];
sx q[3];
rz(-2.0290387) q[3];
sx q[3];
rz(0.15315332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2027407) q[0];
sx q[0];
rz(-0.11883141) q[0];
sx q[0];
rz(0.76549292) q[0];
rz(-1.0368404) q[1];
sx q[1];
rz(-1.2988657) q[1];
sx q[1];
rz(0.11229215) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13829198) q[0];
sx q[0];
rz(-1.1213741) q[0];
sx q[0];
rz(1.9339191) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1767581) q[2];
sx q[2];
rz(-0.8886742) q[2];
sx q[2];
rz(-1.8471931) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4203305) q[1];
sx q[1];
rz(-0.76618505) q[1];
sx q[1];
rz(-1.2797194) q[1];
x q[2];
rz(0.78386098) q[3];
sx q[3];
rz(-1.5376159) q[3];
sx q[3];
rz(2.9110094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3969193) q[2];
sx q[2];
rz(-1.2770709) q[2];
sx q[2];
rz(-3.0199158) q[2];
rz(2.7820898) q[3];
sx q[3];
rz(-0.72327852) q[3];
sx q[3];
rz(0.017688964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.551238) q[0];
sx q[0];
rz(-1.6726765) q[0];
sx q[0];
rz(-1.8400675) q[0];
rz(1.7474668) q[1];
sx q[1];
rz(-1.9448517) q[1];
sx q[1];
rz(1.3762884) q[1];
rz(2.0040705) q[2];
sx q[2];
rz(-2.1480297) q[2];
sx q[2];
rz(0.99110023) q[2];
rz(1.2944503) q[3];
sx q[3];
rz(-2.69097) q[3];
sx q[3];
rz(0.16926058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
