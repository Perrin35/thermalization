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
rz(-1.902154) q[0];
sx q[0];
rz(-1.3286989) q[0];
sx q[0];
rz(0.21188307) q[0];
rz(-1.4172685) q[1];
sx q[1];
rz(-2.6091726) q[1];
sx q[1];
rz(0.37791696) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1013452) q[0];
sx q[0];
rz(-1.5703526) q[0];
sx q[0];
rz(-1.5803278) q[0];
rz(-pi) q[1];
rz(-0.67899668) q[2];
sx q[2];
rz(-1.8583991) q[2];
sx q[2];
rz(1.2305413) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7211044) q[1];
sx q[1];
rz(-1.2961565) q[1];
sx q[1];
rz(2.4824449) q[1];
rz(2.3232949) q[3];
sx q[3];
rz(-0.94486559) q[3];
sx q[3];
rz(-1.4656386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8523031) q[2];
sx q[2];
rz(-2.0630344) q[2];
sx q[2];
rz(-3.0618073) q[2];
rz(-2.1763109) q[3];
sx q[3];
rz(-1.691317) q[3];
sx q[3];
rz(-2.4108346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47736436) q[0];
sx q[0];
rz(-1.0034765) q[0];
sx q[0];
rz(-0.088951237) q[0];
rz(1.7695919) q[1];
sx q[1];
rz(-1.4274495) q[1];
sx q[1];
rz(-2.037183) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8744132) q[0];
sx q[0];
rz(-0.66775087) q[0];
sx q[0];
rz(0.37607583) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78271336) q[2];
sx q[2];
rz(-1.209895) q[2];
sx q[2];
rz(2.1825143) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1762878) q[1];
sx q[1];
rz(-1.4215648) q[1];
sx q[1];
rz(2.01841) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4711963) q[3];
sx q[3];
rz(-2.1079113) q[3];
sx q[3];
rz(2.296534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8770807) q[2];
sx q[2];
rz(-1.1254213) q[2];
sx q[2];
rz(1.6748927) q[2];
rz(0.029021164) q[3];
sx q[3];
rz(-2.1252188) q[3];
sx q[3];
rz(-2.4513054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2324209) q[0];
sx q[0];
rz(-3.0747774) q[0];
sx q[0];
rz(0.27798852) q[0];
rz(1.4311283) q[1];
sx q[1];
rz(-2.2093096) q[1];
sx q[1];
rz(-2.7412282) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4543661) q[0];
sx q[0];
rz(-2.3355744) q[0];
sx q[0];
rz(-2.5044652) q[0];
x q[1];
rz(0.68636559) q[2];
sx q[2];
rz(-1.7698506) q[2];
sx q[2];
rz(-2.484427) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6975721) q[1];
sx q[1];
rz(-2.6409147) q[1];
sx q[1];
rz(-0.33659192) q[1];
rz(-0.34137643) q[3];
sx q[3];
rz(-1.2814643) q[3];
sx q[3];
rz(3.0738664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4572738) q[2];
sx q[2];
rz(-1.6088586) q[2];
sx q[2];
rz(0.45480967) q[2];
rz(-1.8999892) q[3];
sx q[3];
rz(-2.1865215) q[3];
sx q[3];
rz(-0.72030592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1860745) q[0];
sx q[0];
rz(-1.3294514) q[0];
sx q[0];
rz(-0.00051001471) q[0];
rz(-2.5406802) q[1];
sx q[1];
rz(-0.83952236) q[1];
sx q[1];
rz(2.9972163) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3568319) q[0];
sx q[0];
rz(-0.5934754) q[0];
sx q[0];
rz(-1.1952101) q[0];
x q[1];
rz(0.52835502) q[2];
sx q[2];
rz(-2.8502185) q[2];
sx q[2];
rz(1.0519415) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51203242) q[1];
sx q[1];
rz(-1.7981537) q[1];
sx q[1];
rz(0.65452202) q[1];
rz(-1.6779283) q[3];
sx q[3];
rz(-1.0262353) q[3];
sx q[3];
rz(0.0060826172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1986177) q[2];
sx q[2];
rz(-1.696442) q[2];
sx q[2];
rz(0.96251881) q[2];
rz(1.5396384) q[3];
sx q[3];
rz(-1.736015) q[3];
sx q[3];
rz(0.34077728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9050423) q[0];
sx q[0];
rz(-1.6860697) q[0];
sx q[0];
rz(1.9566253) q[0];
rz(0.22625893) q[1];
sx q[1];
rz(-2.2626651) q[1];
sx q[1];
rz(-1.0386946) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13030355) q[0];
sx q[0];
rz(-1.652392) q[0];
sx q[0];
rz(-0.58873436) q[0];
rz(-pi) q[1];
rz(2.4034385) q[2];
sx q[2];
rz(-1.0760436) q[2];
sx q[2];
rz(-0.02502266) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9925527) q[1];
sx q[1];
rz(-0.91616183) q[1];
sx q[1];
rz(0.66019411) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1306954) q[3];
sx q[3];
rz(-1.6301042) q[3];
sx q[3];
rz(-1.1820933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.431939) q[2];
sx q[2];
rz(-2.4434872) q[2];
sx q[2];
rz(-2.8254438) q[2];
rz(-1.6759253) q[3];
sx q[3];
rz(-2.0114055) q[3];
sx q[3];
rz(1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6370711) q[0];
sx q[0];
rz(-0.9032473) q[0];
sx q[0];
rz(1.1035408) q[0];
rz(-0.48031131) q[1];
sx q[1];
rz(-2.4791398) q[1];
sx q[1];
rz(-2.5441817) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7561853) q[0];
sx q[0];
rz(-2.8592337) q[0];
sx q[0];
rz(2.307111) q[0];
x q[1];
rz(0.2603724) q[2];
sx q[2];
rz(-2.3790759) q[2];
sx q[2];
rz(2.0338361) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.061042) q[1];
sx q[1];
rz(-1.4054728) q[1];
sx q[1];
rz(-0.25456659) q[1];
rz(1.8036929) q[3];
sx q[3];
rz(-0.38023708) q[3];
sx q[3];
rz(-1.5986795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9395113) q[2];
sx q[2];
rz(-1.5152405) q[2];
sx q[2];
rz(-2.0651979) q[2];
rz(1.1427897) q[3];
sx q[3];
rz(-0.79311526) q[3];
sx q[3];
rz(2.6514261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6550605) q[0];
sx q[0];
rz(-1.158411) q[0];
sx q[0];
rz(3.1031188) q[0];
rz(3.0768652) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(0.23385349) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7215639) q[0];
sx q[0];
rz(-1.7917243) q[0];
sx q[0];
rz(-1.7608282) q[0];
x q[1];
rz(2.269417) q[2];
sx q[2];
rz(-0.46469122) q[2];
sx q[2];
rz(1.2810436) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0044627) q[1];
sx q[1];
rz(-2.4871768) q[1];
sx q[1];
rz(-0.14738247) q[1];
rz(-pi) q[2];
rz(2.7247808) q[3];
sx q[3];
rz(-2.9997065) q[3];
sx q[3];
rz(2.8519423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6528299) q[2];
sx q[2];
rz(-1.5831999) q[2];
sx q[2];
rz(-2.5433507) q[2];
rz(-0.11387842) q[3];
sx q[3];
rz(-1.3860393) q[3];
sx q[3];
rz(0.84754506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4050196) q[0];
sx q[0];
rz(-0.84252715) q[0];
sx q[0];
rz(-0.2463499) q[0];
rz(1.2975289) q[1];
sx q[1];
rz(-1.1187226) q[1];
sx q[1];
rz(2.6447703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5733872) q[0];
sx q[0];
rz(-0.77952379) q[0];
sx q[0];
rz(2.5318092) q[0];
rz(2.4364528) q[2];
sx q[2];
rz(-1.6821096) q[2];
sx q[2];
rz(-1.4347347) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0441974) q[1];
sx q[1];
rz(-0.5438416) q[1];
sx q[1];
rz(1.103765) q[1];
x q[2];
rz(1.0786177) q[3];
sx q[3];
rz(-1.4340622) q[3];
sx q[3];
rz(1.1059974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7184427) q[2];
sx q[2];
rz(-2.6079874) q[2];
sx q[2];
rz(-1.144484) q[2];
rz(-1.9874969) q[3];
sx q[3];
rz(-1.4242947) q[3];
sx q[3];
rz(0.19788876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.4474739) q[0];
sx q[0];
rz(-2.9070774) q[0];
sx q[0];
rz(-0.06614729) q[0];
rz(-1.1478395) q[1];
sx q[1];
rz(-1.3775974) q[1];
sx q[1];
rz(-2.537421) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3124638) q[0];
sx q[0];
rz(-1.365571) q[0];
sx q[0];
rz(-2.8602296) q[0];
rz(-1.4140903) q[2];
sx q[2];
rz(-2.092166) q[2];
sx q[2];
rz(-0.40796134) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7614038) q[1];
sx q[1];
rz(-2.1288925) q[1];
sx q[1];
rz(1.7978884) q[1];
rz(-pi) q[2];
rz(-0.82138942) q[3];
sx q[3];
rz(-2.2940972) q[3];
sx q[3];
rz(-0.57904527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26312795) q[2];
sx q[2];
rz(-2.2915514) q[2];
sx q[2];
rz(-2.7086332) q[2];
rz(1.912502) q[3];
sx q[3];
rz(-1.9357598) q[3];
sx q[3];
rz(-1.3273201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1935254) q[0];
sx q[0];
rz(-2.0663517) q[0];
sx q[0];
rz(1.8684335) q[0];
rz(-2.676447) q[1];
sx q[1];
rz(-1.3659313) q[1];
sx q[1];
rz(0.21496162) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8037655) q[0];
sx q[0];
rz(-1.318232) q[0];
sx q[0];
rz(-2.4830677) q[0];
x q[1];
rz(0.6647756) q[2];
sx q[2];
rz(-0.6419581) q[2];
sx q[2];
rz(0.73981111) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3167606) q[1];
sx q[1];
rz(-2.2344898) q[1];
sx q[1];
rz(0.65726991) q[1];
rz(-0.56193476) q[3];
sx q[3];
rz(-1.3519545) q[3];
sx q[3];
rz(2.9666025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8421858) q[2];
sx q[2];
rz(-0.81373787) q[2];
sx q[2];
rz(-2.1827533) q[2];
rz(1.1622608) q[3];
sx q[3];
rz(-1.6543038) q[3];
sx q[3];
rz(1.6963262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5398298) q[0];
sx q[0];
rz(-1.6015263) q[0];
sx q[0];
rz(1.482561) q[0];
rz(2.4330347) q[1];
sx q[1];
rz(-0.19150145) q[1];
sx q[1];
rz(2.3932744) q[1];
rz(-2.5790527) q[2];
sx q[2];
rz(-1.8481135) q[2];
sx q[2];
rz(-2.6279966) q[2];
rz(-2.9293168) q[3];
sx q[3];
rz(-0.73321453) q[3];
sx q[3];
rz(-1.1238255) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
