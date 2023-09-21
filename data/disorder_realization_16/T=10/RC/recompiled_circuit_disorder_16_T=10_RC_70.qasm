OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6668532) q[0];
sx q[0];
rz(-2.3119976) q[0];
sx q[0];
rz(-0.15396804) q[0];
rz(0.83377588) q[1];
sx q[1];
rz(4.1339388) q[1];
sx q[1];
rz(9.0864656) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96903893) q[0];
sx q[0];
rz(-1.2533029) q[0];
sx q[0];
rz(2.88455) q[0];
rz(2.2489684) q[2];
sx q[2];
rz(-1.8632338) q[2];
sx q[2];
rz(-0.48722789) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4664647) q[1];
sx q[1];
rz(-2.8505241) q[1];
sx q[1];
rz(-0.18578271) q[1];
x q[2];
rz(1.5428513) q[3];
sx q[3];
rz(-2.6290647) q[3];
sx q[3];
rz(-0.41748369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9989495) q[2];
sx q[2];
rz(-0.34037408) q[2];
sx q[2];
rz(-1.9677229) q[2];
rz(3.0657892) q[3];
sx q[3];
rz(-1.1444164) q[3];
sx q[3];
rz(-0.092806667) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3409815) q[0];
sx q[0];
rz(-1.0656463) q[0];
sx q[0];
rz(3.0766292) q[0];
rz(-2.5669572) q[1];
sx q[1];
rz(-0.42962933) q[1];
sx q[1];
rz(-1.8992791) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64105469) q[0];
sx q[0];
rz(-2.8613052) q[0];
sx q[0];
rz(0.48135249) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9827036) q[2];
sx q[2];
rz(-2.1974568) q[2];
sx q[2];
rz(-1.6275258) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9638483) q[1];
sx q[1];
rz(-1.4114393) q[1];
sx q[1];
rz(-2.6049155) q[1];
rz(-pi) q[2];
rz(-2.5870856) q[3];
sx q[3];
rz(-2.3508361) q[3];
sx q[3];
rz(-2.8046372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.80766455) q[2];
sx q[2];
rz(-1.0753205) q[2];
sx q[2];
rz(0.58369613) q[2];
rz(-2.5675473) q[3];
sx q[3];
rz(-1.1255001) q[3];
sx q[3];
rz(3.0103502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42049256) q[0];
sx q[0];
rz(-2.2476966) q[0];
sx q[0];
rz(2.4131391) q[0];
rz(-1.6473673) q[1];
sx q[1];
rz(-2.7431226) q[1];
sx q[1];
rz(-1.0167936) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8903538) q[0];
sx q[0];
rz(-1.6523424) q[0];
sx q[0];
rz(1.5348877) q[0];
rz(-pi) q[1];
rz(-2.617308) q[2];
sx q[2];
rz(-1.1384083) q[2];
sx q[2];
rz(2.1081032) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6779855) q[1];
sx q[1];
rz(-1.4336168) q[1];
sx q[1];
rz(-0.82171085) q[1];
rz(-1.2891644) q[3];
sx q[3];
rz(-1.6179807) q[3];
sx q[3];
rz(0.55164528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3399405) q[2];
sx q[2];
rz(-1.6093971) q[2];
sx q[2];
rz(2.0920848) q[2];
rz(2.5028051) q[3];
sx q[3];
rz(-2.5103266) q[3];
sx q[3];
rz(1.9558186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8124354) q[0];
sx q[0];
rz(-1.2673459) q[0];
sx q[0];
rz(1.4720434) q[0];
rz(-0.73515785) q[1];
sx q[1];
rz(-0.77886326) q[1];
sx q[1];
rz(0.24681117) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4287764) q[0];
sx q[0];
rz(-1.4220337) q[0];
sx q[0];
rz(-2.2674198) q[0];
rz(-0.47841448) q[2];
sx q[2];
rz(-1.7583499) q[2];
sx q[2];
rz(-2.0699376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.81891638) q[1];
sx q[1];
rz(-2.7671742) q[1];
sx q[1];
rz(0.14426343) q[1];
rz(1.1770583) q[3];
sx q[3];
rz(-0.65440946) q[3];
sx q[3];
rz(-2.8670093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6639158) q[2];
sx q[2];
rz(-1.1185948) q[2];
sx q[2];
rz(-1.6003312) q[2];
rz(2.4345496) q[3];
sx q[3];
rz(-2.0440846) q[3];
sx q[3];
rz(-0.88821205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33070579) q[0];
sx q[0];
rz(-0.72074497) q[0];
sx q[0];
rz(1.8141618) q[0];
rz(1.56303) q[1];
sx q[1];
rz(-0.47416082) q[1];
sx q[1];
rz(0.24838233) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7029593) q[0];
sx q[0];
rz(-2.6002433) q[0];
sx q[0];
rz(-3.0754509) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18647285) q[2];
sx q[2];
rz(-1.4154134) q[2];
sx q[2];
rz(-1.965167) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0548965) q[1];
sx q[1];
rz(-2.5270562) q[1];
sx q[1];
rz(-2.8413248) q[1];
rz(1.3316657) q[3];
sx q[3];
rz(-0.93591792) q[3];
sx q[3];
rz(2.2374416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43626943) q[2];
sx q[2];
rz(-0.98781172) q[2];
sx q[2];
rz(-2.3441337) q[2];
rz(2.752839) q[3];
sx q[3];
rz(-2.5377486) q[3];
sx q[3];
rz(2.6388772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0080863) q[0];
sx q[0];
rz(-0.07645034) q[0];
sx q[0];
rz(-1.3457993) q[0];
rz(-1.0812409) q[1];
sx q[1];
rz(-1.2370279) q[1];
sx q[1];
rz(3.016901) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4437618) q[0];
sx q[0];
rz(-0.19653453) q[0];
sx q[0];
rz(0.4561119) q[0];
rz(-2.5325534) q[2];
sx q[2];
rz(-1.8134724) q[2];
sx q[2];
rz(0.86953029) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2248968) q[1];
sx q[1];
rz(-1.5142913) q[1];
sx q[1];
rz(0.7373666) q[1];
x q[2];
rz(1.7894621) q[3];
sx q[3];
rz(-1.339774) q[3];
sx q[3];
rz(2.4276589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6283915) q[2];
sx q[2];
rz(-1.9548364) q[2];
sx q[2];
rz(0.61895269) q[2];
rz(-1.0533054) q[3];
sx q[3];
rz(-2.9635933) q[3];
sx q[3];
rz(2.2119904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58182794) q[0];
sx q[0];
rz(-1.3904089) q[0];
sx q[0];
rz(-2.0986309) q[0];
rz(0.46328059) q[1];
sx q[1];
rz(-2.0279341) q[1];
sx q[1];
rz(1.0707062) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8585513) q[0];
sx q[0];
rz(-1.0867449) q[0];
sx q[0];
rz(1.5714684) q[0];
x q[1];
rz(2.7878739) q[2];
sx q[2];
rz(-0.36834799) q[2];
sx q[2];
rz(-0.081920953) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5412022) q[1];
sx q[1];
rz(-1.4711079) q[1];
sx q[1];
rz(-2.3669835) q[1];
x q[2];
rz(-1.5911063) q[3];
sx q[3];
rz(-1.2386525) q[3];
sx q[3];
rz(-2.9796245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7730007) q[2];
sx q[2];
rz(-0.7395145) q[2];
sx q[2];
rz(-0.32361844) q[2];
rz(2.1598024) q[3];
sx q[3];
rz(-2.2798645) q[3];
sx q[3];
rz(-1.0872844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48802808) q[0];
sx q[0];
rz(-1.7249148) q[0];
sx q[0];
rz(-1.9352242) q[0];
rz(-1.2127097) q[1];
sx q[1];
rz(-2.2884463) q[1];
sx q[1];
rz(0.94747296) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8038717) q[0];
sx q[0];
rz(-2.5944355) q[0];
sx q[0];
rz(1.7659811) q[0];
rz(-1.7296373) q[2];
sx q[2];
rz(-1.7454595) q[2];
sx q[2];
rz(2.5607002) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.84711134) q[1];
sx q[1];
rz(-1.0080907) q[1];
sx q[1];
rz(0.28114762) q[1];
rz(-pi) q[2];
rz(-0.34668563) q[3];
sx q[3];
rz(-0.83804916) q[3];
sx q[3];
rz(0.29688641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.56132135) q[2];
sx q[2];
rz(-1.3803955) q[2];
sx q[2];
rz(0.46978152) q[2];
rz(1.3011159) q[3];
sx q[3];
rz(-1.7062635) q[3];
sx q[3];
rz(2.8619213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25205055) q[0];
sx q[0];
rz(-0.38953504) q[0];
sx q[0];
rz(-1.8126194) q[0];
rz(2.3503616) q[1];
sx q[1];
rz(-2.8104517) q[1];
sx q[1];
rz(2.9387617) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7751574) q[0];
sx q[0];
rz(-3.0946819) q[0];
sx q[0];
rz(-0.58880083) q[0];
x q[1];
rz(-2.9183396) q[2];
sx q[2];
rz(-1.7749783) q[2];
sx q[2];
rz(0.60165652) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5962332) q[1];
sx q[1];
rz(-1.6714449) q[1];
sx q[1];
rz(-3.0920995) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6796954) q[3];
sx q[3];
rz(-1.7089205) q[3];
sx q[3];
rz(1.4702597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.41708502) q[2];
sx q[2];
rz(-0.26792002) q[2];
sx q[2];
rz(-2.1450796) q[2];
rz(-2.7881682) q[3];
sx q[3];
rz(-2.3963908) q[3];
sx q[3];
rz(-2.3341808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080169454) q[0];
sx q[0];
rz(-0.81289476) q[0];
sx q[0];
rz(2.9598575) q[0];
rz(-0.043047992) q[1];
sx q[1];
rz(-0.64518607) q[1];
sx q[1];
rz(0.28082401) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40598665) q[0];
sx q[0];
rz(-1.1445023) q[0];
sx q[0];
rz(0.50174016) q[0];
rz(-pi) q[1];
rz(0.21364613) q[2];
sx q[2];
rz(-1.3823969) q[2];
sx q[2];
rz(-1.633916) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7933958) q[1];
sx q[1];
rz(-1.5131725) q[1];
sx q[1];
rz(-0.036199526) q[1];
x q[2];
rz(-0.19631581) q[3];
sx q[3];
rz(-1.8972978) q[3];
sx q[3];
rz(1.452009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8250371) q[2];
sx q[2];
rz(-1.2544158) q[2];
sx q[2];
rz(2.5184856) q[2];
rz(2.1394219) q[3];
sx q[3];
rz(-1.3391756) q[3];
sx q[3];
rz(0.56308693) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4939209) q[0];
sx q[0];
rz(-1.5681842) q[0];
sx q[0];
rz(1.6012123) q[0];
rz(-2.2676246) q[1];
sx q[1];
rz(-1.0653492) q[1];
sx q[1];
rz(-3.031562) q[1];
rz(-0.80550823) q[2];
sx q[2];
rz(-1.1136354) q[2];
sx q[2];
rz(-1.8352933) q[2];
rz(-0.86250967) q[3];
sx q[3];
rz(-0.52048341) q[3];
sx q[3];
rz(-2.4181548) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
