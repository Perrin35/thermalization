OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5443213) q[0];
sx q[0];
rz(-2.7379524) q[0];
sx q[0];
rz(0.37024745) q[0];
rz(0.20180841) q[1];
sx q[1];
rz(-1.2528074) q[1];
sx q[1];
rz(1.7226146) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0231853) q[0];
sx q[0];
rz(-1.4010795) q[0];
sx q[0];
rz(0.24766185) q[0];
rz(-0.51214829) q[2];
sx q[2];
rz(-1.3379659) q[2];
sx q[2];
rz(-1.6793959) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.943489) q[1];
sx q[1];
rz(-1.5183105) q[1];
sx q[1];
rz(-2.1211038) q[1];
rz(-pi) q[2];
rz(1.8889514) q[3];
sx q[3];
rz(-2.3080024) q[3];
sx q[3];
rz(2.6300501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7449164) q[2];
sx q[2];
rz(-1.4531206) q[2];
sx q[2];
rz(0.35787004) q[2];
rz(0.19168028) q[3];
sx q[3];
rz(-2.7087757) q[3];
sx q[3];
rz(-0.55309692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-3.1382004) q[0];
sx q[0];
rz(-1.9906185) q[0];
sx q[0];
rz(-0.068280846) q[0];
rz(-1.0066907) q[1];
sx q[1];
rz(-0.12193646) q[1];
sx q[1];
rz(-2.4904747) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7333974) q[0];
sx q[0];
rz(-1.7719643) q[0];
sx q[0];
rz(2.4961619) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5947371) q[2];
sx q[2];
rz(-2.0964453) q[2];
sx q[2];
rz(2.0999694) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.56602851) q[1];
sx q[1];
rz(-2.063942) q[1];
sx q[1];
rz(0.50599392) q[1];
rz(2.7495456) q[3];
sx q[3];
rz(-0.52467504) q[3];
sx q[3];
rz(2.3384561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6887001) q[2];
sx q[2];
rz(-0.15658997) q[2];
sx q[2];
rz(2.2382656) q[2];
rz(2.3870758) q[3];
sx q[3];
rz(-1.6468331) q[3];
sx q[3];
rz(-2.8675458) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84905255) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(-2.6480411) q[0];
rz(-1.4986528) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(-1.0292056) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2632336) q[0];
sx q[0];
rz(-1.2326816) q[0];
sx q[0];
rz(-2.6792206) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35691397) q[2];
sx q[2];
rz(-1.9682353) q[2];
sx q[2];
rz(2.7473161) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5310276) q[1];
sx q[1];
rz(-2.768801) q[1];
sx q[1];
rz(2.5337266) q[1];
x q[2];
rz(2.1152705) q[3];
sx q[3];
rz(-1.4635651) q[3];
sx q[3];
rz(2.2162007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9323953) q[2];
sx q[2];
rz(-2.5019427) q[2];
sx q[2];
rz(-1.2970682) q[2];
rz(2.5629937) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(2.4659618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9349174) q[0];
sx q[0];
rz(-0.6466372) q[0];
sx q[0];
rz(2.7329965) q[0];
rz(-1.714255) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(2.2600007) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5170167) q[0];
sx q[0];
rz(-1.3991742) q[0];
sx q[0];
rz(-0.73715985) q[0];
rz(-pi) q[1];
rz(0.84882952) q[2];
sx q[2];
rz(-1.0366882) q[2];
sx q[2];
rz(1.3918849) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3410014) q[1];
sx q[1];
rz(-2.2644271) q[1];
sx q[1];
rz(-2.4800406) q[1];
rz(-0.41739695) q[3];
sx q[3];
rz(-1.5098803) q[3];
sx q[3];
rz(-3.0792189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0980229) q[2];
sx q[2];
rz(-1.4901525) q[2];
sx q[2];
rz(2.0289452) q[2];
rz(-2.5860795) q[3];
sx q[3];
rz(-1.3534618) q[3];
sx q[3];
rz(-1.6120733) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1495789) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(1.7768815) q[0];
rz(0.31750202) q[1];
sx q[1];
rz(-0.96173871) q[1];
sx q[1];
rz(0.11725765) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1299767) q[0];
sx q[0];
rz(-2.5419309) q[0];
sx q[0];
rz(0.2325124) q[0];
rz(-pi) q[1];
rz(-2.8357382) q[2];
sx q[2];
rz(-0.62539414) q[2];
sx q[2];
rz(1.9276227) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6834516) q[1];
sx q[1];
rz(-1.3740173) q[1];
sx q[1];
rz(1.7276006) q[1];
x q[2];
rz(-1.6860784) q[3];
sx q[3];
rz(-0.59024631) q[3];
sx q[3];
rz(0.78199996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.34510288) q[2];
sx q[2];
rz(-1.3619276) q[2];
sx q[2];
rz(-1.3797181) q[2];
rz(-1.1896677) q[3];
sx q[3];
rz(-0.15933557) q[3];
sx q[3];
rz(3.0670847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9538486) q[0];
sx q[0];
rz(-1.501361) q[0];
sx q[0];
rz(-0.70621079) q[0];
rz(-1.1114936) q[1];
sx q[1];
rz(-0.62875426) q[1];
sx q[1];
rz(-3.0336753) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6861434) q[0];
sx q[0];
rz(-2.3685072) q[0];
sx q[0];
rz(-0.42868741) q[0];
rz(-1.582167) q[2];
sx q[2];
rz(-1.5307384) q[2];
sx q[2];
rz(-0.96206059) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2006827) q[1];
sx q[1];
rz(-2.4001277) q[1];
sx q[1];
rz(0.42592589) q[1];
rz(-0.48072731) q[3];
sx q[3];
rz(-0.93989621) q[3];
sx q[3];
rz(2.7709099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8191021) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(2.0325913) q[2];
rz(1.8479944) q[3];
sx q[3];
rz(-1.785991) q[3];
sx q[3];
rz(-3.0343645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.368822) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(3.1266881) q[0];
rz(2.7203454) q[1];
sx q[1];
rz(-1.0528456) q[1];
sx q[1];
rz(-2.3419535) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5937724) q[0];
sx q[0];
rz(-0.9043588) q[0];
sx q[0];
rz(-1.244546) q[0];
rz(2.4104426) q[2];
sx q[2];
rz(-0.41229782) q[2];
sx q[2];
rz(-1.2286124) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0529424) q[1];
sx q[1];
rz(-0.155641) q[1];
sx q[1];
rz(2.774653) q[1];
x q[2];
rz(-0.13749595) q[3];
sx q[3];
rz(-1.0097479) q[3];
sx q[3];
rz(-1.6784799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.556276) q[2];
sx q[2];
rz(-0.63627807) q[2];
sx q[2];
rz(0.9220534) q[2];
rz(-1.3098035) q[3];
sx q[3];
rz(-1.9366879) q[3];
sx q[3];
rz(-2.183765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9668982) q[0];
sx q[0];
rz(-2.0860724) q[0];
sx q[0];
rz(0.2640557) q[0];
rz(-1.7407725) q[1];
sx q[1];
rz(-1.4512647) q[1];
sx q[1];
rz(-0.31731269) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2469651) q[0];
sx q[0];
rz(-1.6435701) q[0];
sx q[0];
rz(-1.1829681) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9810901) q[2];
sx q[2];
rz(-1.9151033) q[2];
sx q[2];
rz(2.3366994) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.322284) q[1];
sx q[1];
rz(-1.6849736) q[1];
sx q[1];
rz(-1.6627922) q[1];
x q[2];
rz(-1.7371014) q[3];
sx q[3];
rz(-2.2095223) q[3];
sx q[3];
rz(-2.2685662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.43549609) q[2];
sx q[2];
rz(-2.2075682) q[2];
sx q[2];
rz(1.6395817) q[2];
rz(-0.25029415) q[3];
sx q[3];
rz(-1.7264265) q[3];
sx q[3];
rz(-1.9992453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2450927) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(-1.4244351) q[0];
rz(-0.4793438) q[1];
sx q[1];
rz(-1.665325) q[1];
sx q[1];
rz(0.11553484) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.696366) q[0];
sx q[0];
rz(-2.1292902) q[0];
sx q[0];
rz(-2.2541056) q[0];
x q[1];
rz(1.8477511) q[2];
sx q[2];
rz(-0.87522725) q[2];
sx q[2];
rz(0.19687523) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86347307) q[1];
sx q[1];
rz(-1.9895456) q[1];
sx q[1];
rz(-0.31660415) q[1];
rz(-pi) q[2];
rz(2.9980738) q[3];
sx q[3];
rz(-0.70129881) q[3];
sx q[3];
rz(1.6539751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4775548) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(1.4257365) q[2];
rz(1.4005631) q[3];
sx q[3];
rz(-1.0481342) q[3];
sx q[3];
rz(2.8588296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7857159) q[0];
sx q[0];
rz(-2.4079005) q[0];
sx q[0];
rz(2.8724331) q[0];
rz(2.0458938) q[1];
sx q[1];
rz(-0.91047374) q[1];
sx q[1];
rz(-1.7620618) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9513959) q[0];
sx q[0];
rz(-1.405389) q[0];
sx q[0];
rz(-0.10683807) q[0];
rz(-pi) q[1];
rz(-1.3921521) q[2];
sx q[2];
rz(-2.7182455) q[2];
sx q[2];
rz(-1.64738) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.79593149) q[1];
sx q[1];
rz(-0.17594166) q[1];
sx q[1];
rz(-1.7453341) q[1];
rz(-0.33196253) q[3];
sx q[3];
rz(-1.8649857) q[3];
sx q[3];
rz(0.709579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0925838) q[2];
sx q[2];
rz(-2.9276586) q[2];
sx q[2];
rz(1.6258378) q[2];
rz(-1.1670636) q[3];
sx q[3];
rz(-1.556282) q[3];
sx q[3];
rz(1.0313755) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9243069) q[0];
sx q[0];
rz(-1.3047682) q[0];
sx q[0];
rz(-2.7346942) q[0];
rz(-0.45463195) q[1];
sx q[1];
rz(-1.1063207) q[1];
sx q[1];
rz(2.8938821) q[1];
rz(1.0977279) q[2];
sx q[2];
rz(-1.4164783) q[2];
sx q[2];
rz(2.867792) q[2];
rz(3.0804844) q[3];
sx q[3];
rz(-1.8295049) q[3];
sx q[3];
rz(-2.6873333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
