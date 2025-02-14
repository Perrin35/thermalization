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
rz(-1.018723) q[0];
sx q[0];
rz(-0.85917226) q[0];
sx q[0];
rz(2.3265042) q[0];
rz(1.9563142) q[1];
sx q[1];
rz(-1.7307245) q[1];
sx q[1];
rz(-1.0676395) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51781228) q[0];
sx q[0];
rz(-1.5945487) q[0];
sx q[0];
rz(-2.0022814) q[0];
rz(2.6175628) q[2];
sx q[2];
rz(-1.9240148) q[2];
sx q[2];
rz(-2.2729682) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8274535) q[1];
sx q[1];
rz(-1.2541391) q[1];
sx q[1];
rz(-0.10152557) q[1];
x q[2];
rz(-0.88495636) q[3];
sx q[3];
rz(-1.7772632) q[3];
sx q[3];
rz(-1.5625169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.69433576) q[2];
sx q[2];
rz(-0.7434291) q[2];
sx q[2];
rz(1.2747964) q[2];
rz(-2.6796807) q[3];
sx q[3];
rz(-2.4670944) q[3];
sx q[3];
rz(-1.1326724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3502515) q[0];
sx q[0];
rz(-0.30838648) q[0];
sx q[0];
rz(1.8885008) q[0];
rz(-2.99627) q[1];
sx q[1];
rz(-1.3959613) q[1];
sx q[1];
rz(-2.0504418) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50877198) q[0];
sx q[0];
rz(-2.6592386) q[0];
sx q[0];
rz(-1.095196) q[0];
rz(-pi) q[1];
rz(0.45812313) q[2];
sx q[2];
rz(-1.8964185) q[2];
sx q[2];
rz(1.9389012) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.67790612) q[1];
sx q[1];
rz(-0.27707252) q[1];
sx q[1];
rz(0.67986791) q[1];
rz(-pi) q[2];
rz(1.1717779) q[3];
sx q[3];
rz(-1.1091091) q[3];
sx q[3];
rz(0.53129133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.48065177) q[2];
sx q[2];
rz(-1.7512243) q[2];
sx q[2];
rz(-2.6118028) q[2];
rz(-2.3482813) q[3];
sx q[3];
rz(-1.6042234) q[3];
sx q[3];
rz(-2.5180499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48524258) q[0];
sx q[0];
rz(-0.95004496) q[0];
sx q[0];
rz(-0.52870885) q[0];
rz(-0.54620019) q[1];
sx q[1];
rz(-2.1827953) q[1];
sx q[1];
rz(-0.34034696) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2961194) q[0];
sx q[0];
rz(-1.8438135) q[0];
sx q[0];
rz(2.8014328) q[0];
rz(-pi) q[1];
rz(-1.7163926) q[2];
sx q[2];
rz(-1.8032089) q[2];
sx q[2];
rz(0.53047859) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3888549) q[1];
sx q[1];
rz(-1.5675263) q[1];
sx q[1];
rz(-3.0470303) q[1];
rz(-pi) q[2];
rz(-2.3119218) q[3];
sx q[3];
rz(-1.2874914) q[3];
sx q[3];
rz(1.0782575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9693552) q[2];
sx q[2];
rz(-0.41219553) q[2];
sx q[2];
rz(2.7117512) q[2];
rz(-2.3853081) q[3];
sx q[3];
rz(-0.1736621) q[3];
sx q[3];
rz(-2.3156796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7569358) q[0];
sx q[0];
rz(-1.5058368) q[0];
sx q[0];
rz(-1.4935619) q[0];
rz(0.76796302) q[1];
sx q[1];
rz(-2.6710644) q[1];
sx q[1];
rz(0.54642645) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38619216) q[0];
sx q[0];
rz(-1.6419171) q[0];
sx q[0];
rz(3.0767308) q[0];
rz(-pi) q[1];
rz(-2.3713263) q[2];
sx q[2];
rz(-0.67521836) q[2];
sx q[2];
rz(2.9807621) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.029303251) q[1];
sx q[1];
rz(-1.7186856) q[1];
sx q[1];
rz(-1.4724971) q[1];
rz(-pi) q[2];
rz(0.68164556) q[3];
sx q[3];
rz(-2.1053616) q[3];
sx q[3];
rz(-0.93972423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4297318) q[2];
sx q[2];
rz(-1.4873361) q[2];
sx q[2];
rz(-2.8374953) q[2];
rz(1.3652623) q[3];
sx q[3];
rz(-1.920776) q[3];
sx q[3];
rz(2.064866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51233184) q[0];
sx q[0];
rz(-2.6973695) q[0];
sx q[0];
rz(-0.00057922676) q[0];
rz(-2.0594788) q[1];
sx q[1];
rz(-2.6172456) q[1];
sx q[1];
rz(2.2023315) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9094641) q[0];
sx q[0];
rz(-0.8117903) q[0];
sx q[0];
rz(-2.5215197) q[0];
rz(-pi) q[1];
rz(-2.0505191) q[2];
sx q[2];
rz(-0.66893286) q[2];
sx q[2];
rz(0.26593966) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.68785948) q[1];
sx q[1];
rz(-2.5304768) q[1];
sx q[1];
rz(-1.1920209) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95114189) q[3];
sx q[3];
rz(-0.37887805) q[3];
sx q[3];
rz(-0.71228107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1028563) q[2];
sx q[2];
rz(-1.880371) q[2];
sx q[2];
rz(-0.2612513) q[2];
rz(0.86483613) q[3];
sx q[3];
rz(-1.7295001) q[3];
sx q[3];
rz(-2.849546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.29702) q[0];
sx q[0];
rz(-1.7140056) q[0];
sx q[0];
rz(-2.936506) q[0];
rz(-2.0948441) q[1];
sx q[1];
rz(-1.3958684) q[1];
sx q[1];
rz(-2.7632025) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16941026) q[0];
sx q[0];
rz(-0.43146389) q[0];
sx q[0];
rz(-2.0305102) q[0];
x q[1];
rz(-0.47176265) q[2];
sx q[2];
rz(-1.7182351) q[2];
sx q[2];
rz(2.9094537) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7582963) q[1];
sx q[1];
rz(-1.2219011) q[1];
sx q[1];
rz(2.1339244) q[1];
rz(-pi) q[2];
rz(-2.4491485) q[3];
sx q[3];
rz(-1.2298349) q[3];
sx q[3];
rz(2.8657262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.64688524) q[2];
sx q[2];
rz(-1.9834221) q[2];
sx q[2];
rz(-1.0542487) q[2];
rz(0.65822893) q[3];
sx q[3];
rz(-0.64526486) q[3];
sx q[3];
rz(0.48404199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.9996416) q[0];
sx q[0];
rz(-1.9331837) q[0];
sx q[0];
rz(2.5010338) q[0];
rz(-0.41796747) q[1];
sx q[1];
rz(-1.2203981) q[1];
sx q[1];
rz(2.3366065) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98950878) q[0];
sx q[0];
rz(-1.594127) q[0];
sx q[0];
rz(-2.4341466) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80193582) q[2];
sx q[2];
rz(-1.2705497) q[2];
sx q[2];
rz(1.9682457) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5857081) q[1];
sx q[1];
rz(-1.5991365) q[1];
sx q[1];
rz(3.1217087) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42607362) q[3];
sx q[3];
rz(-2.5059627) q[3];
sx q[3];
rz(2.9866708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5927222) q[2];
sx q[2];
rz(-0.9950811) q[2];
sx q[2];
rz(0.76118809) q[2];
rz(-0.74782863) q[3];
sx q[3];
rz(-1.3176094) q[3];
sx q[3];
rz(0.67659155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6250896) q[0];
sx q[0];
rz(-0.28446063) q[0];
sx q[0];
rz(1.6647343) q[0];
rz(-2.6853216) q[1];
sx q[1];
rz(-1.7218593) q[1];
sx q[1];
rz(-0.87108535) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7719771) q[0];
sx q[0];
rz(-0.61423683) q[0];
sx q[0];
rz(-0.02158879) q[0];
x q[1];
rz(-0.41042491) q[2];
sx q[2];
rz(-0.3730118) q[2];
sx q[2];
rz(-2.5665064) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2312517) q[1];
sx q[1];
rz(-1.3518855) q[1];
sx q[1];
rz(-0.33556767) q[1];
rz(0.060163946) q[3];
sx q[3];
rz(-1.9964661) q[3];
sx q[3];
rz(-2.343246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2386834) q[2];
sx q[2];
rz(-2.6148655) q[2];
sx q[2];
rz(-0.61863679) q[2];
rz(-3.1213308) q[3];
sx q[3];
rz(-2.198115) q[3];
sx q[3];
rz(-2.3616135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0778377) q[0];
sx q[0];
rz(-1.5851333) q[0];
sx q[0];
rz(0.42386398) q[0];
rz(-2.9546812) q[1];
sx q[1];
rz(-2.321545) q[1];
sx q[1];
rz(-1.4580457) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2423201) q[0];
sx q[0];
rz(-1.5566711) q[0];
sx q[0];
rz(-0.11798162) q[0];
x q[1];
rz(0.27159043) q[2];
sx q[2];
rz(-1.0729861) q[2];
sx q[2];
rz(-1.9099727) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8193389) q[1];
sx q[1];
rz(-1.5194494) q[1];
sx q[1];
rz(3.0027185) q[1];
rz(-1.4814754) q[3];
sx q[3];
rz(-1.2161939) q[3];
sx q[3];
rz(1.3077298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.75108782) q[2];
sx q[2];
rz(-2.0743399) q[2];
sx q[2];
rz(1.0941774) q[2];
rz(1.68082) q[3];
sx q[3];
rz(-1.7655617) q[3];
sx q[3];
rz(1.8028397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.3222892) q[0];
sx q[0];
rz(-1.3811454) q[0];
sx q[0];
rz(-1.639701) q[0];
rz(-2.8627401) q[1];
sx q[1];
rz(-0.96062213) q[1];
sx q[1];
rz(2.1174812) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49660027) q[0];
sx q[0];
rz(-1.6628213) q[0];
sx q[0];
rz(2.651398) q[0];
x q[1];
rz(1.7488519) q[2];
sx q[2];
rz(-2.8626056) q[2];
sx q[2];
rz(1.6131608) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1361724) q[1];
sx q[1];
rz(-1.7917624) q[1];
sx q[1];
rz(-1.3653838) q[1];
x q[2];
rz(1.2028221) q[3];
sx q[3];
rz(-2.362613) q[3];
sx q[3];
rz(-0.72768962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2037105) q[2];
sx q[2];
rz(-1.8546162) q[2];
sx q[2];
rz(2.6046275) q[2];
rz(-1.9133866) q[3];
sx q[3];
rz(-0.95913404) q[3];
sx q[3];
rz(-2.8502407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.1163597) q[0];
sx q[0];
rz(-1.7740842) q[0];
sx q[0];
rz(-2.0728003) q[0];
rz(2.4346726) q[1];
sx q[1];
rz(-1.128935) q[1];
sx q[1];
rz(3.1405906) q[1];
rz(0.40590103) q[2];
sx q[2];
rz(-1.7465456) q[2];
sx q[2];
rz(-3.0301574) q[2];
rz(1.3523921) q[3];
sx q[3];
rz(-1.7159749) q[3];
sx q[3];
rz(0.80930474) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
