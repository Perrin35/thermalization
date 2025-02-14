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
rz(1.4003657) q[0];
sx q[0];
rz(-0.23973149) q[0];
sx q[0];
rz(-2.8170407) q[0];
rz(5.3311081) q[1];
sx q[1];
rz(6.0573112) q[1];
sx q[1];
rz(4.4499302) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1150779) q[0];
sx q[0];
rz(-0.9008207) q[0];
sx q[0];
rz(1.7448533) q[0];
rz(-1.7003072) q[2];
sx q[2];
rz(-1.0900094) q[2];
sx q[2];
rz(-0.12809424) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1349012) q[1];
sx q[1];
rz(-2.6009702) q[1];
sx q[1];
rz(2.0903793) q[1];
x q[2];
rz(1.4300554) q[3];
sx q[3];
rz(-2.2651197) q[3];
sx q[3];
rz(1.3545881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41210458) q[2];
sx q[2];
rz(-1.8821913) q[2];
sx q[2];
rz(0.35120249) q[2];
rz(1.2900194) q[3];
sx q[3];
rz(-1.131564) q[3];
sx q[3];
rz(-1.2235519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7555162) q[0];
sx q[0];
rz(-1.9367243) q[0];
sx q[0];
rz(-2.2112041) q[0];
rz(1.3049841) q[1];
sx q[1];
rz(-0.84140673) q[1];
sx q[1];
rz(-2.5321541) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6822426) q[0];
sx q[0];
rz(-1.9960023) q[0];
sx q[0];
rz(0.11328477) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8571505) q[2];
sx q[2];
rz(-2.831818) q[2];
sx q[2];
rz(-1.3510493) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7623065) q[1];
sx q[1];
rz(-1.9348782) q[1];
sx q[1];
rz(1.9696139) q[1];
rz(0.16709631) q[3];
sx q[3];
rz(-2.1347988) q[3];
sx q[3];
rz(0.43844863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.096752) q[2];
sx q[2];
rz(-2.1614306) q[2];
sx q[2];
rz(-0.074782221) q[2];
rz(-0.52037248) q[3];
sx q[3];
rz(-2.4986391) q[3];
sx q[3];
rz(-0.61409942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5388913) q[0];
sx q[0];
rz(-0.7170054) q[0];
sx q[0];
rz(0.30211788) q[0];
rz(-0.27613861) q[1];
sx q[1];
rz(-1.8243676) q[1];
sx q[1];
rz(1.709323) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0841504) q[0];
sx q[0];
rz(-1.0330811) q[0];
sx q[0];
rz(2.6938963) q[0];
rz(-pi) q[1];
rz(-0.59467043) q[2];
sx q[2];
rz(-1.7886908) q[2];
sx q[2];
rz(-1.7105107) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.36784962) q[1];
sx q[1];
rz(-2.5605695) q[1];
sx q[1];
rz(1.2817205) q[1];
rz(1.7114221) q[3];
sx q[3];
rz(-0.46984497) q[3];
sx q[3];
rz(1.2805366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8351195) q[2];
sx q[2];
rz(-0.46572954) q[2];
sx q[2];
rz(-1.2960557) q[2];
rz(0.45201388) q[3];
sx q[3];
rz(-2.0343503) q[3];
sx q[3];
rz(-0.38255102) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52962676) q[0];
sx q[0];
rz(-1.1671678) q[0];
sx q[0];
rz(1.3324598) q[0];
rz(2.0924856) q[1];
sx q[1];
rz(-1.0915979) q[1];
sx q[1];
rz(1.5886935) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5063254) q[0];
sx q[0];
rz(-2.679279) q[0];
sx q[0];
rz(2.4145831) q[0];
rz(-2.0715782) q[2];
sx q[2];
rz(-1.5588648) q[2];
sx q[2];
rz(-2.3119761) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5415948) q[1];
sx q[1];
rz(-1.7686971) q[1];
sx q[1];
rz(2.2210414) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53557204) q[3];
sx q[3];
rz(-0.93129292) q[3];
sx q[3];
rz(-0.65062338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4730452) q[2];
sx q[2];
rz(-2.1210402) q[2];
sx q[2];
rz(-0.2505396) q[2];
rz(0.9969095) q[3];
sx q[3];
rz(-2.0018115) q[3];
sx q[3];
rz(-0.38058773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32960358) q[0];
sx q[0];
rz(-2.6530837) q[0];
sx q[0];
rz(1.3478152) q[0];
rz(2.7373121) q[1];
sx q[1];
rz(-0.75899044) q[1];
sx q[1];
rz(1.1291198) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30325748) q[0];
sx q[0];
rz(-1.9651439) q[0];
sx q[0];
rz(-2.3189937) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2363844) q[2];
sx q[2];
rz(-1.7844177) q[2];
sx q[2];
rz(2.7366432) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3394062) q[1];
sx q[1];
rz(-1.4660343) q[1];
sx q[1];
rz(0.89749194) q[1];
rz(-pi) q[2];
rz(-1.681116) q[3];
sx q[3];
rz(-1.8838143) q[3];
sx q[3];
rz(2.4576994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79444844) q[2];
sx q[2];
rz(-0.930154) q[2];
sx q[2];
rz(1.244119) q[2];
rz(2.0813023) q[3];
sx q[3];
rz(-2.5131707) q[3];
sx q[3];
rz(2.8384143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4853972) q[0];
sx q[0];
rz(-3.0363016) q[0];
sx q[0];
rz(-0.42718497) q[0];
rz(-3.1071013) q[1];
sx q[1];
rz(-1.5723615) q[1];
sx q[1];
rz(0.010206612) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9968352) q[0];
sx q[0];
rz(-0.0029276927) q[0];
sx q[0];
rz(-2.5064431) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9811007) q[2];
sx q[2];
rz(-0.94466034) q[2];
sx q[2];
rz(0.4753091) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.806655) q[1];
sx q[1];
rz(-0.68938556) q[1];
sx q[1];
rz(-1.6720119) q[1];
x q[2];
rz(-2.5764546) q[3];
sx q[3];
rz(-2.0083154) q[3];
sx q[3];
rz(1.138162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5222142) q[2];
sx q[2];
rz(-0.41316119) q[2];
sx q[2];
rz(-0.79356066) q[2];
rz(-2.2788952) q[3];
sx q[3];
rz(-1.5746652) q[3];
sx q[3];
rz(0.70393744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4392387) q[0];
sx q[0];
rz(-0.84504253) q[0];
sx q[0];
rz(-2.0080361) q[0];
rz(2.0639065) q[1];
sx q[1];
rz(-2.6262941) q[1];
sx q[1];
rz(-2.8394707) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7493499) q[0];
sx q[0];
rz(-1.8729405) q[0];
sx q[0];
rz(-2.7878615) q[0];
rz(2.0321192) q[2];
sx q[2];
rz(-1.0240842) q[2];
sx q[2];
rz(2.6448768) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6670468) q[1];
sx q[1];
rz(-1.3549558) q[1];
sx q[1];
rz(0.53498909) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68639836) q[3];
sx q[3];
rz(-2.0455615) q[3];
sx q[3];
rz(1.7976185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9083531) q[2];
sx q[2];
rz(-2.2519604) q[2];
sx q[2];
rz(1.3471777) q[2];
rz(1.8917278) q[3];
sx q[3];
rz(-1.9439387) q[3];
sx q[3];
rz(3.0344322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5938479) q[0];
sx q[0];
rz(-2.7084454) q[0];
sx q[0];
rz(-3.0331842) q[0];
rz(-1.0477061) q[1];
sx q[1];
rz(-0.81755081) q[1];
sx q[1];
rz(-1.0188867) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86025809) q[0];
sx q[0];
rz(-1.0509914) q[0];
sx q[0];
rz(-2.5297574) q[0];
rz(-pi) q[1];
rz(-0.0068130612) q[2];
sx q[2];
rz(-2.3802813) q[2];
sx q[2];
rz(0.50466621) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4561676) q[1];
sx q[1];
rz(-1.9677094) q[1];
sx q[1];
rz(0.49938582) q[1];
rz(-pi) q[2];
rz(-1.4298444) q[3];
sx q[3];
rz(-2.3059855) q[3];
sx q[3];
rz(2.6648389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.63742796) q[2];
sx q[2];
rz(-2.1647858) q[2];
sx q[2];
rz(0.43761474) q[2];
rz(2.6484683) q[3];
sx q[3];
rz(-0.92032856) q[3];
sx q[3];
rz(-1.5776618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.86579943) q[0];
sx q[0];
rz(-1.2340622) q[0];
sx q[0];
rz(0.27780521) q[0];
rz(1.5996784) q[1];
sx q[1];
rz(-2.1733687) q[1];
sx q[1];
rz(2.7154198) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0876678) q[0];
sx q[0];
rz(-0.024294596) q[0];
sx q[0];
rz(0.75456516) q[0];
x q[1];
rz(-2.4075899) q[2];
sx q[2];
rz(-0.91911572) q[2];
sx q[2];
rz(0.82254788) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.81282367) q[1];
sx q[1];
rz(-1.4534833) q[1];
sx q[1];
rz(-1.3552865) q[1];
rz(-0.83287333) q[3];
sx q[3];
rz(-0.98372059) q[3];
sx q[3];
rz(0.47490109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6466732) q[2];
sx q[2];
rz(-0.16725954) q[2];
sx q[2];
rz(2.0885928) q[2];
rz(2.7827175) q[3];
sx q[3];
rz(-1.5263298) q[3];
sx q[3];
rz(1.0927965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36295715) q[0];
sx q[0];
rz(-0.21720049) q[0];
sx q[0];
rz(0.92078513) q[0];
rz(0.50865632) q[1];
sx q[1];
rz(-0.76728907) q[1];
sx q[1];
rz(0.42053929) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8201308) q[0];
sx q[0];
rz(-2.252914) q[0];
sx q[0];
rz(2.2338623) q[0];
x q[1];
rz(1.2447717) q[2];
sx q[2];
rz(-2.4851126) q[2];
sx q[2];
rz(-3.1132389) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0691472) q[1];
sx q[1];
rz(-2.0305567) q[1];
sx q[1];
rz(2.2531177) q[1];
rz(-1.4484506) q[3];
sx q[3];
rz(-2.4469355) q[3];
sx q[3];
rz(-2.8784424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4113808) q[2];
sx q[2];
rz(-0.79019848) q[2];
sx q[2];
rz(-0.31663695) q[2];
rz(-0.5591048) q[3];
sx q[3];
rz(-2.5059301) q[3];
sx q[3];
rz(0.81812286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.485514) q[0];
sx q[0];
rz(-2.06388) q[0];
sx q[0];
rz(-2.0583454) q[0];
rz(2.6651233) q[1];
sx q[1];
rz(-2.10119) q[1];
sx q[1];
rz(1.4269921) q[1];
rz(2.3219288) q[2];
sx q[2];
rz(-1.507187) q[2];
sx q[2];
rz(1.0095467) q[2];
rz(1.9540167) q[3];
sx q[3];
rz(-1.6578309) q[3];
sx q[3];
rz(2.1014392) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
