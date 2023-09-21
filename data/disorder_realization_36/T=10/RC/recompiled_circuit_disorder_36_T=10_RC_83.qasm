OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6498123) q[0];
sx q[0];
rz(-0.28591135) q[0];
sx q[0];
rz(-2.6262992) q[0];
rz(4.4858785) q[1];
sx q[1];
rz(2.9872515) q[1];
sx q[1];
rz(6.8607688) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87969765) q[0];
sx q[0];
rz(-1.7356153) q[0];
sx q[0];
rz(2.4796955) q[0];
rz(-pi) q[1];
rz(-0.75391407) q[2];
sx q[2];
rz(-2.4587817) q[2];
sx q[2];
rz(-2.4726601) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2540993) q[1];
sx q[1];
rz(-1.3026397) q[1];
sx q[1];
rz(-1.1521794) q[1];
x q[2];
rz(-0.76831423) q[3];
sx q[3];
rz(-0.81211219) q[3];
sx q[3];
rz(2.1384359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.704533) q[2];
sx q[2];
rz(-1.5455064) q[2];
sx q[2];
rz(-0.68721592) q[2];
rz(-2.1263188) q[3];
sx q[3];
rz(-1.7679368) q[3];
sx q[3];
rz(-0.12250531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9706443) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(-1.8815536) q[0];
rz(-2.1353703) q[1];
sx q[1];
rz(-2.1496014) q[1];
sx q[1];
rz(-0.84567436) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0613522) q[0];
sx q[0];
rz(-0.70525673) q[0];
sx q[0];
rz(-2.4387226) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0078366) q[2];
sx q[2];
rz(-1.4417366) q[2];
sx q[2];
rz(-1.918902) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80026514) q[1];
sx q[1];
rz(-2.6086573) q[1];
sx q[1];
rz(-0.37377263) q[1];
rz(2.2250697) q[3];
sx q[3];
rz(-0.73361165) q[3];
sx q[3];
rz(0.77378002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3530897) q[2];
sx q[2];
rz(-0.22558364) q[2];
sx q[2];
rz(-0.4804002) q[2];
rz(1.7885615) q[3];
sx q[3];
rz(-1.055911) q[3];
sx q[3];
rz(1.9539179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1903494) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(2.3685266) q[0];
rz(-3.0103325) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(-1.0864331) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4341136) q[0];
sx q[0];
rz(-2.0322324) q[0];
sx q[0];
rz(3.1397318) q[0];
x q[1];
rz(1.418872) q[2];
sx q[2];
rz(-0.65290367) q[2];
sx q[2];
rz(-2.8002847) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.72213458) q[1];
sx q[1];
rz(-0.28885435) q[1];
sx q[1];
rz(-0.018317776) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9757189) q[3];
sx q[3];
rz(-1.7305264) q[3];
sx q[3];
rz(-2.3111642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0125668) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(1.3712937) q[2];
rz(0.38315547) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(0.80254054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6894158) q[0];
sx q[0];
rz(-1.2503662) q[0];
sx q[0];
rz(3.0932328) q[0];
rz(-2.9776749) q[1];
sx q[1];
rz(-0.36968958) q[1];
sx q[1];
rz(-1.6960467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7955129) q[0];
sx q[0];
rz(-0.68329408) q[0];
sx q[0];
rz(-0.28738316) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21638685) q[2];
sx q[2];
rz(-1.6587703) q[2];
sx q[2];
rz(2.3073334) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4353664) q[1];
sx q[1];
rz(-1.5063783) q[1];
sx q[1];
rz(2.9219342) q[1];
rz(-pi) q[2];
rz(1.0501782) q[3];
sx q[3];
rz(-1.319066) q[3];
sx q[3];
rz(0.43581918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.066102862) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(2.1172822) q[2];
rz(-1.5284437) q[3];
sx q[3];
rz(-1.5214835) q[3];
sx q[3];
rz(0.23322341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4399453) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(1.8540927) q[0];
rz(-2.8225186) q[1];
sx q[1];
rz(-1.5417475) q[1];
sx q[1];
rz(0.85420001) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4695278) q[0];
sx q[0];
rz(-0.85692353) q[0];
sx q[0];
rz(-0.010849997) q[0];
rz(-pi) q[1];
rz(1.008026) q[2];
sx q[2];
rz(-0.78283435) q[2];
sx q[2];
rz(-2.7340739) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8728719) q[1];
sx q[1];
rz(-0.99209058) q[1];
sx q[1];
rz(-2.242356) q[1];
x q[2];
rz(2.5913127) q[3];
sx q[3];
rz(-2.0645421) q[3];
sx q[3];
rz(-0.20875904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1896818) q[2];
sx q[2];
rz(-2.5482735) q[2];
sx q[2];
rz(-2.5642776) q[2];
rz(-2.632085) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1381056) q[0];
sx q[0];
rz(-2.0697937) q[0];
sx q[0];
rz(-3.0694718) q[0];
rz(-1.1068608) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(-0.12621005) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5835411) q[0];
sx q[0];
rz(-1.7833033) q[0];
sx q[0];
rz(1.538518) q[0];
rz(-pi) q[1];
rz(0.38527617) q[2];
sx q[2];
rz(-0.81309536) q[2];
sx q[2];
rz(0.77020459) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8773552) q[1];
sx q[1];
rz(-0.74062956) q[1];
sx q[1];
rz(-1.2901558) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8993127) q[3];
sx q[3];
rz(-1.7528755) q[3];
sx q[3];
rz(-1.8329221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2386027) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(-2.3664756) q[2];
rz(0.827968) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(-1.1221788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59259748) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(0.098408498) q[0];
rz(1.1920284) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(-2.5820406) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10683051) q[0];
sx q[0];
rz(-0.97495279) q[0];
sx q[0];
rz(-1.4164657) q[0];
rz(-pi) q[1];
rz(-1.5188811) q[2];
sx q[2];
rz(-1.8318818) q[2];
sx q[2];
rz(1.7920997) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.93058649) q[1];
sx q[1];
rz(-1.6582656) q[1];
sx q[1];
rz(1.7006111) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8758994) q[3];
sx q[3];
rz(-0.64722792) q[3];
sx q[3];
rz(-0.94097394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6825535) q[2];
sx q[2];
rz(-1.8457396) q[2];
sx q[2];
rz(-0.34379488) q[2];
rz(-2.5750459) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(-0.47376537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7664117) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(-2.4108316) q[0];
rz(-0.14239755) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(0.87160814) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5985142) q[0];
sx q[0];
rz(-0.075965479) q[0];
sx q[0];
rz(3.024858) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34198728) q[2];
sx q[2];
rz(-0.42765289) q[2];
sx q[2];
rz(-2.4004186) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3155047) q[1];
sx q[1];
rz(-1.9313889) q[1];
sx q[1];
rz(2.4396067) q[1];
rz(-pi) q[2];
rz(1.3685162) q[3];
sx q[3];
rz(-0.61198046) q[3];
sx q[3];
rz(2.4806541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4153851) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(-2.0020206) q[2];
rz(-1.4987882) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(-0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20294872) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(-1.9198445) q[0];
rz(0.16601673) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(-1.6171914) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.411392) q[0];
sx q[0];
rz(-0.087326614) q[0];
sx q[0];
rz(1.2325531) q[0];
rz(-pi) q[1];
rz(-1.379307) q[2];
sx q[2];
rz(-1.101149) q[2];
sx q[2];
rz(-0.42275235) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.4001273) q[1];
sx q[1];
rz(-1.8997846) q[1];
sx q[1];
rz(0.16727438) q[1];
x q[2];
rz(-2.5025326) q[3];
sx q[3];
rz(-1.8250873) q[3];
sx q[3];
rz(1.1921079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1200072) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(-0.35153708) q[2];
rz(-2.0848138) q[3];
sx q[3];
rz(-0.52962279) q[3];
sx q[3];
rz(-0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4979424) q[0];
sx q[0];
rz(-0.90181667) q[0];
sx q[0];
rz(1.3051916) q[0];
rz(2.7611043) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(0.25340733) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2449269) q[0];
sx q[0];
rz(-1.9650808) q[0];
sx q[0];
rz(3.0299597) q[0];
x q[1];
rz(1.2316522) q[2];
sx q[2];
rz(-1.412743) q[2];
sx q[2];
rz(-0.16976419) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3897755) q[1];
sx q[1];
rz(-2.7174065) q[1];
sx q[1];
rz(-0.61702375) q[1];
x q[2];
rz(2.5640423) q[3];
sx q[3];
rz(-2.4452219) q[3];
sx q[3];
rz(1.3814572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5499251) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(-2.6386476) q[2];
rz(-0.89899603) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(1.1635273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1702561) q[0];
sx q[0];
rz(-1.5383056) q[0];
sx q[0];
rz(-2.8785895) q[0];
rz(0.7111711) q[1];
sx q[1];
rz(-2.053459) q[1];
sx q[1];
rz(-1.4278535) q[1];
rz(1.3118369) q[2];
sx q[2];
rz(-1.8428409) q[2];
sx q[2];
rz(0.98696282) q[2];
rz(2.4651299) q[3];
sx q[3];
rz(-2.2623895) q[3];
sx q[3];
rz(-1.1415979) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];