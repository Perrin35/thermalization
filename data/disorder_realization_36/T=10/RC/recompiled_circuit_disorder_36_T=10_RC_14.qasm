OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49178034) q[0];
sx q[0];
rz(3.427504) q[0];
sx q[0];
rz(8.9094845) q[0];
rz(1.3442858) q[1];
sx q[1];
rz(-2.9872515) q[1];
sx q[1];
rz(0.57758346) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87969765) q[0];
sx q[0];
rz(-1.4059773) q[0];
sx q[0];
rz(-0.66189712) q[0];
rz(-pi) q[1];
rz(-1.0628113) q[2];
sx q[2];
rz(-1.0928109) q[2];
sx q[2];
rz(-2.9302772) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88749332) q[1];
sx q[1];
rz(-1.8389529) q[1];
sx q[1];
rz(-1.1521794) q[1];
rz(2.203381) q[3];
sx q[3];
rz(-2.1198366) q[3];
sx q[3];
rz(0.051018056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.704533) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(2.4543767) q[2];
rz(-1.0152738) q[3];
sx q[3];
rz(-1.7679368) q[3];
sx q[3];
rz(-3.0190873) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17094831) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(1.2600391) q[0];
rz(-2.1353703) q[1];
sx q[1];
rz(-2.1496014) q[1];
sx q[1];
rz(-0.84567436) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0613522) q[0];
sx q[0];
rz(-2.4363359) q[0];
sx q[0];
rz(2.4387226) q[0];
rz(-0.13375608) q[2];
sx q[2];
rz(-1.6998561) q[2];
sx q[2];
rz(-1.2226906) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0963124) q[1];
sx q[1];
rz(-1.7573866) q[1];
sx q[1];
rz(-0.50218302) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91652292) q[3];
sx q[3];
rz(-0.73361165) q[3];
sx q[3];
rz(0.77378002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.78850293) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(-2.6611924) q[2];
rz(-1.3530312) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95124328) q[0];
sx q[0];
rz(-0.23878637) q[0];
sx q[0];
rz(-0.7730661) q[0];
rz(0.13126016) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(2.0551596) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.438293) q[0];
sx q[0];
rz(-0.46143954) q[0];
sx q[0];
rz(-1.5745387) q[0];
rz(1.418872) q[2];
sx q[2];
rz(-0.65290367) q[2];
sx q[2];
rz(-2.8002847) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8311027) q[1];
sx q[1];
rz(-1.576014) q[1];
sx q[1];
rz(-0.28880854) q[1];
rz(-pi) q[2];
rz(-1.9757189) q[3];
sx q[3];
rz(-1.7305264) q[3];
sx q[3];
rz(0.83042849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0125668) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(1.770299) q[2];
rz(0.38315547) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(0.80254054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6894158) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(3.0932328) q[0];
rz(0.16391779) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(1.6960467) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7955129) q[0];
sx q[0];
rz(-0.68329408) q[0];
sx q[0];
rz(0.28738316) q[0];
x q[1];
rz(-0.21638685) q[2];
sx q[2];
rz(-1.6587703) q[2];
sx q[2];
rz(-0.83425922) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1452892) q[1];
sx q[1];
rz(-0.22876303) q[1];
sx q[1];
rz(0.28782515) q[1];
rz(-2.0479855) q[3];
sx q[3];
rz(-2.5683937) q[3];
sx q[3];
rz(1.5447865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0754898) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(-2.1172822) q[2];
rz(1.5284437) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(-2.9083692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7016474) q[0];
sx q[0];
rz(-2.3174536) q[0];
sx q[0];
rz(-1.2874999) q[0];
rz(0.31907407) q[1];
sx q[1];
rz(-1.5417475) q[1];
sx q[1];
rz(-2.2873926) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.655495) q[0];
sx q[0];
rz(-0.71394074) q[0];
sx q[0];
rz(1.5833202) q[0];
rz(-2.1335667) q[2];
sx q[2];
rz(-2.3587583) q[2];
sx q[2];
rz(-0.40751878) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2687208) q[1];
sx q[1];
rz(-0.99209058) q[1];
sx q[1];
rz(-0.89923664) q[1];
x q[2];
rz(-2.5913127) q[3];
sx q[3];
rz(-1.0770505) q[3];
sx q[3];
rz(2.9328336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1896818) q[2];
sx q[2];
rz(-2.5482735) q[2];
sx q[2];
rz(-0.577315) q[2];
rz(-2.632085) q[3];
sx q[3];
rz(-2.7039492) q[3];
sx q[3];
rz(-0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1381056) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(3.0694718) q[0];
rz(2.0347319) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(-0.12621005) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.431625) q[0];
sx q[0];
rz(-2.9266848) q[0];
sx q[0];
rz(2.9931195) q[0];
x q[1];
rz(1.1926786) q[2];
sx q[2];
rz(-0.83231229) q[2];
sx q[2];
rz(-2.9044915) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.63657657) q[1];
sx q[1];
rz(-0.86537433) q[1];
sx q[1];
rz(0.24800639) q[1];
x q[2];
rz(-1.24228) q[3];
sx q[3];
rz(-1.7528755) q[3];
sx q[3];
rz(-1.8329221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2386027) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(0.77511707) q[2];
rz(0.827968) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(-1.1221788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5489952) q[0];
sx q[0];
rz(-0.42625517) q[0];
sx q[0];
rz(-0.098408498) q[0];
rz(1.9495643) q[1];
sx q[1];
rz(-1.3339309) q[1];
sx q[1];
rz(0.55955204) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7647117) q[0];
sx q[0];
rz(-1.6983713) q[0];
sx q[0];
rz(-0.60140951) q[0];
rz(-pi) q[1];
rz(-2.8801708) q[2];
sx q[2];
rz(-1.520642) q[2];
sx q[2];
rz(-2.9068771) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4899788) q[1];
sx q[1];
rz(-1.4414806) q[1];
sx q[1];
rz(-3.053385) q[1];
x q[2];
rz(0.63013245) q[3];
sx q[3];
rz(-1.7297941) q[3];
sx q[3];
rz(-0.4160479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6825535) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(2.7977978) q[2];
rz(2.5750459) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(-2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7664117) q[0];
sx q[0];
rz(-1.8090929) q[0];
sx q[0];
rz(2.4108316) q[0];
rz(2.9991951) q[1];
sx q[1];
rz(-1.8715033) q[1];
sx q[1];
rz(-0.87160814) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5985142) q[0];
sx q[0];
rz(-0.075965479) q[0];
sx q[0];
rz(-0.11673467) q[0];
rz(-pi) q[1];
rz(-0.40558221) q[2];
sx q[2];
rz(-1.7103346) q[2];
sx q[2];
rz(1.9987193) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4910482) q[1];
sx q[1];
rz(-0.77495134) q[1];
sx q[1];
rz(-2.6130555) q[1];
x q[2];
rz(-1.3685162) q[3];
sx q[3];
rz(-2.5296122) q[3];
sx q[3];
rz(-0.66093854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4153851) q[2];
sx q[2];
rz(-2.0724847) q[2];
sx q[2];
rz(1.139572) q[2];
rz(1.6428044) q[3];
sx q[3];
rz(-2.747624) q[3];
sx q[3];
rz(0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20294872) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(1.9198445) q[0];
rz(2.9755759) q[1];
sx q[1];
rz(-1.8211726) q[1];
sx q[1];
rz(-1.6171914) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7302007) q[0];
sx q[0];
rz(-3.054266) q[0];
sx q[0];
rz(-1.9090396) q[0];
rz(-2.6644601) q[2];
sx q[2];
rz(-1.7413365) q[2];
sx q[2];
rz(1.0605304) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0601378) q[1];
sx q[1];
rz(-0.36768915) q[1];
sx q[1];
rz(2.0245488) q[1];
rz(-pi) q[2];
rz(-1.2576305) q[3];
sx q[3];
rz(-0.95542613) q[3];
sx q[3];
rz(2.9477011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.021585492) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(-2.7900556) q[2];
rz(1.0567788) q[3];
sx q[3];
rz(-0.52962279) q[3];
sx q[3];
rz(2.3969011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64365023) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(1.836401) q[0];
rz(0.38048831) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(-0.25340733) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2828335) q[0];
sx q[0];
rz(-1.6738322) q[0];
sx q[0];
rz(1.9673002) q[0];
x q[1];
rz(-1.1240187) q[2];
sx q[2];
rz(-0.3728711) q[2];
sx q[2];
rz(2.1602221) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0903783) q[1];
sx q[1];
rz(-1.9131294) q[1];
sx q[1];
rz(1.8263837) q[1];
x q[2];
rz(-0.57755034) q[3];
sx q[3];
rz(-2.4452219) q[3];
sx q[3];
rz(1.3814572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5499251) q[2];
sx q[2];
rz(-2.2558236) q[2];
sx q[2];
rz(-0.5029451) q[2];
rz(2.2425966) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(1.1635273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1702561) q[0];
sx q[0];
rz(-1.6032871) q[0];
sx q[0];
rz(0.26300318) q[0];
rz(-2.4304216) q[1];
sx q[1];
rz(-2.053459) q[1];
sx q[1];
rz(-1.4278535) q[1];
rz(-1.8297557) q[2];
sx q[2];
rz(-1.8428409) q[2];
sx q[2];
rz(0.98696282) q[2];
rz(-0.75541227) q[3];
sx q[3];
rz(-2.073954) q[3];
sx q[3];
rz(3.0975773) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];