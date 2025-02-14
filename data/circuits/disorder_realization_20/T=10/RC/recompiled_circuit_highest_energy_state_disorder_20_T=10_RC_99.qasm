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
rz(-2.1751997) q[0];
sx q[0];
rz(-1.3858495) q[0];
sx q[0];
rz(0.59544271) q[0];
rz(1.2915986) q[1];
sx q[1];
rz(-2.5505677) q[1];
sx q[1];
rz(-0.12330595) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3759252) q[0];
sx q[0];
rz(-1.0320837) q[0];
sx q[0];
rz(-2.6508191) q[0];
x q[1];
rz(-0.78211981) q[2];
sx q[2];
rz(-2.189866) q[2];
sx q[2];
rz(-1.9591046) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1733381) q[1];
sx q[1];
rz(-2.9072793) q[1];
sx q[1];
rz(1.987275) q[1];
x q[2];
rz(-2.7601065) q[3];
sx q[3];
rz(-2.7499928) q[3];
sx q[3];
rz(0.42575437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3143602) q[2];
sx q[2];
rz(-2.0388956) q[2];
sx q[2];
rz(-3.0296791) q[2];
rz(1.3917475) q[3];
sx q[3];
rz(-1.1575969) q[3];
sx q[3];
rz(-0.30645034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4157442) q[0];
sx q[0];
rz(-0.84472504) q[0];
sx q[0];
rz(2.9916905) q[0];
rz(0.28597486) q[1];
sx q[1];
rz(-1.7466702) q[1];
sx q[1];
rz(-1.1289977) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1492831) q[0];
sx q[0];
rz(-0.2899) q[0];
sx q[0];
rz(1.6768181) q[0];
x q[1];
rz(-0.01458077) q[2];
sx q[2];
rz(-1.6304071) q[2];
sx q[2];
rz(2.4042442) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.93176881) q[1];
sx q[1];
rz(-1.9289513) q[1];
sx q[1];
rz(-1.3937373) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5516594) q[3];
sx q[3];
rz(-2.1141756) q[3];
sx q[3];
rz(2.6389299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4631606) q[2];
sx q[2];
rz(-2.2191935) q[2];
sx q[2];
rz(-1.9909667) q[2];
rz(-1.9567418) q[3];
sx q[3];
rz(-1.7246282) q[3];
sx q[3];
rz(-3.1209893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-2.6612369) q[0];
sx q[0];
rz(-0.24599563) q[0];
sx q[0];
rz(-2.7599957) q[0];
rz(1.8474139) q[1];
sx q[1];
rz(-1.1062063) q[1];
sx q[1];
rz(-0.19827422) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0906774) q[0];
sx q[0];
rz(-0.28813513) q[0];
sx q[0];
rz(-0.91403065) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2402693) q[2];
sx q[2];
rz(-1.5076037) q[2];
sx q[2];
rz(-0.92869273) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1041996) q[1];
sx q[1];
rz(-1.8036808) q[1];
sx q[1];
rz(1.3523471) q[1];
rz(-pi) q[2];
rz(1.4573077) q[3];
sx q[3];
rz(-2.9656565) q[3];
sx q[3];
rz(-2.4324377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.76564378) q[2];
sx q[2];
rz(-0.15268923) q[2];
sx q[2];
rz(-2.622733) q[2];
rz(1.8910003) q[3];
sx q[3];
rz(-2.0542681) q[3];
sx q[3];
rz(-0.58108228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64312235) q[0];
sx q[0];
rz(-0.20007087) q[0];
sx q[0];
rz(0.44922391) q[0];
rz(-1.6953702) q[1];
sx q[1];
rz(-2.4999373) q[1];
sx q[1];
rz(-0.62072388) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0064471) q[0];
sx q[0];
rz(-2.2250611) q[0];
sx q[0];
rz(-1.3885137) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16731842) q[2];
sx q[2];
rz(-0.64019055) q[2];
sx q[2];
rz(0.13969914) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6606969) q[1];
sx q[1];
rz(-2.6891853) q[1];
sx q[1];
rz(0.47641944) q[1];
x q[2];
rz(-1.4820547) q[3];
sx q[3];
rz(-0.54026287) q[3];
sx q[3];
rz(-2.1335354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0351403) q[2];
sx q[2];
rz(-1.2713212) q[2];
sx q[2];
rz(-2.980496) q[2];
rz(0.39786878) q[3];
sx q[3];
rz(-1.4784003) q[3];
sx q[3];
rz(-0.14959344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.38254) q[0];
sx q[0];
rz(-1.362514) q[0];
sx q[0];
rz(0.25704849) q[0];
rz(-2.2611179) q[1];
sx q[1];
rz(-2.5011261) q[1];
sx q[1];
rz(-1.8880728) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2855447) q[0];
sx q[0];
rz(-0.040321983) q[0];
sx q[0];
rz(-2.0364385) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3257371) q[2];
sx q[2];
rz(-2.365621) q[2];
sx q[2];
rz(-1.6505591) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67318007) q[1];
sx q[1];
rz(-2.367503) q[1];
sx q[1];
rz(1.5062529) q[1];
rz(-pi) q[2];
rz(1.569031) q[3];
sx q[3];
rz(-1.9179419) q[3];
sx q[3];
rz(2.7070759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4132061) q[2];
sx q[2];
rz(-1.2089968) q[2];
sx q[2];
rz(-1.2665292) q[2];
rz(2.5201216) q[3];
sx q[3];
rz(-0.60756835) q[3];
sx q[3];
rz(2.6004041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40015873) q[0];
sx q[0];
rz(-2.3024237) q[0];
sx q[0];
rz(-0.025064502) q[0];
rz(-1.9301682) q[1];
sx q[1];
rz(-1.8823267) q[1];
sx q[1];
rz(2.8866344) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099440558) q[0];
sx q[0];
rz(-1.6201311) q[0];
sx q[0];
rz(-1.5758324) q[0];
rz(-pi) q[1];
rz(0.58664257) q[2];
sx q[2];
rz(-1.1544747) q[2];
sx q[2];
rz(-3.0468009) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7073696) q[1];
sx q[1];
rz(-2.1277782) q[1];
sx q[1];
rz(-2.2051478) q[1];
rz(1.2473769) q[3];
sx q[3];
rz(-2.5445523) q[3];
sx q[3];
rz(0.14123329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.24594626) q[2];
sx q[2];
rz(-2.0993555) q[2];
sx q[2];
rz(-0.45905217) q[2];
rz(1.4970655) q[3];
sx q[3];
rz(-0.11681695) q[3];
sx q[3];
rz(-0.12115255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6727305) q[0];
sx q[0];
rz(-0.48968306) q[0];
sx q[0];
rz(0.086061867) q[0];
rz(2.22279) q[1];
sx q[1];
rz(-1.1163534) q[1];
sx q[1];
rz(1.930621) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0491517) q[0];
sx q[0];
rz(-1.3685797) q[0];
sx q[0];
rz(-0.33009712) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1384905) q[2];
sx q[2];
rz(-1.1922115) q[2];
sx q[2];
rz(2.1673983) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1001491) q[1];
sx q[1];
rz(-0.38310928) q[1];
sx q[1];
rz(-0.99975296) q[1];
rz(-0.36465148) q[3];
sx q[3];
rz(-1.9130008) q[3];
sx q[3];
rz(-2.3221644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.81898895) q[2];
sx q[2];
rz(-0.49590597) q[2];
sx q[2];
rz(-0.71713478) q[2];
rz(2.1142193) q[3];
sx q[3];
rz(-1.4077978) q[3];
sx q[3];
rz(0.43924847) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9913919) q[0];
sx q[0];
rz(-2.1450277) q[0];
sx q[0];
rz(3.126934) q[0];
rz(2.390059) q[1];
sx q[1];
rz(-1.9207759) q[1];
sx q[1];
rz(1.6709447) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96318564) q[0];
sx q[0];
rz(-1.3826332) q[0];
sx q[0];
rz(1.7117731) q[0];
x q[1];
rz(-2.884955) q[2];
sx q[2];
rz(-1.4932639) q[2];
sx q[2];
rz(-0.77285779) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4336509) q[1];
sx q[1];
rz(-1.6835064) q[1];
sx q[1];
rz(-1.1658843) q[1];
rz(-1.2257158) q[3];
sx q[3];
rz(-1.5672641) q[3];
sx q[3];
rz(0.47560233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8810001) q[2];
sx q[2];
rz(-1.1129817) q[2];
sx q[2];
rz(2.8680958) q[2];
rz(-1.986844) q[3];
sx q[3];
rz(-2.4879849) q[3];
sx q[3];
rz(2.4912513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038789373) q[0];
sx q[0];
rz(-1.1253072) q[0];
sx q[0];
rz(2.0661085) q[0];
rz(-0.12818809) q[1];
sx q[1];
rz(-0.85103858) q[1];
sx q[1];
rz(2.7959965) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81886473) q[0];
sx q[0];
rz(-2.4562307) q[0];
sx q[0];
rz(1.3166053) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32175154) q[2];
sx q[2];
rz(-2.1158127) q[2];
sx q[2];
rz(0.36361748) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4491475) q[1];
sx q[1];
rz(-1.0303823) q[1];
sx q[1];
rz(-1.497333) q[1];
x q[2];
rz(0.2407705) q[3];
sx q[3];
rz(-2.5298389) q[3];
sx q[3];
rz(1.4021378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0547611) q[2];
sx q[2];
rz(-1.0826449) q[2];
sx q[2];
rz(1.5209939) q[2];
rz(1.7704376) q[3];
sx q[3];
rz(-0.75826472) q[3];
sx q[3];
rz(2.7267406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3138251) q[0];
sx q[0];
rz(-0.96243745) q[0];
sx q[0];
rz(-2.9058822) q[0];
rz(-0.57304263) q[1];
sx q[1];
rz(-2.133281) q[1];
sx q[1];
rz(-1.7726353) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1397676) q[0];
sx q[0];
rz(-0.84472504) q[0];
sx q[0];
rz(0.98066179) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1718352) q[2];
sx q[2];
rz(-0.98669334) q[2];
sx q[2];
rz(-2.7764708) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6073709) q[1];
sx q[1];
rz(-2.1335619) q[1];
sx q[1];
rz(2.9396179) q[1];
rz(-pi) q[2];
x q[2];
rz(1.664613) q[3];
sx q[3];
rz(-1.929885) q[3];
sx q[3];
rz(0.94993587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3706751) q[2];
sx q[2];
rz(-1.3266027) q[2];
sx q[2];
rz(3.0885922) q[2];
rz(0.75183374) q[3];
sx q[3];
rz(-1.0337318) q[3];
sx q[3];
rz(-2.2161765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62405217) q[0];
sx q[0];
rz(-1.0160099) q[0];
sx q[0];
rz(2.1852063) q[0];
rz(-1.3507631) q[1];
sx q[1];
rz(-2.7443934) q[1];
sx q[1];
rz(2.9846356) q[1];
rz(0.15539697) q[2];
sx q[2];
rz(-2.4383895) q[2];
sx q[2];
rz(-0.40021605) q[2];
rz(-2.5637473) q[3];
sx q[3];
rz(-2.2457849) q[3];
sx q[3];
rz(-0.59413322) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
