OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.98439944) q[0];
sx q[0];
rz(-1.1097044) q[0];
sx q[0];
rz(1.00534) q[0];
rz(-2.4069064) q[1];
sx q[1];
rz(-0.35597304) q[1];
sx q[1];
rz(0.89207831) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9514314) q[0];
sx q[0];
rz(-2.3888458) q[0];
sx q[0];
rz(-2.997082) q[0];
rz(1.35688) q[2];
sx q[2];
rz(-0.61442539) q[2];
sx q[2];
rz(0.015903552) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0301117) q[1];
sx q[1];
rz(-1.4731543) q[1];
sx q[1];
rz(-2.6565115) q[1];
rz(-pi) q[2];
rz(-1.0451116) q[3];
sx q[3];
rz(-2.1108004) q[3];
sx q[3];
rz(1.2837376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72508183) q[2];
sx q[2];
rz(-1.9981013) q[2];
sx q[2];
rz(1.7007281) q[2];
rz(0.95669389) q[3];
sx q[3];
rz(-1.1172833) q[3];
sx q[3];
rz(2.8982437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6956534) q[0];
sx q[0];
rz(-2.74701) q[0];
sx q[0];
rz(-1.12895) q[0];
rz(2.8961862) q[1];
sx q[1];
rz(-1.8281728) q[1];
sx q[1];
rz(-0.66171563) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11479325) q[0];
sx q[0];
rz(-2.4519992) q[0];
sx q[0];
rz(-2.6845776) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7595046) q[2];
sx q[2];
rz(-1.7799062) q[2];
sx q[2];
rz(1.2651625) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7364517) q[1];
sx q[1];
rz(-1.5390522) q[1];
sx q[1];
rz(-2.870138) q[1];
rz(-pi) q[2];
rz(-1.2564959) q[3];
sx q[3];
rz(-2.8366025) q[3];
sx q[3];
rz(-1.5105607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5439862) q[2];
sx q[2];
rz(-2.643955) q[2];
sx q[2];
rz(1.0428766) q[2];
rz(-1.6889702) q[3];
sx q[3];
rz(-0.85113227) q[3];
sx q[3];
rz(-2.0378621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76688981) q[0];
sx q[0];
rz(-2.4536528) q[0];
sx q[0];
rz(-1.8324628) q[0];
rz(-0.4492999) q[1];
sx q[1];
rz(-0.6895014) q[1];
sx q[1];
rz(-1.7151394) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5912522) q[0];
sx q[0];
rz(-0.36446291) q[0];
sx q[0];
rz(-2.3191207) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9932991) q[2];
sx q[2];
rz(-1.587217) q[2];
sx q[2];
rz(-2.8016479) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0589185) q[1];
sx q[1];
rz(-2.2737496) q[1];
sx q[1];
rz(1.6723796) q[1];
rz(-pi) q[2];
rz(-1.8168114) q[3];
sx q[3];
rz(-1.6520471) q[3];
sx q[3];
rz(-1.2556835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7766777) q[2];
sx q[2];
rz(-1.9675576) q[2];
sx q[2];
rz(3.0752227) q[2];
rz(-0.55473173) q[3];
sx q[3];
rz(-1.6603671) q[3];
sx q[3];
rz(1.5167282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1944815) q[0];
sx q[0];
rz(-0.28488657) q[0];
sx q[0];
rz(-2.3696005) q[0];
rz(-2.4783065) q[1];
sx q[1];
rz(-2.3978077) q[1];
sx q[1];
rz(-3.0139121) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6948815) q[0];
sx q[0];
rz(-0.90411964) q[0];
sx q[0];
rz(-1.5645909) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.588932) q[2];
sx q[2];
rz(-1.2465806) q[2];
sx q[2];
rz(0.49190258) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.046317958) q[1];
sx q[1];
rz(-1.1409586) q[1];
sx q[1];
rz(2.6916987) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7112682) q[3];
sx q[3];
rz(-1.7146304) q[3];
sx q[3];
rz(-2.6324684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.53106236) q[2];
sx q[2];
rz(-2.1805111) q[2];
sx q[2];
rz(2.5430211) q[2];
rz(0.5851723) q[3];
sx q[3];
rz(-1.7681311) q[3];
sx q[3];
rz(3.0857871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7912306) q[0];
sx q[0];
rz(-1.6329916) q[0];
sx q[0];
rz(0.97306657) q[0];
rz(2.9452501) q[1];
sx q[1];
rz(-1.1390431) q[1];
sx q[1];
rz(1.6729209) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.136841) q[0];
sx q[0];
rz(-0.84293619) q[0];
sx q[0];
rz(2.1086295) q[0];
rz(0.25428793) q[2];
sx q[2];
rz(-1.0788222) q[2];
sx q[2];
rz(-0.51153431) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.10188516) q[1];
sx q[1];
rz(-0.3488003) q[1];
sx q[1];
rz(-1.0451911) q[1];
rz(-1.0319388) q[3];
sx q[3];
rz(-1.5232289) q[3];
sx q[3];
rz(-0.69456929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68957442) q[2];
sx q[2];
rz(-1.8147899) q[2];
sx q[2];
rz(-2.9791974) q[2];
rz(1.1299805) q[3];
sx q[3];
rz(-2.2584848) q[3];
sx q[3];
rz(-2.0067154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3910386) q[0];
sx q[0];
rz(-0.51893187) q[0];
sx q[0];
rz(3.0911875) q[0];
rz(2.9734036) q[1];
sx q[1];
rz(-1.5805809) q[1];
sx q[1];
rz(-2.2680297) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45437434) q[0];
sx q[0];
rz(-0.66965398) q[0];
sx q[0];
rz(-2.7424916) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8040405) q[2];
sx q[2];
rz(-2.9961259) q[2];
sx q[2];
rz(1.625753) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7616523) q[1];
sx q[1];
rz(-2.4177175) q[1];
sx q[1];
rz(0.70751247) q[1];
rz(-2.843552) q[3];
sx q[3];
rz(-1.7665157) q[3];
sx q[3];
rz(0.51575553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0198589) q[2];
sx q[2];
rz(-0.21616082) q[2];
sx q[2];
rz(2.7609694) q[2];
rz(0.56626433) q[3];
sx q[3];
rz(-1.4004613) q[3];
sx q[3];
rz(0.80999058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64067632) q[0];
sx q[0];
rz(-0.2114507) q[0];
sx q[0];
rz(-0.40263116) q[0];
rz(1.7598049) q[1];
sx q[1];
rz(-0.80843061) q[1];
sx q[1];
rz(3.0852539) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0005971) q[0];
sx q[0];
rz(-2.6133446) q[0];
sx q[0];
rz(-0.88185593) q[0];
rz(-pi) q[1];
rz(-1.1001141) q[2];
sx q[2];
rz(-1.6073445) q[2];
sx q[2];
rz(1.4135897) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.90625396) q[1];
sx q[1];
rz(-1.5829931) q[1];
sx q[1];
rz(-2.1349299) q[1];
x q[2];
rz(-2.5272156) q[3];
sx q[3];
rz(-0.77765948) q[3];
sx q[3];
rz(-1.3953502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.23400433) q[2];
sx q[2];
rz(-1.3434709) q[2];
sx q[2];
rz(2.6410356) q[2];
rz(0.060700011) q[3];
sx q[3];
rz(-1.7874291) q[3];
sx q[3];
rz(2.6589656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57182264) q[0];
sx q[0];
rz(-1.3615384) q[0];
sx q[0];
rz(-1.3609591) q[0];
rz(-2.3979777) q[1];
sx q[1];
rz(-1.7230325) q[1];
sx q[1];
rz(-0.13882151) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9567989) q[0];
sx q[0];
rz(-2.0763872) q[0];
sx q[0];
rz(-0.57539815) q[0];
x q[1];
rz(-2.1651504) q[2];
sx q[2];
rz(-1.9658372) q[2];
sx q[2];
rz(0.76752418) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3070345) q[1];
sx q[1];
rz(-1.2456248) q[1];
sx q[1];
rz(-2.9284655) q[1];
rz(-1.9622063) q[3];
sx q[3];
rz(-1.3782805) q[3];
sx q[3];
rz(-2.8716536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6737785) q[2];
sx q[2];
rz(-2.0970586) q[2];
sx q[2];
rz(-2.4274801) q[2];
rz(-0.95947391) q[3];
sx q[3];
rz(-1.6474612) q[3];
sx q[3];
rz(-0.97904557) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4176843) q[0];
sx q[0];
rz(-2.4651616) q[0];
sx q[0];
rz(-0.12829256) q[0];
rz(0.3872321) q[1];
sx q[1];
rz(-2.5435244) q[1];
sx q[1];
rz(1.8181575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4856199) q[0];
sx q[0];
rz(-1.4733629) q[0];
sx q[0];
rz(1.3809588) q[0];
rz(0.59201957) q[2];
sx q[2];
rz(-0.97416234) q[2];
sx q[2];
rz(-0.48812619) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1897498) q[1];
sx q[1];
rz(-1.7147168) q[1];
sx q[1];
rz(3.1274432) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3714482) q[3];
sx q[3];
rz(-2.2888765) q[3];
sx q[3];
rz(1.1465286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6748176) q[2];
sx q[2];
rz(-1.520949) q[2];
sx q[2];
rz(2.9368994) q[2];
rz(-1.8681059) q[3];
sx q[3];
rz(-0.73015648) q[3];
sx q[3];
rz(-1.5654806) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8650763) q[0];
sx q[0];
rz(-0.62224329) q[0];
sx q[0];
rz(-2.9283071) q[0];
rz(-0.22008303) q[1];
sx q[1];
rz(-1.8890231) q[1];
sx q[1];
rz(1.782104) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8934568) q[0];
sx q[0];
rz(-1.5611959) q[0];
sx q[0];
rz(1.5548156) q[0];
x q[1];
rz(-2.67138) q[2];
sx q[2];
rz(-0.5222975) q[2];
sx q[2];
rz(-0.14005113) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1759739) q[1];
sx q[1];
rz(-0.63116628) q[1];
sx q[1];
rz(-2.1793773) q[1];
rz(-2.1207934) q[3];
sx q[3];
rz(-2.3925943) q[3];
sx q[3];
rz(-1.4326381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.517211) q[2];
sx q[2];
rz(-1.5456079) q[2];
sx q[2];
rz(-0.33779302) q[2];
rz(1.7289915) q[3];
sx q[3];
rz(-1.1510886) q[3];
sx q[3];
rz(2.3201449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74078858) q[0];
sx q[0];
rz(-2.3176226) q[0];
sx q[0];
rz(1.2299706) q[0];
rz(-2.7578655) q[1];
sx q[1];
rz(-1.4839254) q[1];
sx q[1];
rz(0.76212777) q[1];
rz(1.9202283) q[2];
sx q[2];
rz(-1.952662) q[2];
sx q[2];
rz(-1.9145415) q[2];
rz(1.9894133) q[3];
sx q[3];
rz(-0.55057303) q[3];
sx q[3];
rz(0.84921992) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
