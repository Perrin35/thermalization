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
rz(-0.083813341) q[0];
sx q[0];
rz(-2.7901791) q[0];
sx q[0];
rz(2.0961528) q[0];
rz(0.84402973) q[1];
sx q[1];
rz(-2.5764155) q[1];
sx q[1];
rz(1.3619818) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8979826) q[0];
sx q[0];
rz(-1.9145537) q[0];
sx q[0];
rz(2.5512098) q[0];
rz(-pi) q[1];
rz(0.092081618) q[2];
sx q[2];
rz(-1.9972587) q[2];
sx q[2];
rz(1.895176) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1782383) q[1];
sx q[1];
rz(-0.46666086) q[1];
sx q[1];
rz(-0.51161029) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0718189) q[3];
sx q[3];
rz(-1.8706027) q[3];
sx q[3];
rz(0.48619798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7094946) q[2];
sx q[2];
rz(-2.7420036) q[2];
sx q[2];
rz(2.438365) q[2];
rz(2.382459) q[3];
sx q[3];
rz(-0.96733624) q[3];
sx q[3];
rz(0.67559344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21509875) q[0];
sx q[0];
rz(-0.7248942) q[0];
sx q[0];
rz(-2.3469927) q[0];
rz(-2.3130747) q[1];
sx q[1];
rz(-1.0732032) q[1];
sx q[1];
rz(0.4812831) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0006162) q[0];
sx q[0];
rz(-3.0191026) q[0];
sx q[0];
rz(0.23167892) q[0];
rz(-pi) q[1];
rz(-2.0688624) q[2];
sx q[2];
rz(-2.2553372) q[2];
sx q[2];
rz(-1.3985776) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0161288) q[1];
sx q[1];
rz(-1.4939346) q[1];
sx q[1];
rz(1.3400588) q[1];
rz(-pi) q[2];
rz(0.68010211) q[3];
sx q[3];
rz(-2.3643593) q[3];
sx q[3];
rz(-0.2118563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.98722297) q[2];
sx q[2];
rz(-0.96174806) q[2];
sx q[2];
rz(-2.6340384) q[2];
rz(0.65732035) q[3];
sx q[3];
rz(-2.1828914) q[3];
sx q[3];
rz(-1.9067732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9364612) q[0];
sx q[0];
rz(-1.5345804) q[0];
sx q[0];
rz(-0.24612799) q[0];
rz(0.71324619) q[1];
sx q[1];
rz(-1.5245707) q[1];
sx q[1];
rz(-0.86457843) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.789398) q[0];
sx q[0];
rz(-1.3653132) q[0];
sx q[0];
rz(-1.2557538) q[0];
x q[1];
rz(0.53094338) q[2];
sx q[2];
rz(-1.0304255) q[2];
sx q[2];
rz(-2.2744389) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44251014) q[1];
sx q[1];
rz(-1.1103295) q[1];
sx q[1];
rz(-0.090126474) q[1];
x q[2];
rz(-1.4186612) q[3];
sx q[3];
rz(-0.78923038) q[3];
sx q[3];
rz(-1.8735739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.41512179) q[2];
sx q[2];
rz(-0.39602009) q[2];
sx q[2];
rz(-1.2064639) q[2];
rz(0.11780277) q[3];
sx q[3];
rz(-2.466187) q[3];
sx q[3];
rz(-0.12731586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29663157) q[0];
sx q[0];
rz(-1.6806108) q[0];
sx q[0];
rz(2.6764181) q[0];
rz(-0.90452114) q[1];
sx q[1];
rz(-2.1280839) q[1];
sx q[1];
rz(-1.4314338) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44712191) q[0];
sx q[0];
rz(-1.1621124) q[0];
sx q[0];
rz(-2.7848836) q[0];
rz(1.935237) q[2];
sx q[2];
rz(-2.6195457) q[2];
sx q[2];
rz(-2.4625157) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.91421724) q[1];
sx q[1];
rz(-2.3721116) q[1];
sx q[1];
rz(0.51336594) q[1];
rz(-pi) q[2];
rz(1.9172098) q[3];
sx q[3];
rz(-0.87408057) q[3];
sx q[3];
rz(2.5908366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5005834) q[2];
sx q[2];
rz(-2.9312129) q[2];
sx q[2];
rz(-1.6453936) q[2];
rz(3.1032041) q[3];
sx q[3];
rz(-1.8206785) q[3];
sx q[3];
rz(-2.3827609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.1372304) q[0];
sx q[0];
rz(-0.69843233) q[0];
sx q[0];
rz(-1.6001562) q[0];
rz(-1.6372708) q[1];
sx q[1];
rz(-1.5613104) q[1];
sx q[1];
rz(-1.6391594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2509643) q[0];
sx q[0];
rz(-0.67442964) q[0];
sx q[0];
rz(-0.78309632) q[0];
rz(-1.4203625) q[2];
sx q[2];
rz(-1.022136) q[2];
sx q[2];
rz(2.7799431) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1790601) q[1];
sx q[1];
rz(-2.0174562) q[1];
sx q[1];
rz(-2.0455749) q[1];
x q[2];
rz(1.5030902) q[3];
sx q[3];
rz(-1.6793402) q[3];
sx q[3];
rz(-2.6084765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0043103546) q[2];
sx q[2];
rz(-0.31108019) q[2];
sx q[2];
rz(2.7610371) q[2];
rz(0.54008326) q[3];
sx q[3];
rz(-1.249908) q[3];
sx q[3];
rz(-1.165747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49749097) q[0];
sx q[0];
rz(-3.0971425) q[0];
sx q[0];
rz(-2.7569726) q[0];
rz(-0.050845536) q[1];
sx q[1];
rz(-2.6652002) q[1];
sx q[1];
rz(-1.6772038) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8576524) q[0];
sx q[0];
rz(-0.63125718) q[0];
sx q[0];
rz(0.75729482) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88846859) q[2];
sx q[2];
rz(-1.1184449) q[2];
sx q[2];
rz(2.7301228) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6932043) q[1];
sx q[1];
rz(-0.12388661) q[1];
sx q[1];
rz(2.3883567) q[1];
rz(2.4568632) q[3];
sx q[3];
rz(-1.7312418) q[3];
sx q[3];
rz(-1.1477276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9836318) q[2];
sx q[2];
rz(-2.3775358) q[2];
sx q[2];
rz(-0.91510406) q[2];
rz(-0.30465952) q[3];
sx q[3];
rz(-2.4200078) q[3];
sx q[3];
rz(1.6662395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7810998) q[0];
sx q[0];
rz(-2.2486794) q[0];
sx q[0];
rz(1.1156981) q[0];
rz(-2.5005493) q[1];
sx q[1];
rz(-1.2572925) q[1];
sx q[1];
rz(1.8152274) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4531739) q[0];
sx q[0];
rz(-0.50473112) q[0];
sx q[0];
rz(-0.27413989) q[0];
rz(2.7671779) q[2];
sx q[2];
rz(-1.7257938) q[2];
sx q[2];
rz(-2.5520476) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4600582) q[1];
sx q[1];
rz(-1.0435972) q[1];
sx q[1];
rz(2.3982226) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2828045) q[3];
sx q[3];
rz(-0.988171) q[3];
sx q[3];
rz(-1.7395541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.556584) q[2];
sx q[2];
rz(-0.98238397) q[2];
sx q[2];
rz(-1.8741685) q[2];
rz(-3.0804539) q[3];
sx q[3];
rz(-1.8784524) q[3];
sx q[3];
rz(0.88501969) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0201482) q[0];
sx q[0];
rz(-0.76524884) q[0];
sx q[0];
rz(0.52566093) q[0];
rz(1.6540182) q[1];
sx q[1];
rz(-1.1143538) q[1];
sx q[1];
rz(2.9590327) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8480685) q[0];
sx q[0];
rz(-1.7901359) q[0];
sx q[0];
rz(0.42869403) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0261954) q[2];
sx q[2];
rz(-1.8517074) q[2];
sx q[2];
rz(1.2942435) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.1614496) q[1];
sx q[1];
rz(-1.399446) q[1];
sx q[1];
rz(-2.2049987) q[1];
rz(-pi) q[2];
rz(0.28589272) q[3];
sx q[3];
rz(-2.0276311) q[3];
sx q[3];
rz(1.735812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.93099) q[2];
sx q[2];
rz(-1.8068204) q[2];
sx q[2];
rz(0.1087428) q[2];
rz(-3.0336618) q[3];
sx q[3];
rz(-2.7863672) q[3];
sx q[3];
rz(1.3817374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050960798) q[0];
sx q[0];
rz(-1.2010295) q[0];
sx q[0];
rz(-2.4884124) q[0];
rz(1.2921035) q[1];
sx q[1];
rz(-2.8958246) q[1];
sx q[1];
rz(-0.076315708) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6298593) q[0];
sx q[0];
rz(-2.7329067) q[0];
sx q[0];
rz(0.9291533) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7387975) q[2];
sx q[2];
rz(-1.6036878) q[2];
sx q[2];
rz(0.55247441) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.46901033) q[1];
sx q[1];
rz(-0.82090506) q[1];
sx q[1];
rz(1.3806848) q[1];
rz(-pi) q[2];
rz(-2.7870922) q[3];
sx q[3];
rz(-0.97149847) q[3];
sx q[3];
rz(1.3196303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.896758) q[2];
sx q[2];
rz(-1.9568169) q[2];
sx q[2];
rz(0.99311382) q[2];
rz(1.7874329) q[3];
sx q[3];
rz(-2.6024151) q[3];
sx q[3];
rz(1.7206934) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83574522) q[0];
sx q[0];
rz(-1.5033686) q[0];
sx q[0];
rz(2.848023) q[0];
rz(-2.1096032) q[1];
sx q[1];
rz(-2.3753765) q[1];
sx q[1];
rz(2.8242677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62516934) q[0];
sx q[0];
rz(-0.95808555) q[0];
sx q[0];
rz(1.4132383) q[0];
rz(-pi) q[1];
rz(2.7950049) q[2];
sx q[2];
rz(-0.58243299) q[2];
sx q[2];
rz(-2.1726441) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7723497) q[1];
sx q[1];
rz(-2.0017821) q[1];
sx q[1];
rz(1.950077) q[1];
rz(-2.9568697) q[3];
sx q[3];
rz(-2.8392042) q[3];
sx q[3];
rz(2.8222366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.62852922) q[2];
sx q[2];
rz(-0.3231914) q[2];
sx q[2];
rz(-0.31488669) q[2];
rz(-0.033673938) q[3];
sx q[3];
rz(-2.0458872) q[3];
sx q[3];
rz(-1.4382188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2122129) q[0];
sx q[0];
rz(-1.0502945) q[0];
sx q[0];
rz(2.1631277) q[0];
rz(0.20044151) q[1];
sx q[1];
rz(-2.0574175) q[1];
sx q[1];
rz(2.5331694) q[1];
rz(0.54666109) q[2];
sx q[2];
rz(-3.0849059) q[2];
sx q[2];
rz(1.0824587) q[2];
rz(1.325812) q[3];
sx q[3];
rz(-1.0272588) q[3];
sx q[3];
rz(2.5486265) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
