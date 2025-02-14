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
rz(3.0577793) q[0];
sx q[0];
rz(-0.35141355) q[0];
sx q[0];
rz(-2.0961528) q[0];
rz(0.84402973) q[1];
sx q[1];
rz(-2.5764155) q[1];
sx q[1];
rz(1.3619818) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24361006) q[0];
sx q[0];
rz(-1.9145537) q[0];
sx q[0];
rz(0.59038281) q[0];
rz(1.9988598) q[2];
sx q[2];
rz(-1.4869823) q[2];
sx q[2];
rz(0.36255896) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2842362) q[1];
sx q[1];
rz(-1.3487089) q[1];
sx q[1];
rz(-0.41389334) q[1];
rz(-pi) q[2];
rz(-0.33884873) q[3];
sx q[3];
rz(-2.04755) q[3];
sx q[3];
rz(-0.92425215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.43209806) q[2];
sx q[2];
rz(-0.39958909) q[2];
sx q[2];
rz(-0.70322767) q[2];
rz(-0.75913366) q[3];
sx q[3];
rz(-2.1742564) q[3];
sx q[3];
rz(-0.67559344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9264939) q[0];
sx q[0];
rz(-0.7248942) q[0];
sx q[0];
rz(2.3469927) q[0];
rz(-0.82851797) q[1];
sx q[1];
rz(-2.0683894) q[1];
sx q[1];
rz(-2.6603096) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1409765) q[0];
sx q[0];
rz(-0.12249002) q[0];
sx q[0];
rz(-2.9099137) q[0];
x q[1];
rz(1.0727302) q[2];
sx q[2];
rz(-2.2553372) q[2];
sx q[2];
rz(1.7430151) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.46336922) q[1];
sx q[1];
rz(-1.3407523) q[1];
sx q[1];
rz(-0.078945625) q[1];
rz(1.0167502) q[3];
sx q[3];
rz(-2.1475128) q[3];
sx q[3];
rz(-2.5050636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98722297) q[2];
sx q[2];
rz(-0.96174806) q[2];
sx q[2];
rz(0.50755429) q[2];
rz(0.65732035) q[3];
sx q[3];
rz(-2.1828914) q[3];
sx q[3];
rz(1.2348194) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2051314) q[0];
sx q[0];
rz(-1.6070123) q[0];
sx q[0];
rz(2.8954647) q[0];
rz(0.71324619) q[1];
sx q[1];
rz(-1.617022) q[1];
sx q[1];
rz(0.86457843) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35219463) q[0];
sx q[0];
rz(-1.3653132) q[0];
sx q[0];
rz(-1.8858389) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86979308) q[2];
sx q[2];
rz(-2.403069) q[2];
sx q[2];
rz(-3.125762) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1684226) q[1];
sx q[1];
rz(-1.4900786) q[1];
sx q[1];
rz(-2.0328841) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4186612) q[3];
sx q[3];
rz(-0.78923038) q[3];
sx q[3];
rz(-1.8735739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41512179) q[2];
sx q[2];
rz(-0.39602009) q[2];
sx q[2];
rz(-1.2064639) q[2];
rz(0.11780277) q[3];
sx q[3];
rz(-2.466187) q[3];
sx q[3];
rz(3.0142768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8449611) q[0];
sx q[0];
rz(-1.6806108) q[0];
sx q[0];
rz(-2.6764181) q[0];
rz(-0.90452114) q[1];
sx q[1];
rz(-1.0135087) q[1];
sx q[1];
rz(-1.7101589) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9412044) q[0];
sx q[0];
rz(-2.6058307) q[0];
sx q[0];
rz(2.2493811) q[0];
rz(2.9393497) q[2];
sx q[2];
rz(-2.0554539) q[2];
sx q[2];
rz(1.09367) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5618973) q[1];
sx q[1];
rz(-0.919678) q[1];
sx q[1];
rz(1.1267594) q[1];
rz(1.2243829) q[3];
sx q[3];
rz(-0.87408057) q[3];
sx q[3];
rz(-2.5908366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5005834) q[2];
sx q[2];
rz(-2.9312129) q[2];
sx q[2];
rz(-1.496199) q[2];
rz(-3.1032041) q[3];
sx q[3];
rz(-1.3209141) q[3];
sx q[3];
rz(0.7588318) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1372304) q[0];
sx q[0];
rz(-0.69843233) q[0];
sx q[0];
rz(-1.6001562) q[0];
rz(-1.5043219) q[1];
sx q[1];
rz(-1.5802822) q[1];
sx q[1];
rz(-1.6391594) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.980687) q[0];
sx q[0];
rz(-1.1146) q[0];
sx q[0];
rz(-0.51551694) q[0];
rz(-0.55372766) q[2];
sx q[2];
rz(-1.6990176) q[2];
sx q[2];
rz(2.011337) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0906265) q[1];
sx q[1];
rz(-2.5016666) q[1];
sx q[1];
rz(-0.76211318) q[1];
rz(-pi) q[2];
rz(-0.10879119) q[3];
sx q[3];
rz(-1.5034893) q[3];
sx q[3];
rz(-1.045026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1372823) q[2];
sx q[2];
rz(-0.31108019) q[2];
sx q[2];
rz(2.7610371) q[2];
rz(0.54008326) q[3];
sx q[3];
rz(-1.8916847) q[3];
sx q[3];
rz(-1.9758457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49749097) q[0];
sx q[0];
rz(-0.044450132) q[0];
sx q[0];
rz(-0.38462001) q[0];
rz(-3.0907471) q[1];
sx q[1];
rz(-0.47639242) q[1];
sx q[1];
rz(1.4643889) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58005262) q[0];
sx q[0];
rz(-1.1275575) q[0];
sx q[0];
rz(-1.1053941) q[0];
x q[1];
rz(0.55942499) q[2];
sx q[2];
rz(-0.96772268) q[2];
sx q[2];
rz(-0.81808264) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6932043) q[1];
sx q[1];
rz(-3.017706) q[1];
sx q[1];
rz(-0.75323592) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8910849) q[3];
sx q[3];
rz(-2.4412812) q[3];
sx q[3];
rz(0.61628646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.1579608) q[2];
sx q[2];
rz(-0.76405683) q[2];
sx q[2];
rz(0.91510406) q[2];
rz(2.8369331) q[3];
sx q[3];
rz(-0.72158486) q[3];
sx q[3];
rz(-1.6662395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3604928) q[0];
sx q[0];
rz(-0.89291328) q[0];
sx q[0];
rz(-2.0258946) q[0];
rz(2.5005493) q[1];
sx q[1];
rz(-1.2572925) q[1];
sx q[1];
rz(-1.8152274) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64102252) q[0];
sx q[0];
rz(-1.7020853) q[0];
sx q[0];
rz(-0.48878756) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37441476) q[2];
sx q[2];
rz(-1.4157989) q[2];
sx q[2];
rz(2.5520476) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5976482) q[1];
sx q[1];
rz(-0.94606384) q[1];
sx q[1];
rz(0.90170699) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.556584) q[2];
sx q[2];
rz(-2.1592087) q[2];
sx q[2];
rz(-1.8741685) q[2];
rz(-3.0804539) q[3];
sx q[3];
rz(-1.8784524) q[3];
sx q[3];
rz(-2.256573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0201482) q[0];
sx q[0];
rz(-0.76524884) q[0];
sx q[0];
rz(2.6159317) q[0];
rz(1.6540182) q[1];
sx q[1];
rz(-1.1143538) q[1];
sx q[1];
rz(2.9590327) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2935242) q[0];
sx q[0];
rz(-1.7901359) q[0];
sx q[0];
rz(2.7128986) q[0];
x q[1];
rz(2.8307249) q[2];
sx q[2];
rz(-1.1345049) q[2];
sx q[2];
rz(2.7300961) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5044587) q[1];
sx q[1];
rz(-0.65385039) q[1];
sx q[1];
rz(-1.8549396) q[1];
rz(1.0973516) q[3];
sx q[3];
rz(-1.8266738) q[3];
sx q[3];
rz(0.29395834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21060264) q[2];
sx q[2];
rz(-1.3347722) q[2];
sx q[2];
rz(3.0328499) q[2];
rz(0.10793081) q[3];
sx q[3];
rz(-2.7863672) q[3];
sx q[3];
rz(1.3817374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050960798) q[0];
sx q[0];
rz(-1.9405631) q[0];
sx q[0];
rz(-0.65318024) q[0];
rz(1.2921035) q[1];
sx q[1];
rz(-2.8958246) q[1];
sx q[1];
rz(-0.076315708) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51173333) q[0];
sx q[0];
rz(-0.40868592) q[0];
sx q[0];
rz(0.9291533) q[0];
rz(-pi) q[1];
rz(1.4027951) q[2];
sx q[2];
rz(-1.5379049) q[2];
sx q[2];
rz(-0.55247441) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19382754) q[1];
sx q[1];
rz(-0.76903957) q[1];
sx q[1];
rz(0.2001708) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35450046) q[3];
sx q[3];
rz(-0.97149847) q[3];
sx q[3];
rz(1.8219624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.896758) q[2];
sx q[2];
rz(-1.9568169) q[2];
sx q[2];
rz(-0.99311382) q[2];
rz(1.3541597) q[3];
sx q[3];
rz(-0.5391776) q[3];
sx q[3];
rz(1.7206934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83574522) q[0];
sx q[0];
rz(-1.638224) q[0];
sx q[0];
rz(2.848023) q[0];
rz(-2.1096032) q[1];
sx q[1];
rz(-2.3753765) q[1];
sx q[1];
rz(2.8242677) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2870796) q[0];
sx q[0];
rz(-1.6995158) q[0];
sx q[0];
rz(-0.61858701) q[0];
rz(-pi) q[1];
rz(-1.3506838) q[2];
sx q[2];
rz(-2.1145472) q[2];
sx q[2];
rz(-1.3771283) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.60734487) q[1];
sx q[1];
rz(-2.5754654) q[1];
sx q[1];
rz(0.67791636) q[1];
rz(-pi) q[2];
rz(-1.6280319) q[3];
sx q[3];
rz(-1.2737107) q[3];
sx q[3];
rz(0.12606584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5130634) q[2];
sx q[2];
rz(-2.8184012) q[2];
sx q[2];
rz(0.31488669) q[2];
rz(-3.1079187) q[3];
sx q[3];
rz(-2.0458872) q[3];
sx q[3];
rz(-1.7033738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2122129) q[0];
sx q[0];
rz(-2.0912981) q[0];
sx q[0];
rz(-0.97846497) q[0];
rz(2.9411511) q[1];
sx q[1];
rz(-1.0841752) q[1];
sx q[1];
rz(-0.60842327) q[1];
rz(-1.5413054) q[2];
sx q[2];
rz(-1.5223778) q[2];
sx q[2];
rz(1.6298339) q[2];
rz(-1.325812) q[3];
sx q[3];
rz(-2.1143338) q[3];
sx q[3];
rz(-0.59296617) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
