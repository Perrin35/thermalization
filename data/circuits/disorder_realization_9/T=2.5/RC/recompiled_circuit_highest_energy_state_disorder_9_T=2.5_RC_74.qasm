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
rz(2.7901791) q[0];
sx q[0];
rz(8.3793381) q[0];
rz(-5.4391556) q[1];
sx q[1];
rz(0.56517711) q[1];
sx q[1];
rz(11.204389) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7931853) q[0];
sx q[0];
rz(-0.67272609) q[0];
sx q[0];
rz(-0.5714709) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7704813) q[2];
sx q[2];
rz(-2.7058995) q[2];
sx q[2];
rz(1.0267804) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.85735649) q[1];
sx q[1];
rz(-1.7928837) q[1];
sx q[1];
rz(2.7276993) q[1];
rz(1.0697738) q[3];
sx q[3];
rz(-1.27099) q[3];
sx q[3];
rz(0.48619798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7094946) q[2];
sx q[2];
rz(-2.7420036) q[2];
sx q[2];
rz(2.438365) q[2];
rz(-0.75913366) q[3];
sx q[3];
rz(-0.96733624) q[3];
sx q[3];
rz(0.67559344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9264939) q[0];
sx q[0];
rz(-2.4166985) q[0];
sx q[0];
rz(2.3469927) q[0];
rz(-2.3130747) q[1];
sx q[1];
rz(-2.0683894) q[1];
sx q[1];
rz(2.6603096) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7672509) q[0];
sx q[0];
rz(-1.4515948) q[0];
sx q[0];
rz(-1.5425372) q[0];
x q[1];
rz(2.3929446) q[2];
sx q[2];
rz(-1.1916797) q[2];
sx q[2];
rz(2.9826255) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.46336922) q[1];
sx q[1];
rz(-1.8008403) q[1];
sx q[1];
rz(3.062647) q[1];
x q[2];
rz(1.0167502) q[3];
sx q[3];
rz(-2.1475128) q[3];
sx q[3];
rz(0.63652906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.98722297) q[2];
sx q[2];
rz(-2.1798446) q[2];
sx q[2];
rz(0.50755429) q[2];
rz(-2.4842723) q[3];
sx q[3];
rz(-2.1828914) q[3];
sx q[3];
rz(1.2348194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9364612) q[0];
sx q[0];
rz(-1.6070123) q[0];
sx q[0];
rz(-0.24612799) q[0];
rz(2.4283465) q[1];
sx q[1];
rz(-1.617022) q[1];
sx q[1];
rz(-0.86457843) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9893892) q[0];
sx q[0];
rz(-1.2626022) q[0];
sx q[0];
rz(2.9257923) q[0];
rz(0.96295348) q[2];
sx q[2];
rz(-2.0199482) q[2];
sx q[2];
rz(-2.1445865) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6990825) q[1];
sx q[1];
rz(-2.0312632) q[1];
sx q[1];
rz(0.090126474) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7229314) q[3];
sx q[3];
rz(-0.78923038) q[3];
sx q[3];
rz(1.2680188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7264709) q[2];
sx q[2];
rz(-2.7455726) q[2];
sx q[2];
rz(-1.2064639) q[2];
rz(-3.0237899) q[3];
sx q[3];
rz(-0.67540568) q[3];
sx q[3];
rz(-3.0142768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29663157) q[0];
sx q[0];
rz(-1.6806108) q[0];
sx q[0];
rz(-0.46517459) q[0];
rz(-0.90452114) q[1];
sx q[1];
rz(-2.1280839) q[1];
sx q[1];
rz(1.7101589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44712191) q[0];
sx q[0];
rz(-1.9794802) q[0];
sx q[0];
rz(2.7848836) q[0];
rz(-pi) q[1];
rz(-0.202243) q[2];
sx q[2];
rz(-1.0861388) q[2];
sx q[2];
rz(-1.09367) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5796954) q[1];
sx q[1];
rz(-0.919678) q[1];
sx q[1];
rz(-2.0148333) q[1];
x q[2];
rz(-0.38549417) q[3];
sx q[3];
rz(-0.76498308) q[3];
sx q[3];
rz(-1.0631595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64100921) q[2];
sx q[2];
rz(-0.21037978) q[2];
sx q[2];
rz(1.6453936) q[2];
rz(-0.03838852) q[3];
sx q[3];
rz(-1.3209141) q[3];
sx q[3];
rz(2.3827609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.1372304) q[0];
sx q[0];
rz(-0.69843233) q[0];
sx q[0];
rz(1.6001562) q[0];
rz(-1.5043219) q[1];
sx q[1];
rz(-1.5802822) q[1];
sx q[1];
rz(-1.6391594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.980687) q[0];
sx q[0];
rz(-1.1146) q[0];
sx q[0];
rz(-2.6260757) q[0];
x q[1];
rz(-0.2404332) q[2];
sx q[2];
rz(-2.5747262) q[2];
sx q[2];
rz(0.64449233) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0509661) q[1];
sx q[1];
rz(-0.6399261) q[1];
sx q[1];
rz(-0.76211318) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10879119) q[3];
sx q[3];
rz(-1.5034893) q[3];
sx q[3];
rz(-1.045026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.0043103546) q[2];
sx q[2];
rz(-0.31108019) q[2];
sx q[2];
rz(-2.7610371) q[2];
rz(-2.6015094) q[3];
sx q[3];
rz(-1.249908) q[3];
sx q[3];
rz(-1.165747) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49749097) q[0];
sx q[0];
rz(-0.044450132) q[0];
sx q[0];
rz(-0.38462001) q[0];
rz(0.050845536) q[1];
sx q[1];
rz(-0.47639242) q[1];
sx q[1];
rz(-1.6772038) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.56154) q[0];
sx q[0];
rz(-2.0140352) q[0];
sx q[0];
rz(1.1053941) q[0];
rz(-pi) q[1];
rz(0.55942499) q[2];
sx q[2];
rz(-2.17387) q[2];
sx q[2];
rz(-2.32351) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.44838833) q[1];
sx q[1];
rz(-0.12388661) q[1];
sx q[1];
rz(-2.3883567) q[1];
x q[2];
rz(-0.25050779) q[3];
sx q[3];
rz(-0.70031149) q[3];
sx q[3];
rz(2.5253062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.1579608) q[2];
sx q[2];
rz(-0.76405683) q[2];
sx q[2];
rz(-2.2264886) q[2];
rz(-2.8369331) q[3];
sx q[3];
rz(-0.72158486) q[3];
sx q[3];
rz(1.6662395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3604928) q[0];
sx q[0];
rz(-0.89291328) q[0];
sx q[0];
rz(1.1156981) q[0];
rz(-2.5005493) q[1];
sx q[1];
rz(-1.8843001) q[1];
sx q[1];
rz(-1.8152274) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68841877) q[0];
sx q[0];
rz(-2.6368615) q[0];
sx q[0];
rz(-0.27413989) q[0];
rz(0.37441476) q[2];
sx q[2];
rz(-1.7257938) q[2];
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
rz(-1.4600582) q[1];
sx q[1];
rz(-2.0979954) q[1];
sx q[1];
rz(-2.3982226) q[1];
x q[2];
rz(-2.3604403) q[3];
sx q[3];
rz(-2.255126) q[3];
sx q[3];
rz(0.39881166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.556584) q[2];
sx q[2];
rz(-0.98238397) q[2];
sx q[2];
rz(1.8741685) q[2];
rz(0.061138717) q[3];
sx q[3];
rz(-1.2631402) q[3];
sx q[3];
rz(2.256573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1214445) q[0];
sx q[0];
rz(-2.3763438) q[0];
sx q[0];
rz(2.6159317) q[0];
rz(-1.6540182) q[1];
sx q[1];
rz(-2.0272389) q[1];
sx q[1];
rz(2.9590327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4199969) q[0];
sx q[0];
rz(-2.6631626) q[0];
sx q[0];
rz(-0.49225505) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0261954) q[2];
sx q[2];
rz(-1.2898852) q[2];
sx q[2];
rz(1.2942435) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6371339) q[1];
sx q[1];
rz(-0.65385039) q[1];
sx q[1];
rz(1.2866531) q[1];
rz(-2.0917039) q[3];
sx q[3];
rz(-2.6080797) q[3];
sx q[3];
rz(-0.81797879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.93099) q[2];
sx q[2];
rz(-1.8068204) q[2];
sx q[2];
rz(0.1087428) q[2];
rz(0.10793081) q[3];
sx q[3];
rz(-0.35522541) q[3];
sx q[3];
rz(-1.3817374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0906319) q[0];
sx q[0];
rz(-1.9405631) q[0];
sx q[0];
rz(-0.65318024) q[0];
rz(-1.2921035) q[1];
sx q[1];
rz(-2.8958246) q[1];
sx q[1];
rz(-3.0652769) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6600403) q[0];
sx q[0];
rz(-1.8109492) q[0];
sx q[0];
rz(-1.9047407) q[0];
x q[1];
rz(1.7387975) q[2];
sx q[2];
rz(-1.5379049) q[2];
sx q[2];
rz(-2.5891182) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2321978) q[1];
sx q[1];
rz(-1.7095209) q[1];
sx q[1];
rz(0.75896427) q[1];
rz(-pi) q[2];
rz(-0.35450046) q[3];
sx q[3];
rz(-2.1700942) q[3];
sx q[3];
rz(-1.8219624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2448347) q[2];
sx q[2];
rz(-1.1847757) q[2];
sx q[2];
rz(-0.99311382) q[2];
rz(1.7874329) q[3];
sx q[3];
rz(-2.6024151) q[3];
sx q[3];
rz(1.7206934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3058474) q[0];
sx q[0];
rz(-1.5033686) q[0];
sx q[0];
rz(-2.848023) q[0];
rz(1.0319895) q[1];
sx q[1];
rz(-0.76621619) q[1];
sx q[1];
rz(-2.8242677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.246884) q[0];
sx q[0];
rz(-2.5114759) q[0];
sx q[0];
rz(2.9219841) q[0];
x q[1];
rz(0.34658771) q[2];
sx q[2];
rz(-0.58243299) q[2];
sx q[2];
rz(2.1726441) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60734487) q[1];
sx q[1];
rz(-0.56612724) q[1];
sx q[1];
rz(-0.67791636) q[1];
rz(-0.18472291) q[3];
sx q[3];
rz(-2.8392042) q[3];
sx q[3];
rz(0.31935605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5130634) q[2];
sx q[2];
rz(-0.3231914) q[2];
sx q[2];
rz(2.826706) q[2];
rz(-3.1079187) q[3];
sx q[3];
rz(-1.0957054) q[3];
sx q[3];
rz(1.7033738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92937975) q[0];
sx q[0];
rz(-2.0912981) q[0];
sx q[0];
rz(-0.97846497) q[0];
rz(-2.9411511) q[1];
sx q[1];
rz(-2.0574175) q[1];
sx q[1];
rz(2.5331694) q[1];
rz(-0.54666109) q[2];
sx q[2];
rz(-0.056686747) q[2];
sx q[2];
rz(-2.0591339) q[2];
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
