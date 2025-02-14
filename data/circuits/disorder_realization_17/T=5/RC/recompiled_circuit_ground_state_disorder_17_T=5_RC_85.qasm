OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.021304) q[0];
sx q[0];
rz(2.6074183) q[0];
sx q[0];
rz(8.3745126) q[0];
rz(1.8771111) q[1];
sx q[1];
rz(-0.85327947) q[1];
sx q[1];
rz(-2.5929911) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48196822) q[0];
sx q[0];
rz(-2.7197356) q[0];
sx q[0];
rz(0.78849383) q[0];
rz(-pi) q[1];
rz(1.5496375) q[2];
sx q[2];
rz(-0.97351626) q[2];
sx q[2];
rz(-2.6279272) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9160936) q[1];
sx q[1];
rz(-2.0830983) q[1];
sx q[1];
rz(2.6698547) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24031432) q[3];
sx q[3];
rz(-1.8776456) q[3];
sx q[3];
rz(3.0003529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41142622) q[2];
sx q[2];
rz(-0.41894087) q[2];
sx q[2];
rz(2.1298998) q[2];
rz(2.2569979) q[3];
sx q[3];
rz(-1.9832289) q[3];
sx q[3];
rz(-1.8959034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0678134) q[0];
sx q[0];
rz(-0.99869204) q[0];
sx q[0];
rz(2.3538537) q[0];
rz(0.9785606) q[1];
sx q[1];
rz(-1.9906882) q[1];
sx q[1];
rz(1.2501134) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7594371) q[0];
sx q[0];
rz(-0.50001745) q[0];
sx q[0];
rz(-1.3445205) q[0];
rz(-pi) q[1];
rz(2.1627126) q[2];
sx q[2];
rz(-0.47431163) q[2];
sx q[2];
rz(1.7537376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.77630328) q[1];
sx q[1];
rz(-0.76500101) q[1];
sx q[1];
rz(1.6056152) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7911508) q[3];
sx q[3];
rz(-1.7856132) q[3];
sx q[3];
rz(-0.99220105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0222187) q[2];
sx q[2];
rz(-1.4639857) q[2];
sx q[2];
rz(3.019943) q[2];
rz(2.4308448) q[3];
sx q[3];
rz(-0.23356479) q[3];
sx q[3];
rz(0.60230437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0372593) q[0];
sx q[0];
rz(-2.4132044) q[0];
sx q[0];
rz(2.4816568) q[0];
rz(2.624699) q[1];
sx q[1];
rz(-0.70534244) q[1];
sx q[1];
rz(-1.0924115) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7832062) q[0];
sx q[0];
rz(-0.45166812) q[0];
sx q[0];
rz(2.2048401) q[0];
x q[1];
rz(0.87320341) q[2];
sx q[2];
rz(-1.3371468) q[2];
sx q[2];
rz(0.2327118) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1221083) q[1];
sx q[1];
rz(-2.1017764) q[1];
sx q[1];
rz(0.89096689) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64869439) q[3];
sx q[3];
rz(-1.3645384) q[3];
sx q[3];
rz(-2.1163396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2163781) q[2];
sx q[2];
rz(-0.34319147) q[2];
sx q[2];
rz(-1.7197616) q[2];
rz(0.4839932) q[3];
sx q[3];
rz(-1.657594) q[3];
sx q[3];
rz(-1.7234195) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5928818) q[0];
sx q[0];
rz(-3.0573248) q[0];
sx q[0];
rz(-1.3053869) q[0];
rz(0.0072172324) q[1];
sx q[1];
rz(-2.9199298) q[1];
sx q[1];
rz(0.73297393) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3586377) q[0];
sx q[0];
rz(-1.452139) q[0];
sx q[0];
rz(1.7041185) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4568366) q[2];
sx q[2];
rz(-1.2780398) q[2];
sx q[2];
rz(1.2639015) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8165255) q[1];
sx q[1];
rz(-0.96615929) q[1];
sx q[1];
rz(-2.6139392) q[1];
x q[2];
rz(-2.4184138) q[3];
sx q[3];
rz(-1.6537602) q[3];
sx q[3];
rz(1.3845598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8492154) q[2];
sx q[2];
rz(-1.7376309) q[2];
sx q[2];
rz(2.2461069) q[2];
rz(2.605947) q[3];
sx q[3];
rz(-1.5786542) q[3];
sx q[3];
rz(3.079788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8931005) q[0];
sx q[0];
rz(-0.19344261) q[0];
sx q[0];
rz(0.50791159) q[0];
rz(-1.6912564) q[1];
sx q[1];
rz(-2.129887) q[1];
sx q[1];
rz(1.3571665) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0944509) q[0];
sx q[0];
rz(-3.019382) q[0];
sx q[0];
rz(2.727174) q[0];
rz(-1.9903899) q[2];
sx q[2];
rz(-1.9997678) q[2];
sx q[2];
rz(1.8544514) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0101945) q[1];
sx q[1];
rz(-1.6784188) q[1];
sx q[1];
rz(-1.1788998) q[1];
x q[2];
rz(1.0413707) q[3];
sx q[3];
rz(-0.66877194) q[3];
sx q[3];
rz(-0.032162766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.42673972) q[2];
sx q[2];
rz(-1.667495) q[2];
sx q[2];
rz(-1.1023785) q[2];
rz(2.0057996) q[3];
sx q[3];
rz(-1.2165242) q[3];
sx q[3];
rz(0.77146012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.256846) q[0];
sx q[0];
rz(-0.9698292) q[0];
sx q[0];
rz(-0.92887512) q[0];
rz(-0.38201395) q[1];
sx q[1];
rz(-2.1451352) q[1];
sx q[1];
rz(1.4487723) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3253052) q[0];
sx q[0];
rz(-2.0465809) q[0];
sx q[0];
rz(1.2289117) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5097627) q[2];
sx q[2];
rz(-2.054913) q[2];
sx q[2];
rz(-2.8060437) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0121921) q[1];
sx q[1];
rz(-1.6794278) q[1];
sx q[1];
rz(1.9700178) q[1];
rz(-1.9881416) q[3];
sx q[3];
rz(-2.7355237) q[3];
sx q[3];
rz(1.690133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.53943071) q[2];
sx q[2];
rz(-0.346589) q[2];
sx q[2];
rz(2.3740785) q[2];
rz(1.2933939) q[3];
sx q[3];
rz(-0.69197217) q[3];
sx q[3];
rz(-0.20492157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2061737) q[0];
sx q[0];
rz(-0.17429166) q[0];
sx q[0];
rz(-0.25099227) q[0];
rz(-0.17768606) q[1];
sx q[1];
rz(-1.6975941) q[1];
sx q[1];
rz(-0.65690717) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11372317) q[0];
sx q[0];
rz(-1.8060469) q[0];
sx q[0];
rz(0.21224888) q[0];
rz(-pi) q[1];
rz(-2.4877704) q[2];
sx q[2];
rz(-1.6966239) q[2];
sx q[2];
rz(0.74761151) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.71437876) q[1];
sx q[1];
rz(-1.2717112) q[1];
sx q[1];
rz(-0.40451357) q[1];
rz(2.5890089) q[3];
sx q[3];
rz(-2.4173418) q[3];
sx q[3];
rz(0.16952215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3333007) q[2];
sx q[2];
rz(-0.57922816) q[2];
sx q[2];
rz(2.2283238) q[2];
rz(-2.463786) q[3];
sx q[3];
rz(-1.6401688) q[3];
sx q[3];
rz(-2.138413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-2.102757) q[0];
sx q[0];
rz(-1.1433733) q[0];
sx q[0];
rz(-2.5586149) q[0];
rz(1.673117) q[1];
sx q[1];
rz(-2.1491094) q[1];
sx q[1];
rz(2.800422) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.720168) q[0];
sx q[0];
rz(-2.9377794) q[0];
sx q[0];
rz(-1.3595306) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2874574) q[2];
sx q[2];
rz(-1.8036799) q[2];
sx q[2];
rz(0.25668677) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.22108378) q[1];
sx q[1];
rz(-2.608641) q[1];
sx q[1];
rz(2.3531662) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.736015) q[3];
sx q[3];
rz(-2.2080407) q[3];
sx q[3];
rz(3.0831856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.94656452) q[2];
sx q[2];
rz(-2.3174758) q[2];
sx q[2];
rz(-0.80671802) q[2];
rz(0.45754704) q[3];
sx q[3];
rz(-2.1541607) q[3];
sx q[3];
rz(0.88361067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6954527) q[0];
sx q[0];
rz(-1.0574295) q[0];
sx q[0];
rz(0.17644185) q[0];
rz(0.31013075) q[1];
sx q[1];
rz(-1.6519494) q[1];
sx q[1];
rz(2.2055221) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1520278) q[0];
sx q[0];
rz(-1.7195722) q[0];
sx q[0];
rz(2.8418733) q[0];
x q[1];
rz(2.0240729) q[2];
sx q[2];
rz(-2.1195115) q[2];
sx q[2];
rz(-1.1191776) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99451274) q[1];
sx q[1];
rz(-1.6638866) q[1];
sx q[1];
rz(1.2017815) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0447254) q[3];
sx q[3];
rz(-1.4634261) q[3];
sx q[3];
rz(2.3176258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.04756847) q[2];
sx q[2];
rz(-1.8393643) q[2];
sx q[2];
rz(0.0085208323) q[2];
rz(2.408037) q[3];
sx q[3];
rz(-0.80500427) q[3];
sx q[3];
rz(-0.12695299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9200639) q[0];
sx q[0];
rz(-0.41291741) q[0];
sx q[0];
rz(-0.76989663) q[0];
rz(3.0072615) q[1];
sx q[1];
rz(-0.7901935) q[1];
sx q[1];
rz(1.92164) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36496057) q[0];
sx q[0];
rz(-1.3377995) q[0];
sx q[0];
rz(-2.8316336) q[0];
rz(0.6069746) q[2];
sx q[2];
rz(-1.3262981) q[2];
sx q[2];
rz(2.2331657) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.091948) q[1];
sx q[1];
rz(-1.94749) q[1];
sx q[1];
rz(-2.6833181) q[1];
rz(-pi) q[2];
rz(2.0073529) q[3];
sx q[3];
rz(-0.16664342) q[3];
sx q[3];
rz(-1.327001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.28025815) q[2];
sx q[2];
rz(-0.65797776) q[2];
sx q[2];
rz(-0.78557837) q[2];
rz(-0.33513364) q[3];
sx q[3];
rz(-0.98888713) q[3];
sx q[3];
rz(-0.82211632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32421865) q[0];
sx q[0];
rz(-0.86295177) q[0];
sx q[0];
rz(-1.2865768) q[0];
rz(2.4556976) q[1];
sx q[1];
rz(-2.5527725) q[1];
sx q[1];
rz(1.7935161) q[1];
rz(-2.1519312) q[2];
sx q[2];
rz(-1.5753395) q[2];
sx q[2];
rz(0.01103845) q[2];
rz(1.8635345) q[3];
sx q[3];
rz(-1.6061898) q[3];
sx q[3];
rz(1.3502179) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
