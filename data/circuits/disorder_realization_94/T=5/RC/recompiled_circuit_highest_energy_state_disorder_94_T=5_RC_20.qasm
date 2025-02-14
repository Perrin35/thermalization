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
rz(-2.2686181) q[0];
sx q[0];
rz(-2.8277446) q[0];
sx q[0];
rz(-2.0883972) q[0];
rz(1.6246417) q[1];
sx q[1];
rz(2.8487974) q[1];
sx q[1];
rz(8.4881633) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7760843) q[0];
sx q[0];
rz(-1.062724) q[0];
sx q[0];
rz(0.020678542) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83990015) q[2];
sx q[2];
rz(-1.9790302) q[2];
sx q[2];
rz(-0.74980914) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0055252) q[1];
sx q[1];
rz(-2.4930291) q[1];
sx q[1];
rz(2.6678209) q[1];
x q[2];
rz(-0.4406919) q[3];
sx q[3];
rz(-2.7134136) q[3];
sx q[3];
rz(2.6994914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8006353) q[2];
sx q[2];
rz(-0.95658797) q[2];
sx q[2];
rz(-2.0802278) q[2];
rz(-0.64570767) q[3];
sx q[3];
rz(-2.784745) q[3];
sx q[3];
rz(-2.1319353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89479947) q[0];
sx q[0];
rz(-0.49870393) q[0];
sx q[0];
rz(2.6820768) q[0];
rz(1.7543606) q[1];
sx q[1];
rz(-0.51001716) q[1];
sx q[1];
rz(-1.0646819) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6094006) q[0];
sx q[0];
rz(-1.7292892) q[0];
sx q[0];
rz(-2.81937) q[0];
rz(-3.0783976) q[2];
sx q[2];
rz(-2.8738351) q[2];
sx q[2];
rz(1.4227941) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.83512703) q[1];
sx q[1];
rz(-1.4773932) q[1];
sx q[1];
rz(-2.8207835) q[1];
x q[2];
rz(1.9078988) q[3];
sx q[3];
rz(-1.8300984) q[3];
sx q[3];
rz(-1.269066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77949828) q[2];
sx q[2];
rz(-2.4694314) q[2];
sx q[2];
rz(-0.54074311) q[2];
rz(0.27124673) q[3];
sx q[3];
rz(-1.2416779) q[3];
sx q[3];
rz(1.1312243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3830477) q[0];
sx q[0];
rz(-0.31065148) q[0];
sx q[0];
rz(2.7591822) q[0];
rz(-0.27579871) q[1];
sx q[1];
rz(-2.661992) q[1];
sx q[1];
rz(-0.84024215) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58038515) q[0];
sx q[0];
rz(-1.9044283) q[0];
sx q[0];
rz(-1.946627) q[0];
rz(-pi) q[1];
rz(-2.1539168) q[2];
sx q[2];
rz(-2.7842369) q[2];
sx q[2];
rz(-0.59229702) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.099681252) q[1];
sx q[1];
rz(-2.1800092) q[1];
sx q[1];
rz(-0.64263338) q[1];
rz(-pi) q[2];
rz(0.11798604) q[3];
sx q[3];
rz(-2.6528721) q[3];
sx q[3];
rz(1.5262469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8648839) q[2];
sx q[2];
rz(-1.0397006) q[2];
sx q[2];
rz(-3.0806105) q[2];
rz(1.3649155) q[3];
sx q[3];
rz(-0.18326062) q[3];
sx q[3];
rz(-2.0257559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0404496) q[0];
sx q[0];
rz(-1.018486) q[0];
sx q[0];
rz(-2.2818991) q[0];
rz(2.1172093) q[1];
sx q[1];
rz(-0.59737098) q[1];
sx q[1];
rz(0.76622564) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47536248) q[0];
sx q[0];
rz(-1.6138617) q[0];
sx q[0];
rz(-1.5353139) q[0];
rz(-pi) q[1];
rz(-0.17373093) q[2];
sx q[2];
rz(-2.0933009) q[2];
sx q[2];
rz(-1.3050543) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7304036) q[1];
sx q[1];
rz(-1.0497739) q[1];
sx q[1];
rz(-1.26544) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2732417) q[3];
sx q[3];
rz(-1.5159226) q[3];
sx q[3];
rz(-3.1153278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.031294558) q[2];
sx q[2];
rz(-2.2790907) q[2];
sx q[2];
rz(-0.026570126) q[2];
rz(-2.8734112) q[3];
sx q[3];
rz(-0.53332204) q[3];
sx q[3];
rz(2.4466483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37102315) q[0];
sx q[0];
rz(-2.4257648) q[0];
sx q[0];
rz(0.41719607) q[0];
rz(0.5873276) q[1];
sx q[1];
rz(-0.53343499) q[1];
sx q[1];
rz(2.3951098) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5076818) q[0];
sx q[0];
rz(-0.48134781) q[0];
sx q[0];
rz(-0.7397299) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3801489) q[2];
sx q[2];
rz(-2.1816744) q[2];
sx q[2];
rz(2.2100984) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6986148) q[1];
sx q[1];
rz(-1.1931975) q[1];
sx q[1];
rz(-2.5966132) q[1];
rz(-1.0899441) q[3];
sx q[3];
rz(-2.0084511) q[3];
sx q[3];
rz(-0.57695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.31009659) q[2];
sx q[2];
rz(-2.8079171) q[2];
sx q[2];
rz(0.92158544) q[2];
rz(1.0725675) q[3];
sx q[3];
rz(-1.8662063) q[3];
sx q[3];
rz(-2.5785562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6801341) q[0];
sx q[0];
rz(-0.38668329) q[0];
sx q[0];
rz(2.4899546) q[0];
rz(-0.98546511) q[1];
sx q[1];
rz(-2.5990504) q[1];
sx q[1];
rz(-0.14269565) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9938719) q[0];
sx q[0];
rz(-1.8447478) q[0];
sx q[0];
rz(0.065263211) q[0];
rz(-pi) q[1];
rz(2.9551836) q[2];
sx q[2];
rz(-1.6244955) q[2];
sx q[2];
rz(2.3153564) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0662996) q[1];
sx q[1];
rz(-2.651804) q[1];
sx q[1];
rz(-2.879786) q[1];
rz(1.296259) q[3];
sx q[3];
rz(-0.46911383) q[3];
sx q[3];
rz(-0.20524552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5720713) q[2];
sx q[2];
rz(-0.60939747) q[2];
sx q[2];
rz(-2.209668) q[2];
rz(2.7052687) q[3];
sx q[3];
rz(-2.6659129) q[3];
sx q[3];
rz(1.9816192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49509224) q[0];
sx q[0];
rz(-2.2252872) q[0];
sx q[0];
rz(-1.5402933) q[0];
rz(-2.1401999) q[1];
sx q[1];
rz(-1.5903558) q[1];
sx q[1];
rz(2.1839949) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0315899) q[0];
sx q[0];
rz(-2.6837807) q[0];
sx q[0];
rz(2.1643967) q[0];
rz(-2.7805336) q[2];
sx q[2];
rz(-2.8734837) q[2];
sx q[2];
rz(-2.1092871) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.56768306) q[1];
sx q[1];
rz(-2.0459818) q[1];
sx q[1];
rz(-0.0077933776) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0882629) q[3];
sx q[3];
rz(-2.5825273) q[3];
sx q[3];
rz(-2.6757765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6340948) q[2];
sx q[2];
rz(-1.8562506) q[2];
sx q[2];
rz(2.0095339) q[2];
rz(-0.13103983) q[3];
sx q[3];
rz(-1.9877501) q[3];
sx q[3];
rz(1.8989782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2496654) q[0];
sx q[0];
rz(-0.51814336) q[0];
sx q[0];
rz(0.074935496) q[0];
rz(-2.9451008) q[1];
sx q[1];
rz(-2.7939929) q[1];
sx q[1];
rz(0.67952716) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4970447) q[0];
sx q[0];
rz(-2.1247852) q[0];
sx q[0];
rz(-1.333167) q[0];
rz(-pi) q[1];
x q[1];
rz(0.029506186) q[2];
sx q[2];
rz(-1.1669807) q[2];
sx q[2];
rz(2.4798415) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36489242) q[1];
sx q[1];
rz(-0.51799315) q[1];
sx q[1];
rz(-0.97229506) q[1];
rz(-0.51107589) q[3];
sx q[3];
rz(-1.5580439) q[3];
sx q[3];
rz(1.8796569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3139265) q[2];
sx q[2];
rz(-2.1570666) q[2];
sx q[2];
rz(-2.0619681) q[2];
rz(-2.6880654) q[3];
sx q[3];
rz(-2.270416) q[3];
sx q[3];
rz(-0.96550226) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6249348) q[0];
sx q[0];
rz(-2.9624532) q[0];
sx q[0];
rz(3.0058885) q[0];
rz(0.16595674) q[1];
sx q[1];
rz(-0.9981007) q[1];
sx q[1];
rz(0.96806324) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0128202) q[0];
sx q[0];
rz(-2.6261683) q[0];
sx q[0];
rz(1.5933871) q[0];
rz(-pi) q[1];
rz(-2.9918475) q[2];
sx q[2];
rz(-2.1080842) q[2];
sx q[2];
rz(-1.4531236) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.49478) q[1];
sx q[1];
rz(-0.54807263) q[1];
sx q[1];
rz(0.41278028) q[1];
rz(-pi) q[2];
rz(-1.3126722) q[3];
sx q[3];
rz(-1.2069697) q[3];
sx q[3];
rz(-2.8175333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82219899) q[2];
sx q[2];
rz(-2.2182756) q[2];
sx q[2];
rz(-0.53259069) q[2];
rz(-0.79832625) q[3];
sx q[3];
rz(-0.39755487) q[3];
sx q[3];
rz(-3.047191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7867197) q[0];
sx q[0];
rz(-0.69364554) q[0];
sx q[0];
rz(0.68674809) q[0];
rz(-1.2314388) q[1];
sx q[1];
rz(-1.8309007) q[1];
sx q[1];
rz(-3.0795857) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4486478) q[0];
sx q[0];
rz(-0.8668859) q[0];
sx q[0];
rz(-1.6567898) q[0];
rz(-pi) q[1];
rz(-0.19066452) q[2];
sx q[2];
rz(-1.4989509) q[2];
sx q[2];
rz(-2.004625) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0289284) q[1];
sx q[1];
rz(-0.24953574) q[1];
sx q[1];
rz(-2.1965532) q[1];
rz(0.047895821) q[3];
sx q[3];
rz(-2.8433099) q[3];
sx q[3];
rz(3.0259231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.20971458) q[2];
sx q[2];
rz(-0.96499062) q[2];
sx q[2];
rz(-2.9689201) q[2];
rz(-0.73575819) q[3];
sx q[3];
rz(-2.5685205) q[3];
sx q[3];
rz(0.47634038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.944735) q[0];
sx q[0];
rz(-1.6898962) q[0];
sx q[0];
rz(1.663399) q[0];
rz(-0.860515) q[1];
sx q[1];
rz(-1.1970701) q[1];
sx q[1];
rz(1.7477716) q[1];
rz(-1.6780268) q[2];
sx q[2];
rz(-1.6926395) q[2];
sx q[2];
rz(-2.0933685) q[2];
rz(2.2279578) q[3];
sx q[3];
rz(-1.1362094) q[3];
sx q[3];
rz(2.1528578) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
