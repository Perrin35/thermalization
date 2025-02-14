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
rz(0.87297451) q[0];
sx q[0];
rz(-0.31384808) q[0];
sx q[0];
rz(2.0883972) q[0];
rz(1.6246417) q[1];
sx q[1];
rz(2.8487974) q[1];
sx q[1];
rz(8.4881633) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7760843) q[0];
sx q[0];
rz(-2.0788686) q[0];
sx q[0];
rz(-3.1209141) q[0];
rz(-pi) q[1];
rz(0.99586113) q[2];
sx q[2];
rz(-2.323192) q[2];
sx q[2];
rz(-0.40413293) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9642311) q[1];
sx q[1];
rz(-1.8500016) q[1];
sx q[1];
rz(2.5482168) q[1];
rz(2.7500912) q[3];
sx q[3];
rz(-1.3927407) q[3];
sx q[3];
rz(2.4181929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34095731) q[2];
sx q[2];
rz(-0.95658797) q[2];
sx q[2];
rz(2.0802278) q[2];
rz(0.64570767) q[3];
sx q[3];
rz(-2.784745) q[3];
sx q[3];
rz(2.1319353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89479947) q[0];
sx q[0];
rz(-0.49870393) q[0];
sx q[0];
rz(-2.6820768) q[0];
rz(1.3872321) q[1];
sx q[1];
rz(-0.51001716) q[1];
sx q[1];
rz(1.0646819) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48021218) q[0];
sx q[0];
rz(-2.7837232) q[0];
sx q[0];
rz(0.46741875) q[0];
x q[1];
rz(3.0783976) q[2];
sx q[2];
rz(-0.26775751) q[2];
sx q[2];
rz(-1.7187985) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.766651) q[1];
sx q[1];
rz(-1.2514352) q[1];
sx q[1];
rz(1.4724031) q[1];
x q[2];
rz(-1.2336938) q[3];
sx q[3];
rz(-1.3114942) q[3];
sx q[3];
rz(1.269066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3620944) q[2];
sx q[2];
rz(-2.4694314) q[2];
sx q[2];
rz(-2.6008495) q[2];
rz(-0.27124673) q[3];
sx q[3];
rz(-1.8999148) q[3];
sx q[3];
rz(1.1312243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75854492) q[0];
sx q[0];
rz(-0.31065148) q[0];
sx q[0];
rz(-0.3824105) q[0];
rz(-0.27579871) q[1];
sx q[1];
rz(-0.47960061) q[1];
sx q[1];
rz(0.84024215) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29768794) q[0];
sx q[0];
rz(-0.49722222) q[0];
sx q[0];
rz(-2.3275359) q[0];
rz(0.20277299) q[2];
sx q[2];
rz(-1.2744858) q[2];
sx q[2];
rz(-1.9358097) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0419114) q[1];
sx q[1];
rz(-0.96158342) q[1];
sx q[1];
rz(0.64263338) q[1];
x q[2];
rz(0.48583416) q[3];
sx q[3];
rz(-1.6260901) q[3];
sx q[3];
rz(0.059730949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8648839) q[2];
sx q[2];
rz(-2.101892) q[2];
sx q[2];
rz(3.0806105) q[2];
rz(1.3649155) q[3];
sx q[3];
rz(-0.18326062) q[3];
sx q[3];
rz(1.1158367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0404496) q[0];
sx q[0];
rz(-1.018486) q[0];
sx q[0];
rz(-2.2818991) q[0];
rz(1.0243833) q[1];
sx q[1];
rz(-2.5442217) q[1];
sx q[1];
rz(0.76622564) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9767155) q[0];
sx q[0];
rz(-3.0857997) q[0];
sx q[0];
rz(-0.68875046) q[0];
x q[1];
rz(-1.0417074) q[2];
sx q[2];
rz(-1.4204362) q[2];
sx q[2];
rz(2.9632115) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1376201) q[1];
sx q[1];
rz(-1.8345873) q[1];
sx q[1];
rz(0.54171087) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4859824) q[3];
sx q[3];
rz(-2.437371) q[3];
sx q[3];
rz(-1.4798284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.031294558) q[2];
sx q[2];
rz(-0.86250192) q[2];
sx q[2];
rz(3.1150225) q[2];
rz(-2.8734112) q[3];
sx q[3];
rz(-2.6082706) q[3];
sx q[3];
rz(0.69494438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7705695) q[0];
sx q[0];
rz(-2.4257648) q[0];
sx q[0];
rz(-0.41719607) q[0];
rz(2.5542651) q[1];
sx q[1];
rz(-0.53343499) q[1];
sx q[1];
rz(0.74648285) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4338845) q[0];
sx q[0];
rz(-1.9198155) q[0];
sx q[0];
rz(-1.9093348) q[0];
rz(0.79277798) q[2];
sx q[2];
rz(-2.2054892) q[2];
sx q[2];
rz(1.9605876) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7938817) q[1];
sx q[1];
rz(-2.073596) q[1];
sx q[1];
rz(-2.0050843) q[1];
x q[2];
rz(-2.0516486) q[3];
sx q[3];
rz(-2.0084511) q[3];
sx q[3];
rz(-2.5646427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.31009659) q[2];
sx q[2];
rz(-2.8079171) q[2];
sx q[2];
rz(-0.92158544) q[2];
rz(-1.0725675) q[3];
sx q[3];
rz(-1.8662063) q[3];
sx q[3];
rz(-0.5630365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614586) q[0];
sx q[0];
rz(-0.38668329) q[0];
sx q[0];
rz(0.65163809) q[0];
rz(-2.1561275) q[1];
sx q[1];
rz(-2.5990504) q[1];
sx q[1];
rz(-2.998897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9938719) q[0];
sx q[0];
rz(-1.8447478) q[0];
sx q[0];
rz(-3.0763294) q[0];
rz(2.8593117) q[2];
sx q[2];
rz(-2.9476894) q[2];
sx q[2];
rz(1.0218203) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4139) q[1];
sx q[1];
rz(-1.6928612) q[1];
sx q[1];
rz(-0.47553606) q[1];
rz(3.0050395) q[3];
sx q[3];
rz(-1.1205744) q[3];
sx q[3];
rz(-0.10060752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5695213) q[2];
sx q[2];
rz(-0.60939747) q[2];
sx q[2];
rz(-0.93192464) q[2];
rz(-2.7052687) q[3];
sx q[3];
rz(-0.47567979) q[3];
sx q[3];
rz(-1.1599734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6465004) q[0];
sx q[0];
rz(-0.91630542) q[0];
sx q[0];
rz(1.6012993) q[0];
rz(2.1401999) q[1];
sx q[1];
rz(-1.5903558) q[1];
sx q[1];
rz(0.95759773) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0580826) q[0];
sx q[0];
rz(-1.8206114) q[0];
sx q[0];
rz(1.183038) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4740491) q[2];
sx q[2];
rz(-1.3203586) q[2];
sx q[2];
rz(1.7360404) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1349138) q[1];
sx q[1];
rz(-1.5638664) q[1];
sx q[1];
rz(2.0459941) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6041338) q[3];
sx q[3];
rz(-1.0126202) q[3];
sx q[3];
rz(-0.5287002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5074978) q[2];
sx q[2];
rz(-1.2853421) q[2];
sx q[2];
rz(-2.0095339) q[2];
rz(0.13103983) q[3];
sx q[3];
rz(-1.1538426) q[3];
sx q[3];
rz(-1.2426144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89192724) q[0];
sx q[0];
rz(-2.6234493) q[0];
sx q[0];
rz(-0.074935496) q[0];
rz(2.9451008) q[1];
sx q[1];
rz(-0.34759977) q[1];
sx q[1];
rz(0.67952716) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20048902) q[0];
sx q[0];
rz(-1.3692432) q[0];
sx q[0];
rz(2.5747712) q[0];
rz(-pi) q[1];
rz(0.029506186) q[2];
sx q[2];
rz(-1.1669807) q[2];
sx q[2];
rz(-0.6617512) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36489242) q[1];
sx q[1];
rz(-0.51799315) q[1];
sx q[1];
rz(0.97229506) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6305168) q[3];
sx q[3];
rz(-1.5580439) q[3];
sx q[3];
rz(-1.2619357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.8276662) q[2];
sx q[2];
rz(-0.9845261) q[2];
sx q[2];
rz(1.0796245) q[2];
rz(-2.6880654) q[3];
sx q[3];
rz(-2.270416) q[3];
sx q[3];
rz(-0.96550226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6249348) q[0];
sx q[0];
rz(-2.9624532) q[0];
sx q[0];
rz(-3.0058885) q[0];
rz(-0.16595674) q[1];
sx q[1];
rz(-0.9981007) q[1];
sx q[1];
rz(2.1735294) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5639599) q[0];
sx q[0];
rz(-1.5819307) q[0];
sx q[0];
rz(-1.0554814) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14974515) q[2];
sx q[2];
rz(-1.0335084) q[2];
sx q[2];
rz(-1.6884691) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.49478) q[1];
sx q[1];
rz(-0.54807263) q[1];
sx q[1];
rz(-0.41278028) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3126722) q[3];
sx q[3];
rz(-1.934623) q[3];
sx q[3];
rz(2.8175333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3193937) q[2];
sx q[2];
rz(-2.2182756) q[2];
sx q[2];
rz(0.53259069) q[2];
rz(0.79832625) q[3];
sx q[3];
rz(-2.7440378) q[3];
sx q[3];
rz(-3.047191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3548729) q[0];
sx q[0];
rz(-0.69364554) q[0];
sx q[0];
rz(2.4548446) q[0];
rz(1.2314388) q[1];
sx q[1];
rz(-1.8309007) q[1];
sx q[1];
rz(3.0795857) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8221164) q[0];
sx q[0];
rz(-1.5052759) q[0];
sx q[0];
rz(0.70573745) q[0];
rz(-2.9509281) q[2];
sx q[2];
rz(-1.4989509) q[2];
sx q[2];
rz(2.004625) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1126642) q[1];
sx q[1];
rz(-2.8920569) q[1];
sx q[1];
rz(2.1965532) q[1];
rz(-0.047895821) q[3];
sx q[3];
rz(-2.8433099) q[3];
sx q[3];
rz(-3.0259231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9318781) q[2];
sx q[2];
rz(-0.96499062) q[2];
sx q[2];
rz(-2.9689201) q[2];
rz(0.73575819) q[3];
sx q[3];
rz(-0.57307214) q[3];
sx q[3];
rz(-2.6652523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
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
rz(-2.2810777) q[1];
sx q[1];
rz(-1.9445226) q[1];
sx q[1];
rz(-1.3938211) q[1];
rz(-1.4635659) q[2];
sx q[2];
rz(-1.4489531) q[2];
sx q[2];
rz(1.0482241) q[2];
rz(0.53027897) q[3];
sx q[3];
rz(-2.1580631) q[3];
sx q[3];
rz(-2.8736339) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
