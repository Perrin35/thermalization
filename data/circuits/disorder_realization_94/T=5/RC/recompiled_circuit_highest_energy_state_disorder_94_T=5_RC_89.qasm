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
rz(-0.2927953) q[1];
sx q[1];
rz(-2.204978) q[1];
x q[2];
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
rz(0.52626558) q[2];
sx q[2];
rz(-0.91134763) q[2];
sx q[2];
rz(1.1629205) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.17736152) q[1];
sx q[1];
rz(-1.8500016) q[1];
sx q[1];
rz(2.5482168) q[1];
x q[2];
rz(-1.7630834) q[3];
sx q[3];
rz(-1.1858127) q[3];
sx q[3];
rz(-0.92038233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.34095731) q[2];
sx q[2];
rz(-2.1850047) q[2];
sx q[2];
rz(-1.0613649) q[2];
rz(-0.64570767) q[3];
sx q[3];
rz(-0.3568477) q[3];
sx q[3];
rz(2.1319353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89479947) q[0];
sx q[0];
rz(-2.6428887) q[0];
sx q[0];
rz(-0.45951581) q[0];
rz(1.3872321) q[1];
sx q[1];
rz(-2.6315755) q[1];
sx q[1];
rz(-1.0646819) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1275528) q[0];
sx q[0];
rz(-1.2527553) q[0];
sx q[0];
rz(-1.4038588) q[0];
rz(-pi) q[1];
rz(-3.0783976) q[2];
sx q[2];
rz(-2.8738351) q[2];
sx q[2];
rz(1.4227941) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.766651) q[1];
sx q[1];
rz(-1.2514352) q[1];
sx q[1];
rz(1.4724031) q[1];
rz(-2.2467733) q[3];
sx q[3];
rz(-2.7193391) q[3];
sx q[3];
rz(-2.811712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77949828) q[2];
sx q[2];
rz(-2.4694314) q[2];
sx q[2];
rz(-0.54074311) q[2];
rz(-0.27124673) q[3];
sx q[3];
rz(-1.2416779) q[3];
sx q[3];
rz(-1.1312243) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3830477) q[0];
sx q[0];
rz(-0.31065148) q[0];
sx q[0];
rz(-2.7591822) q[0];
rz(2.8657939) q[1];
sx q[1];
rz(-0.47960061) q[1];
sx q[1];
rz(0.84024215) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5612075) q[0];
sx q[0];
rz(-1.9044283) q[0];
sx q[0];
rz(-1.946627) q[0];
rz(-pi) q[1];
rz(-2.9388197) q[2];
sx q[2];
rz(-1.2744858) q[2];
sx q[2];
rz(1.2057829) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0663755) q[1];
sx q[1];
rz(-2.0845958) q[1];
sx q[1];
rz(-0.85388501) q[1];
rz(-pi) q[2];
rz(-0.48583416) q[3];
sx q[3];
rz(-1.5155025) q[3];
sx q[3];
rz(-3.0818617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27670878) q[2];
sx q[2];
rz(-1.0397006) q[2];
sx q[2];
rz(3.0806105) q[2];
rz(1.3649155) q[3];
sx q[3];
rz(-2.958332) q[3];
sx q[3];
rz(-1.1158367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0404496) q[0];
sx q[0];
rz(-2.1231066) q[0];
sx q[0];
rz(2.2818991) q[0];
rz(2.1172093) q[1];
sx q[1];
rz(-0.59737098) q[1];
sx q[1];
rz(0.76622564) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9767155) q[0];
sx q[0];
rz(-0.055792965) q[0];
sx q[0];
rz(-0.68875046) q[0];
x q[1];
rz(1.8623975) q[2];
sx q[2];
rz(-2.5935141) q[2];
sx q[2];
rz(-1.6432135) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0039725) q[1];
sx q[1];
rz(-1.3070053) q[1];
sx q[1];
rz(0.54171087) q[1];
x q[2];
rz(-3.0697508) q[3];
sx q[3];
rz(-0.86962442) q[3];
sx q[3];
rz(1.5909242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.031294558) q[2];
sx q[2];
rz(-2.2790907) q[2];
sx q[2];
rz(0.026570126) q[2];
rz(-0.26818141) q[3];
sx q[3];
rz(-2.6082706) q[3];
sx q[3];
rz(-0.69494438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37102315) q[0];
sx q[0];
rz(-0.71582782) q[0];
sx q[0];
rz(-2.7243966) q[0];
rz(-0.5873276) q[1];
sx q[1];
rz(-2.6081577) q[1];
sx q[1];
rz(-0.74648285) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74325753) q[0];
sx q[0];
rz(-1.2534089) q[0];
sx q[0];
rz(2.7733735) q[0];
x q[1];
rz(-2.339614) q[2];
sx q[2];
rz(-2.1714513) q[2];
sx q[2];
rz(3.0026312) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.44297781) q[1];
sx q[1];
rz(-1.9483951) q[1];
sx q[1];
rz(0.54497949) q[1];
rz(-pi) q[2];
rz(2.3619805) q[3];
sx q[3];
rz(-0.63842183) q[3];
sx q[3];
rz(-0.31159952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.31009659) q[2];
sx q[2];
rz(-2.8079171) q[2];
sx q[2];
rz(-2.2200072) q[2];
rz(2.0690252) q[3];
sx q[3];
rz(-1.2753863) q[3];
sx q[3];
rz(0.5630365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6801341) q[0];
sx q[0];
rz(-2.7549094) q[0];
sx q[0];
rz(-0.65163809) q[0];
rz(0.98546511) q[1];
sx q[1];
rz(-0.54254222) q[1];
sx q[1];
rz(2.998897) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0522766) q[0];
sx q[0];
rz(-0.28142774) q[0];
sx q[0];
rz(-1.3427585) q[0];
rz(-pi) q[1];
rz(1.6254403) q[2];
sx q[2];
rz(-1.7569336) q[2];
sx q[2];
rz(-2.4071549) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.922077) q[1];
sx q[1];
rz(-1.0990881) q[1];
sx q[1];
rz(-1.4336777) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1168994) q[3];
sx q[3];
rz(-1.6936692) q[3];
sx q[3];
rz(1.6116796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5695213) q[2];
sx q[2];
rz(-0.60939747) q[2];
sx q[2];
rz(2.209668) q[2];
rz(-0.43632397) q[3];
sx q[3];
rz(-2.6659129) q[3];
sx q[3];
rz(-1.1599734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(0.49509224) q[0];
sx q[0];
rz(-0.91630542) q[0];
sx q[0];
rz(-1.5402933) q[0];
rz(-1.0013927) q[1];
sx q[1];
rz(-1.5512369) q[1];
sx q[1];
rz(2.1839949) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0580826) q[0];
sx q[0];
rz(-1.8206114) q[0];
sx q[0];
rz(-1.183038) q[0];
rz(-pi) q[1];
rz(-1.4740491) q[2];
sx q[2];
rz(-1.3203586) q[2];
sx q[2];
rz(1.4055523) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5739096) q[1];
sx q[1];
rz(-2.0459818) q[1];
sx q[1];
rz(3.1337993) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.05332975) q[3];
sx q[3];
rz(-0.55906534) q[3];
sx q[3];
rz(0.4658162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6340948) q[2];
sx q[2];
rz(-1.8562506) q[2];
sx q[2];
rz(-1.1320587) q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(3.0666572) q[0];
rz(2.9451008) q[1];
sx q[1];
rz(-2.7939929) q[1];
sx q[1];
rz(2.4620655) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.644548) q[0];
sx q[0];
rz(-2.1247852) q[0];
sx q[0];
rz(1.8084256) q[0];
x q[1];
rz(-1.9747693) q[2];
sx q[2];
rz(-1.543664) q[2];
sx q[2];
rz(-2.2441442) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7767002) q[1];
sx q[1];
rz(-0.51799315) q[1];
sx q[1];
rz(0.97229506) q[1];
x q[2];
rz(-0.026067928) q[3];
sx q[3];
rz(-2.6303718) q[3];
sx q[3];
rz(2.8554684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.8276662) q[2];
sx q[2];
rz(-2.1570666) q[2];
sx q[2];
rz(1.0796245) q[2];
rz(-0.45352724) q[3];
sx q[3];
rz(-2.270416) q[3];
sx q[3];
rz(-2.1760904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51665783) q[0];
sx q[0];
rz(-0.17913945) q[0];
sx q[0];
rz(3.0058885) q[0];
rz(0.16595674) q[1];
sx q[1];
rz(-0.9981007) q[1];
sx q[1];
rz(0.96806324) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5639599) q[0];
sx q[0];
rz(-1.5819307) q[0];
sx q[0];
rz(2.0861113) q[0];
rz(1.8161681) q[2];
sx q[2];
rz(-0.55578631) q[2];
sx q[2];
rz(-1.1664558) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6468127) q[1];
sx q[1];
rz(-2.59352) q[1];
sx q[1];
rz(-0.41278028) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37517199) q[3];
sx q[3];
rz(-1.3299156) q[3];
sx q[3];
rz(1.3404121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82219899) q[2];
sx q[2];
rz(-2.2182756) q[2];
sx q[2];
rz(-0.53259069) q[2];
rz(0.79832625) q[3];
sx q[3];
rz(-2.7440378) q[3];
sx q[3];
rz(-3.047191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3548729) q[0];
sx q[0];
rz(-2.4479471) q[0];
sx q[0];
rz(-2.4548446) q[0];
rz(-1.2314388) q[1];
sx q[1];
rz(-1.310692) q[1];
sx q[1];
rz(-0.062006921) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82536316) q[0];
sx q[0];
rz(-2.4333409) q[0];
sx q[0];
rz(0.10082074) q[0];
x q[1];
rz(-0.36293928) q[2];
sx q[2];
rz(-0.20359765) q[2];
sx q[2];
rz(-0.789895) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.15276423) q[1];
sx q[1];
rz(-1.4256434) q[1];
sx q[1];
rz(1.3671021) q[1];
rz(-pi) q[2];
rz(2.8436321) q[3];
sx q[3];
rz(-1.5567257) q[3];
sx q[3];
rz(1.5009106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.20971458) q[2];
sx q[2];
rz(-0.96499062) q[2];
sx q[2];
rz(-2.9689201) q[2];
rz(-0.73575819) q[3];
sx q[3];
rz(-2.5685205) q[3];
sx q[3];
rz(-2.6652523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.944735) q[0];
sx q[0];
rz(-1.4516964) q[0];
sx q[0];
rz(-1.4781937) q[0];
rz(2.2810777) q[1];
sx q[1];
rz(-1.1970701) q[1];
sx q[1];
rz(1.7477716) q[1];
rz(0.12254006) q[2];
sx q[2];
rz(-1.6772288) q[2];
sx q[2];
rz(2.605938) q[2];
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
