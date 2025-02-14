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
rz(-1.516951) q[1];
sx q[1];
rz(-2.8487974) q[1];
sx q[1];
rz(2.204978) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36550831) q[0];
sx q[0];
rz(-2.0788686) q[0];
sx q[0];
rz(3.1209141) q[0];
rz(0.52626558) q[2];
sx q[2];
rz(-2.230245) q[2];
sx q[2];
rz(1.9786722) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0055252) q[1];
sx q[1];
rz(-2.4930291) q[1];
sx q[1];
rz(-0.4737718) q[1];
rz(-2.7500912) q[3];
sx q[3];
rz(-1.3927407) q[3];
sx q[3];
rz(-2.4181929) q[3];
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
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2467932) q[0];
sx q[0];
rz(-0.49870393) q[0];
sx q[0];
rz(-0.45951581) q[0];
rz(1.7543606) q[1];
sx q[1];
rz(-0.51001716) q[1];
sx q[1];
rz(-1.0646819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6613805) q[0];
sx q[0];
rz(-0.35786942) q[0];
sx q[0];
rz(-0.46741875) q[0];
rz(-pi) q[1];
rz(1.5534723) q[2];
sx q[2];
rz(-1.8380062) q[2];
sx q[2];
rz(1.4883177) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.83512703) q[1];
sx q[1];
rz(-1.4773932) q[1];
sx q[1];
rz(0.3208092) q[1];
rz(-pi) q[2];
rz(-1.9078988) q[3];
sx q[3];
rz(-1.8300984) q[3];
sx q[3];
rz(-1.8725267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.77949828) q[2];
sx q[2];
rz(-2.4694314) q[2];
sx q[2];
rz(2.6008495) q[2];
rz(-0.27124673) q[3];
sx q[3];
rz(-1.8999148) q[3];
sx q[3];
rz(-2.0103683) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75854492) q[0];
sx q[0];
rz(-2.8309412) q[0];
sx q[0];
rz(0.3824105) q[0];
rz(-0.27579871) q[1];
sx q[1];
rz(-2.661992) q[1];
sx q[1];
rz(-0.84024215) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8439047) q[0];
sx q[0];
rz(-2.6443704) q[0];
sx q[0];
rz(2.3275359) q[0];
rz(0.98767583) q[2];
sx q[2];
rz(-0.35735574) q[2];
sx q[2];
rz(-2.5492956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0178595) q[1];
sx q[1];
rz(-2.2869733) q[1];
sx q[1];
rz(-2.2804428) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48583416) q[3];
sx q[3];
rz(-1.5155025) q[3];
sx q[3];
rz(-0.059730949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27670878) q[2];
sx q[2];
rz(-2.101892) q[2];
sx q[2];
rz(-3.0806105) q[2];
rz(1.7766772) q[3];
sx q[3];
rz(-0.18326062) q[3];
sx q[3];
rz(2.0257559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.10114305) q[0];
sx q[0];
rz(-2.1231066) q[0];
sx q[0];
rz(0.85969353) q[0];
rz(-1.0243833) q[1];
sx q[1];
rz(-0.59737098) q[1];
sx q[1];
rz(0.76622564) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0939056) q[0];
sx q[0];
rz(-1.5353468) q[0];
sx q[0];
rz(-3.0985002) q[0];
x q[1];
rz(-1.2791951) q[2];
sx q[2];
rz(-2.5935141) q[2];
sx q[2];
rz(1.4983791) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0039725) q[1];
sx q[1];
rz(-1.3070053) q[1];
sx q[1];
rz(0.54171087) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4859824) q[3];
sx q[3];
rz(-0.70422164) q[3];
sx q[3];
rz(1.6617642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.031294558) q[2];
sx q[2];
rz(-0.86250192) q[2];
sx q[2];
rz(0.026570126) q[2];
rz(-2.8734112) q[3];
sx q[3];
rz(-0.53332204) q[3];
sx q[3];
rz(2.4466483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37102315) q[0];
sx q[0];
rz(-2.4257648) q[0];
sx q[0];
rz(2.7243966) q[0];
rz(-0.5873276) q[1];
sx q[1];
rz(-0.53343499) q[1];
sx q[1];
rz(-2.3951098) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70770812) q[0];
sx q[0];
rz(-1.9198155) q[0];
sx q[0];
rz(1.2322578) q[0];
rz(0.76144371) q[2];
sx q[2];
rz(-0.9599182) q[2];
sx q[2];
rz(-0.93149422) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5601859) q[1];
sx q[1];
rz(-2.489632) q[1];
sx q[1];
rz(-0.65309872) q[1];
x q[2];
rz(1.0899441) q[3];
sx q[3];
rz(-2.0084511) q[3];
sx q[3];
rz(-2.5646427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8314961) q[2];
sx q[2];
rz(-0.3336755) q[2];
sx q[2];
rz(-0.92158544) q[2];
rz(2.0690252) q[3];
sx q[3];
rz(-1.8662063) q[3];
sx q[3];
rz(2.5785562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6801341) q[0];
sx q[0];
rz(-0.38668329) q[0];
sx q[0];
rz(-2.4899546) q[0];
rz(-0.98546511) q[1];
sx q[1];
rz(-0.54254222) q[1];
sx q[1];
rz(0.14269565) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9938719) q[0];
sx q[0];
rz(-1.2968449) q[0];
sx q[0];
rz(-0.065263211) q[0];
x q[1];
rz(2.8593117) q[2];
sx q[2];
rz(-0.19390324) q[2];
sx q[2];
rz(2.1197723) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.922077) q[1];
sx q[1];
rz(-2.0425046) q[1];
sx q[1];
rz(-1.4336777) q[1];
rz(-pi) q[2];
rz(-1.296259) q[3];
sx q[3];
rz(-2.6724788) q[3];
sx q[3];
rz(-0.20524552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5695213) q[2];
sx q[2];
rz(-2.5321952) q[2];
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
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49509224) q[0];
sx q[0];
rz(-2.2252872) q[0];
sx q[0];
rz(-1.6012993) q[0];
rz(2.1401999) q[1];
sx q[1];
rz(-1.5903558) q[1];
sx q[1];
rz(0.95759773) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0580826) q[0];
sx q[0];
rz(-1.3209813) q[0];
sx q[0];
rz(1.183038) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7805336) q[2];
sx q[2];
rz(-2.8734837) q[2];
sx q[2];
rz(1.0323056) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5909429) q[1];
sx q[1];
rz(-0.47524449) q[1];
sx q[1];
rz(-1.5859423) q[1];
rz(3.0882629) q[3];
sx q[3];
rz(-0.55906534) q[3];
sx q[3];
rz(-2.6757765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6340948) q[2];
sx q[2];
rz(-1.8562506) q[2];
sx q[2];
rz(1.1320587) q[2];
rz(-0.13103983) q[3];
sx q[3];
rz(-1.9877501) q[3];
sx q[3];
rz(-1.2426144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89192724) q[0];
sx q[0];
rz(-0.51814336) q[0];
sx q[0];
rz(3.0666572) q[0];
rz(-0.19649188) q[1];
sx q[1];
rz(-0.34759977) q[1];
sx q[1];
rz(0.67952716) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0760114) q[0];
sx q[0];
rz(-0.59787321) q[0];
sx q[0];
rz(-0.36361106) q[0];
rz(-pi) q[1];
rz(-1.9747693) q[2];
sx q[2];
rz(-1.543664) q[2];
sx q[2];
rz(-2.2441442) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4008182) q[1];
sx q[1];
rz(-1.8535103) q[1];
sx q[1];
rz(-1.1307471) q[1];
rz(0.026067928) q[3];
sx q[3];
rz(-0.51122087) q[3];
sx q[3];
rz(2.8554684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3139265) q[2];
sx q[2];
rz(-0.9845261) q[2];
sx q[2];
rz(1.0796245) q[2];
rz(2.6880654) q[3];
sx q[3];
rz(-0.87117666) q[3];
sx q[3];
rz(-0.96550226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51665783) q[0];
sx q[0];
rz(-2.9624532) q[0];
sx q[0];
rz(3.0058885) q[0];
rz(-0.16595674) q[1];
sx q[1];
rz(-2.143492) q[1];
sx q[1];
rz(0.96806324) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5639599) q[0];
sx q[0];
rz(-1.559662) q[0];
sx q[0];
rz(1.0554814) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8161681) q[2];
sx q[2];
rz(-2.5858063) q[2];
sx q[2];
rz(-1.1664558) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6468127) q[1];
sx q[1];
rz(-2.59352) q[1];
sx q[1];
rz(0.41278028) q[1];
rz(-pi) q[2];
rz(2.7664207) q[3];
sx q[3];
rz(-1.811677) q[3];
sx q[3];
rz(1.3404121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3193937) q[2];
sx q[2];
rz(-2.2182756) q[2];
sx q[2];
rz(-2.609002) q[2];
rz(0.79832625) q[3];
sx q[3];
rz(-0.39755487) q[3];
sx q[3];
rz(3.047191) q[3];
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
rz(pi/2) q[2];
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
rz(-1.310692) q[1];
sx q[1];
rz(-3.0795857) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3162295) q[0];
sx q[0];
rz(-0.70825173) q[0];
sx q[0];
rz(3.0407719) q[0];
rz(-pi) q[1];
rz(-1.4976296) q[2];
sx q[2];
rz(-1.760963) q[2];
sx q[2];
rz(-2.721618) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15276423) q[1];
sx q[1];
rz(-1.4256434) q[1];
sx q[1];
rz(1.7744906) q[1];
x q[2];
rz(-0.047895821) q[3];
sx q[3];
rz(-2.8433099) q[3];
sx q[3];
rz(0.1156696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9318781) q[2];
sx q[2];
rz(-0.96499062) q[2];
sx q[2];
rz(2.9689201) q[2];
rz(-2.4058345) q[3];
sx q[3];
rz(-2.5685205) q[3];
sx q[3];
rz(2.6652523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
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
