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
rz(-0.2927953) q[1];
sx q[1];
rz(-2.204978) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40799403) q[0];
sx q[0];
rz(-2.6331365) q[0];
sx q[0];
rz(-1.6079128) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6153271) q[2];
sx q[2];
rz(-2.230245) q[2];
sx q[2];
rz(1.1629205) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.577212) q[1];
sx q[1];
rz(-1.0033157) q[1];
sx q[1];
rz(1.903731) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39150146) q[3];
sx q[3];
rz(-1.3927407) q[3];
sx q[3];
rz(-2.4181929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8006353) q[2];
sx q[2];
rz(-0.95658797) q[2];
sx q[2];
rz(-2.0802278) q[2];
rz(2.495885) q[3];
sx q[3];
rz(-0.3568477) q[3];
sx q[3];
rz(2.1319353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89479947) q[0];
sx q[0];
rz(-2.6428887) q[0];
sx q[0];
rz(2.6820768) q[0];
rz(-1.7543606) q[1];
sx q[1];
rz(-0.51001716) q[1];
sx q[1];
rz(1.0646819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014039847) q[0];
sx q[0];
rz(-1.8888374) q[0];
sx q[0];
rz(-1.4038588) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0783976) q[2];
sx q[2];
rz(-0.26775751) q[2];
sx q[2];
rz(1.4227941) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.766651) q[1];
sx q[1];
rz(-1.8901575) q[1];
sx q[1];
rz(-1.4724031) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8675686) q[3];
sx q[3];
rz(-1.245387) q[3];
sx q[3];
rz(-2.7502378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3620944) q[2];
sx q[2];
rz(-2.4694314) q[2];
sx q[2];
rz(-0.54074311) q[2];
rz(0.27124673) q[3];
sx q[3];
rz(-1.8999148) q[3];
sx q[3];
rz(-1.1312243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3830477) q[0];
sx q[0];
rz(-0.31065148) q[0];
sx q[0];
rz(-0.3824105) q[0];
rz(-2.8657939) q[1];
sx q[1];
rz(-0.47960061) q[1];
sx q[1];
rz(2.3013505) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0226757) q[0];
sx q[0];
rz(-1.9249601) q[0];
sx q[0];
rz(0.35665956) q[0];
rz(2.1539168) q[2];
sx q[2];
rz(-0.35735574) q[2];
sx q[2];
rz(2.5492956) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0663755) q[1];
sx q[1];
rz(-1.0569968) q[1];
sx q[1];
rz(-2.2877076) q[1];
rz(2.6557585) q[3];
sx q[3];
rz(-1.5155025) q[3];
sx q[3];
rz(-3.0818617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.27670878) q[2];
sx q[2];
rz(-1.0397006) q[2];
sx q[2];
rz(-3.0806105) q[2];
rz(1.3649155) q[3];
sx q[3];
rz(-0.18326062) q[3];
sx q[3];
rz(1.1158367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0404496) q[0];
sx q[0];
rz(-1.018486) q[0];
sx q[0];
rz(-0.85969353) q[0];
rz(2.1172093) q[1];
sx q[1];
rz(-0.59737098) q[1];
sx q[1];
rz(-2.375367) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1648772) q[0];
sx q[0];
rz(-3.0857997) q[0];
sx q[0];
rz(0.68875046) q[0];
x q[1];
rz(-2.0998852) q[2];
sx q[2];
rz(-1.7211564) q[2];
sx q[2];
rz(-0.1783812) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1376201) q[1];
sx q[1];
rz(-1.8345873) q[1];
sx q[1];
rz(0.54171087) q[1];
x q[2];
rz(1.6556103) q[3];
sx q[3];
rz(-2.437371) q[3];
sx q[3];
rz(1.4798284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1102981) q[2];
sx q[2];
rz(-2.2790907) q[2];
sx q[2];
rz(-0.026570126) q[2];
rz(2.8734112) q[3];
sx q[3];
rz(-0.53332204) q[3];
sx q[3];
rz(0.69494438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37102315) q[0];
sx q[0];
rz(-0.71582782) q[0];
sx q[0];
rz(2.7243966) q[0];
rz(2.5542651) q[1];
sx q[1];
rz(-2.6081577) q[1];
sx q[1];
rz(2.3951098) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4338845) q[0];
sx q[0];
rz(-1.2217772) q[0];
sx q[0];
rz(1.2322578) q[0];
rz(-pi) q[1];
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
rz(-1.1931975) q[1];
sx q[1];
rz(-0.54497949) q[1];
rz(2.6559791) q[3];
sx q[3];
rz(-1.1385185) q[3];
sx q[3];
rz(1.9301445) q[3];
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
rz(-1.0725675) q[3];
sx q[3];
rz(-1.8662063) q[3];
sx q[3];
rz(-0.5630365) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6801341) q[0];
sx q[0];
rz(-0.38668329) q[0];
sx q[0];
rz(-0.65163809) q[0];
rz(-2.1561275) q[1];
sx q[1];
rz(-0.54254222) q[1];
sx q[1];
rz(-0.14269565) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14772077) q[0];
sx q[0];
rz(-1.8447478) q[0];
sx q[0];
rz(0.065263211) q[0];
rz(-pi) q[1];
rz(-2.8593117) q[2];
sx q[2];
rz(-2.9476894) q[2];
sx q[2];
rz(2.1197723) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0662996) q[1];
sx q[1];
rz(-0.48978862) q[1];
sx q[1];
rz(-2.879786) q[1];
rz(-pi) q[2];
rz(1.296259) q[3];
sx q[3];
rz(-2.6724788) q[3];
sx q[3];
rz(-2.9363471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5695213) q[2];
sx q[2];
rz(-2.5321952) q[2];
sx q[2];
rz(2.209668) q[2];
rz(0.43632397) q[3];
sx q[3];
rz(-0.47567979) q[3];
sx q[3];
rz(-1.1599734) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0315899) q[0];
sx q[0];
rz(-0.45781198) q[0];
sx q[0];
rz(0.97719595) q[0];
rz(-pi) q[1];
rz(1.4740491) q[2];
sx q[2];
rz(-1.821234) q[2];
sx q[2];
rz(1.4055523) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56768306) q[1];
sx q[1];
rz(-1.0956109) q[1];
sx q[1];
rz(-3.1337993) q[1];
rz(-0.55842583) q[3];
sx q[3];
rz(-1.5990725) q[3];
sx q[3];
rz(-2.0818349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5074978) q[2];
sx q[2];
rz(-1.2853421) q[2];
sx q[2];
rz(-2.0095339) q[2];
rz(-3.0105528) q[3];
sx q[3];
rz(-1.9877501) q[3];
sx q[3];
rz(1.2426144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89192724) q[0];
sx q[0];
rz(-2.6234493) q[0];
sx q[0];
rz(3.0666572) q[0];
rz(0.19649188) q[1];
sx q[1];
rz(-2.7939929) q[1];
sx q[1];
rz(-2.4620655) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20048902) q[0];
sx q[0];
rz(-1.7723494) q[0];
sx q[0];
rz(-2.5747712) q[0];
x q[1];
rz(-0.029506186) q[2];
sx q[2];
rz(-1.974612) q[2];
sx q[2];
rz(2.4798415) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30057762) q[1];
sx q[1];
rz(-1.1493582) q[1];
sx q[1];
rz(-0.31068641) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5854168) q[3];
sx q[3];
rz(-2.0818266) q[3];
sx q[3];
rz(2.8255812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3139265) q[2];
sx q[2];
rz(-0.9845261) q[2];
sx q[2];
rz(-1.0796245) q[2];
rz(-0.45352724) q[3];
sx q[3];
rz(-0.87117666) q[3];
sx q[3];
rz(2.1760904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.6249348) q[0];
sx q[0];
rz(-0.17913945) q[0];
sx q[0];
rz(-0.13570413) q[0];
rz(-2.9756359) q[1];
sx q[1];
rz(-2.143492) q[1];
sx q[1];
rz(-0.96806324) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5639599) q[0];
sx q[0];
rz(-1.5819307) q[0];
sx q[0];
rz(-2.0861113) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8161681) q[2];
sx q[2];
rz(-0.55578631) q[2];
sx q[2];
rz(1.1664558) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.28169808) q[1];
sx q[1];
rz(-1.3602232) q[1];
sx q[1];
rz(-0.50986664) q[1];
rz(-pi) q[2];
rz(1.8289205) q[3];
sx q[3];
rz(-1.934623) q[3];
sx q[3];
rz(2.8175333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.82219899) q[2];
sx q[2];
rz(-0.92331702) q[2];
sx q[2];
rz(-2.609002) q[2];
rz(-2.3432664) q[3];
sx q[3];
rz(-2.7440378) q[3];
sx q[3];
rz(0.09440162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.7867197) q[0];
sx q[0];
rz(-0.69364554) q[0];
sx q[0];
rz(2.4548446) q[0];
rz(-1.2314388) q[1];
sx q[1];
rz(-1.8309007) q[1];
sx q[1];
rz(0.062006921) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8221164) q[0];
sx q[0];
rz(-1.6363167) q[0];
sx q[0];
rz(-2.4358552) q[0];
rz(1.4976296) q[2];
sx q[2];
rz(-1.760963) q[2];
sx q[2];
rz(-0.41997465) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9888284) q[1];
sx q[1];
rz(-1.7159492) q[1];
sx q[1];
rz(-1.3671021) q[1];
rz(-pi) q[2];
rz(2.8436321) q[3];
sx q[3];
rz(-1.584867) q[3];
sx q[3];
rz(1.640682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9318781) q[2];
sx q[2];
rz(-0.96499062) q[2];
sx q[2];
rz(2.9689201) q[2];
rz(0.73575819) q[3];
sx q[3];
rz(-2.5685205) q[3];
sx q[3];
rz(2.6652523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19685766) q[0];
sx q[0];
rz(-1.6898962) q[0];
sx q[0];
rz(1.663399) q[0];
rz(2.2810777) q[1];
sx q[1];
rz(-1.1970701) q[1];
sx q[1];
rz(1.7477716) q[1];
rz(-0.12254006) q[2];
sx q[2];
rz(-1.4643639) q[2];
sx q[2];
rz(-0.53565462) q[2];
rz(-2.6113137) q[3];
sx q[3];
rz(-2.1580631) q[3];
sx q[3];
rz(-2.8736339) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
