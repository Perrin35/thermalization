OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.175617) q[0];
sx q[0];
rz(4.7369851) q[0];
sx q[0];
rz(10.935258) q[0];
rz(1.0881967) q[1];
sx q[1];
rz(-1.3047941) q[1];
sx q[1];
rz(2.3545177) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0233106) q[0];
sx q[0];
rz(-2.2459789) q[0];
sx q[0];
rz(-1.7143296) q[0];
x q[1];
rz(-1.6252615) q[2];
sx q[2];
rz(-1.174841) q[2];
sx q[2];
rz(1.8725644) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0144498) q[1];
sx q[1];
rz(-2.9269955) q[1];
sx q[1];
rz(-2.9693274) q[1];
rz(-0.6179469) q[3];
sx q[3];
rz(-2.0701346) q[3];
sx q[3];
rz(-1.7955347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8481019) q[2];
sx q[2];
rz(-0.11152554) q[2];
sx q[2];
rz(0.21609406) q[2];
rz(0.50348336) q[3];
sx q[3];
rz(-1.8899906) q[3];
sx q[3];
rz(3.0751198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9965432) q[0];
sx q[0];
rz(-1.1798877) q[0];
sx q[0];
rz(0.37522069) q[0];
rz(-2.5253865) q[1];
sx q[1];
rz(-2.1245978) q[1];
sx q[1];
rz(-0.18888758) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.045563) q[0];
sx q[0];
rz(-0.83635533) q[0];
sx q[0];
rz(-2.1406853) q[0];
x q[1];
rz(-2.6134256) q[2];
sx q[2];
rz(-0.91120992) q[2];
sx q[2];
rz(0.41925493) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9465883) q[1];
sx q[1];
rz(-2.431181) q[1];
sx q[1];
rz(1.2571774) q[1];
x q[2];
rz(2.5152605) q[3];
sx q[3];
rz(-2.5181558) q[3];
sx q[3];
rz(1.745519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6836493) q[2];
sx q[2];
rz(-1.318149) q[2];
sx q[2];
rz(1.9981492) q[2];
rz(-0.6360561) q[3];
sx q[3];
rz(-0.053074107) q[3];
sx q[3];
rz(1.599357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77883333) q[0];
sx q[0];
rz(-0.23985671) q[0];
sx q[0];
rz(1.973961) q[0];
rz(0.12655839) q[1];
sx q[1];
rz(-1.5153706) q[1];
sx q[1];
rz(0.57356858) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9605947) q[0];
sx q[0];
rz(-0.74686909) q[0];
sx q[0];
rz(-1.4372212) q[0];
rz(1.838348) q[2];
sx q[2];
rz(-1.9256136) q[2];
sx q[2];
rz(1.6521542) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3237649) q[1];
sx q[1];
rz(-0.43356248) q[1];
sx q[1];
rz(1.6852412) q[1];
rz(-pi) q[2];
rz(-0.54053791) q[3];
sx q[3];
rz(-0.73523587) q[3];
sx q[3];
rz(1.0395887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6557189) q[2];
sx q[2];
rz(-1.8334917) q[2];
sx q[2];
rz(-1.5029079) q[2];
rz(-1.2949299) q[3];
sx q[3];
rz(-0.27532268) q[3];
sx q[3];
rz(-0.82726971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5591705) q[0];
sx q[0];
rz(-0.14635135) q[0];
sx q[0];
rz(-0.91648066) q[0];
rz(0.16042635) q[1];
sx q[1];
rz(-1.285099) q[1];
sx q[1];
rz(2.4442576) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9458185) q[0];
sx q[0];
rz(-1.2504638) q[0];
sx q[0];
rz(2.8813502) q[0];
rz(-2.4795697) q[2];
sx q[2];
rz(-1.3274091) q[2];
sx q[2];
rz(-3.0070729) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5494255) q[1];
sx q[1];
rz(-0.87199713) q[1];
sx q[1];
rz(-2.0023842) q[1];
rz(-1.55682) q[3];
sx q[3];
rz(-1.7009354) q[3];
sx q[3];
rz(-1.1443477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9852898) q[2];
sx q[2];
rz(-2.6142945) q[2];
sx q[2];
rz(-1.1227597) q[2];
rz(-1.5491693) q[3];
sx q[3];
rz(-1.4607818) q[3];
sx q[3];
rz(2.7482225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8321946) q[0];
sx q[0];
rz(-3.0351312) q[0];
sx q[0];
rz(1.1291809) q[0];
rz(-0.87942266) q[1];
sx q[1];
rz(-1.8303454) q[1];
sx q[1];
rz(2.5130491) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56282069) q[0];
sx q[0];
rz(-1.4537174) q[0];
sx q[0];
rz(-1.4463918) q[0];
x q[1];
rz(-2.1206003) q[2];
sx q[2];
rz(-2.1537184) q[2];
sx q[2];
rz(1.1958808) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4820833) q[1];
sx q[1];
rz(-1.0190367) q[1];
sx q[1];
rz(3.0943046) q[1];
rz(-0.11793187) q[3];
sx q[3];
rz(-2.3068256) q[3];
sx q[3];
rz(2.6835172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2200372) q[2];
sx q[2];
rz(-2.3735789) q[2];
sx q[2];
rz(-1.7680232) q[2];
rz(0.31275648) q[3];
sx q[3];
rz(-2.0650568) q[3];
sx q[3];
rz(0.086418644) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3304928) q[0];
sx q[0];
rz(-2.436315) q[0];
sx q[0];
rz(-2.9789441) q[0];
rz(-1.2341518) q[1];
sx q[1];
rz(-1.1470497) q[1];
sx q[1];
rz(-2.2208234) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7013848) q[0];
sx q[0];
rz(-0.21344859) q[0];
sx q[0];
rz(1.0769072) q[0];
rz(-pi) q[1];
rz(0.92451325) q[2];
sx q[2];
rz(-1.1072031) q[2];
sx q[2];
rz(1.6789503) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.858135) q[1];
sx q[1];
rz(-0.75443903) q[1];
sx q[1];
rz(-3.0467176) q[1];
rz(1.7397171) q[3];
sx q[3];
rz(-0.79577765) q[3];
sx q[3];
rz(-0.9267048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.84772253) q[2];
sx q[2];
rz(-2.4203478) q[2];
sx q[2];
rz(-3.0976307) q[2];
rz(1.7196767) q[3];
sx q[3];
rz(-0.43216643) q[3];
sx q[3];
rz(3.1066084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1304355) q[0];
sx q[0];
rz(-2.3373117) q[0];
sx q[0];
rz(-0.093611896) q[0];
rz(1.0253819) q[1];
sx q[1];
rz(-1.5954115) q[1];
sx q[1];
rz(2.3536918) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2219395) q[0];
sx q[0];
rz(-1.0353695) q[0];
sx q[0];
rz(0.8485283) q[0];
rz(-pi) q[1];
rz(-1.7406932) q[2];
sx q[2];
rz(-2.133103) q[2];
sx q[2];
rz(3.135709) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8389143) q[1];
sx q[1];
rz(-0.60015772) q[1];
sx q[1];
rz(0.072235302) q[1];
rz(3.0209846) q[3];
sx q[3];
rz(-1.7123619) q[3];
sx q[3];
rz(-2.0742576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0660144) q[2];
sx q[2];
rz(-2.5983073) q[2];
sx q[2];
rz(-0.96599609) q[2];
rz(-0.14723369) q[3];
sx q[3];
rz(-1.6817663) q[3];
sx q[3];
rz(0.60432965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0886993) q[0];
sx q[0];
rz(-1.7558782) q[0];
sx q[0];
rz(-1.8889486) q[0];
rz(3.0839651) q[1];
sx q[1];
rz(-1.3725955) q[1];
sx q[1];
rz(2.7431814) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1712043) q[0];
sx q[0];
rz(-2.1792767) q[0];
sx q[0];
rz(1.7826016) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4008825) q[2];
sx q[2];
rz(-1.0771829) q[2];
sx q[2];
rz(2.5371671) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0115464) q[1];
sx q[1];
rz(-0.1999487) q[1];
sx q[1];
rz(1.5910831) q[1];
x q[2];
rz(0.1755821) q[3];
sx q[3];
rz(-0.55003234) q[3];
sx q[3];
rz(1.8635141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.56732279) q[2];
sx q[2];
rz(-1.2719354) q[2];
sx q[2];
rz(0.92180139) q[2];
rz(-2.4343893) q[3];
sx q[3];
rz(-2.0965529) q[3];
sx q[3];
rz(0.31064492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1929753) q[0];
sx q[0];
rz(-0.12117584) q[0];
sx q[0];
rz(2.2344672) q[0];
rz(2.6616197) q[1];
sx q[1];
rz(-2.0860806) q[1];
sx q[1];
rz(2.5352246) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7146695) q[0];
sx q[0];
rz(-1.4703317) q[0];
sx q[0];
rz(2.8313374) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0125514) q[2];
sx q[2];
rz(-2.0336359) q[2];
sx q[2];
rz(1.7985345) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1923101) q[1];
sx q[1];
rz(-0.81337608) q[1];
sx q[1];
rz(-1.9793012) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4828047) q[3];
sx q[3];
rz(-1.6571991) q[3];
sx q[3];
rz(-0.98443809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2839462) q[2];
sx q[2];
rz(-1.5479156) q[2];
sx q[2];
rz(-0.66961163) q[2];
rz(-0.56704632) q[3];
sx q[3];
rz(-2.0958869) q[3];
sx q[3];
rz(-1.4001575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6794716) q[0];
sx q[0];
rz(-1.3422048) q[0];
sx q[0];
rz(-0.8959499) q[0];
rz(-0.040291928) q[1];
sx q[1];
rz(-0.94172421) q[1];
sx q[1];
rz(-1.5774073) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5524254) q[0];
sx q[0];
rz(-1.5660813) q[0];
sx q[0];
rz(1.5472359) q[0];
rz(2.4559458) q[2];
sx q[2];
rz(-1.2276936) q[2];
sx q[2];
rz(2.1941663) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2103232) q[1];
sx q[1];
rz(-1.0627373) q[1];
sx q[1];
rz(2.6711051) q[1];
x q[2];
rz(0.96503105) q[3];
sx q[3];
rz(-1.26767) q[3];
sx q[3];
rz(-0.81699634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8335874) q[2];
sx q[2];
rz(-2.8218125) q[2];
sx q[2];
rz(1.3099111) q[2];
rz(-1.1854019) q[3];
sx q[3];
rz(-0.88806051) q[3];
sx q[3];
rz(-1.6185224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1401405) q[0];
sx q[0];
rz(-1.8403213) q[0];
sx q[0];
rz(2.1885827) q[0];
rz(-0.078527191) q[1];
sx q[1];
rz(-2.0546866) q[1];
sx q[1];
rz(2.3436117) q[1];
rz(-2.8817202) q[2];
sx q[2];
rz(-2.2287579) q[2];
sx q[2];
rz(-1.8884434) q[2];
rz(-2.4918741) q[3];
sx q[3];
rz(-2.8856554) q[3];
sx q[3];
rz(1.9323991) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
