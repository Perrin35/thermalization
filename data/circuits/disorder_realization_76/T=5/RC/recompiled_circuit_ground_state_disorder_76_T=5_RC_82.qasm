OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6948833) q[0];
sx q[0];
rz(-1.0816242) q[0];
sx q[0];
rz(2.4021436) q[0];
rz(2.3820355) q[1];
sx q[1];
rz(-1.3146725) q[1];
sx q[1];
rz(-1.7159599) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65042881) q[0];
sx q[0];
rz(-0.61507498) q[0];
sx q[0];
rz(1.0577202) q[0];
rz(-2.6665568) q[2];
sx q[2];
rz(-2.5940707) q[2];
sx q[2];
rz(0.56625596) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1358429) q[1];
sx q[1];
rz(-1.4228361) q[1];
sx q[1];
rz(-3.0632076) q[1];
rz(-pi) q[2];
rz(1.6207148) q[3];
sx q[3];
rz(-0.71943362) q[3];
sx q[3];
rz(2.0572544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.90929675) q[2];
sx q[2];
rz(-0.84315073) q[2];
sx q[2];
rz(-0.24093957) q[2];
rz(3.1124034) q[3];
sx q[3];
rz(-1.3386644) q[3];
sx q[3];
rz(1.1791112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643243) q[0];
sx q[0];
rz(-1.9376396) q[0];
sx q[0];
rz(0.48467317) q[0];
rz(-2.4618705) q[1];
sx q[1];
rz(-1.8599963) q[1];
sx q[1];
rz(1.1601123) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74014464) q[0];
sx q[0];
rz(-0.30456802) q[0];
sx q[0];
rz(-0.061621678) q[0];
rz(-pi) q[1];
rz(-0.71573513) q[2];
sx q[2];
rz(-1.01075) q[2];
sx q[2];
rz(2.3079688) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.45534652) q[1];
sx q[1];
rz(-0.97963154) q[1];
sx q[1];
rz(-1.8704198) q[1];
rz(-0.5228225) q[3];
sx q[3];
rz(-1.2827875) q[3];
sx q[3];
rz(0.097974591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.29740563) q[2];
sx q[2];
rz(-2.83941) q[2];
sx q[2];
rz(-0.45787946) q[2];
rz(1.1229905) q[3];
sx q[3];
rz(-2.1419958) q[3];
sx q[3];
rz(-0.56639731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8726525) q[0];
sx q[0];
rz(-0.70916969) q[0];
sx q[0];
rz(-0.031524468) q[0];
rz(2.8541376) q[1];
sx q[1];
rz(-2.2645686) q[1];
sx q[1];
rz(1.9256176) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8376802) q[0];
sx q[0];
rz(-1.0420351) q[0];
sx q[0];
rz(-0.75334511) q[0];
rz(-0.38162614) q[2];
sx q[2];
rz(-0.51593057) q[2];
sx q[2];
rz(0.53781539) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.029473631) q[1];
sx q[1];
rz(-0.39327565) q[1];
sx q[1];
rz(0.78177364) q[1];
rz(-2.3787556) q[3];
sx q[3];
rz(-0.82218542) q[3];
sx q[3];
rz(-1.1824106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.03881255) q[2];
sx q[2];
rz(-1.5993885) q[2];
sx q[2];
rz(0.83596027) q[2];
rz(1.3577667) q[3];
sx q[3];
rz(-0.72503763) q[3];
sx q[3];
rz(0.74762216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(1.8589856) q[0];
sx q[0];
rz(-1.7361807) q[0];
sx q[0];
rz(3.0614241) q[0];
rz(1.0824341) q[1];
sx q[1];
rz(-2.9083462) q[1];
sx q[1];
rz(1.8128043) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4445368) q[0];
sx q[0];
rz(-1.8955243) q[0];
sx q[0];
rz(-0.53931196) q[0];
x q[1];
rz(0.32281876) q[2];
sx q[2];
rz(-1.7344966) q[2];
sx q[2];
rz(1.9893008) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.13212407) q[1];
sx q[1];
rz(-1.4259725) q[1];
sx q[1];
rz(2.0028466) q[1];
rz(-pi) q[2];
rz(0.53502632) q[3];
sx q[3];
rz(-0.38540977) q[3];
sx q[3];
rz(-2.7148394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7666011) q[2];
sx q[2];
rz(-0.98321715) q[2];
sx q[2];
rz(0.90744606) q[2];
rz(0.67029101) q[3];
sx q[3];
rz(-2.3136316) q[3];
sx q[3];
rz(1.7214187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2403253) q[0];
sx q[0];
rz(-2.2690052) q[0];
sx q[0];
rz(-2.7101044) q[0];
rz(1.1314499) q[1];
sx q[1];
rz(-1.9624458) q[1];
sx q[1];
rz(-2.7010837) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90135306) q[0];
sx q[0];
rz(-0.47751891) q[0];
sx q[0];
rz(1.3130472) q[0];
rz(1.9856057) q[2];
sx q[2];
rz(-1.7602663) q[2];
sx q[2];
rz(-2.1304325) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7311598) q[1];
sx q[1];
rz(-1.9010547) q[1];
sx q[1];
rz(-1.6124658) q[1];
rz(-pi) q[2];
rz(2.2045361) q[3];
sx q[3];
rz(-2.0645294) q[3];
sx q[3];
rz(-2.1449094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0195007) q[2];
sx q[2];
rz(-0.92635218) q[2];
sx q[2];
rz(2.7276373) q[2];
rz(1.7533938) q[3];
sx q[3];
rz(-1.5940758) q[3];
sx q[3];
rz(3.0164914) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51467657) q[0];
sx q[0];
rz(-1.7358945) q[0];
sx q[0];
rz(2.2077014) q[0];
rz(-2.4712708) q[1];
sx q[1];
rz(-1.2443845) q[1];
sx q[1];
rz(-1.3353039) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1191171) q[0];
sx q[0];
rz(-2.5316105) q[0];
sx q[0];
rz(1.7458292) q[0];
x q[1];
rz(0.1576169) q[2];
sx q[2];
rz(-2.981957) q[2];
sx q[2];
rz(0.0052513382) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0508479) q[1];
sx q[1];
rz(-0.11493348) q[1];
sx q[1];
rz(1.3728218) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1298529) q[3];
sx q[3];
rz(-0.63705963) q[3];
sx q[3];
rz(0.54367346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2495217) q[2];
sx q[2];
rz(-2.8688909) q[2];
sx q[2];
rz(1.8708694) q[2];
rz(1.3696085) q[3];
sx q[3];
rz(-1.9000051) q[3];
sx q[3];
rz(0.52186596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18913604) q[0];
sx q[0];
rz(-0.58681762) q[0];
sx q[0];
rz(-2.9879046) q[0];
rz(-2.9391089) q[1];
sx q[1];
rz(-1.2362365) q[1];
sx q[1];
rz(2.1930146) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5178793) q[0];
sx q[0];
rz(-2.1800632) q[0];
sx q[0];
rz(-1.3246956) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93164753) q[2];
sx q[2];
rz(-1.8338565) q[2];
sx q[2];
rz(-1.5342889) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2413509) q[1];
sx q[1];
rz(-0.35120041) q[1];
sx q[1];
rz(-2.3022149) q[1];
x q[2];
rz(-2.2205686) q[3];
sx q[3];
rz(-1.463784) q[3];
sx q[3];
rz(-1.8821074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.004868) q[2];
sx q[2];
rz(-1.0978881) q[2];
sx q[2];
rz(-0.0014121545) q[2];
rz(-0.062156113) q[3];
sx q[3];
rz(-1.7319738) q[3];
sx q[3];
rz(0.96127659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75981265) q[0];
sx q[0];
rz(-2.3842922) q[0];
sx q[0];
rz(1.2148452) q[0];
rz(-0.75792056) q[1];
sx q[1];
rz(-0.52752033) q[1];
sx q[1];
rz(2.5835999) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0601301) q[0];
sx q[0];
rz(-1.234982) q[0];
sx q[0];
rz(-2.8804638) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55155374) q[2];
sx q[2];
rz(-1.4467014) q[2];
sx q[2];
rz(-2.6089904) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2633725) q[1];
sx q[1];
rz(-1.0910631) q[1];
sx q[1];
rz(-2.8271476) q[1];
rz(-pi) q[2];
rz(-3.1147546) q[3];
sx q[3];
rz(-1.2937814) q[3];
sx q[3];
rz(0.42624302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5672292) q[2];
sx q[2];
rz(-1.9476674) q[2];
sx q[2];
rz(0.26926678) q[2];
rz(-1.9783798) q[3];
sx q[3];
rz(-2.5389157) q[3];
sx q[3];
rz(0.24063024) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6729386) q[0];
sx q[0];
rz(-2.1042295) q[0];
sx q[0];
rz(-1.0733676) q[0];
rz(-2.0138373) q[1];
sx q[1];
rz(-2.914371) q[1];
sx q[1];
rz(-0.55666298) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4024324) q[0];
sx q[0];
rz(-0.35050979) q[0];
sx q[0];
rz(-1.1465766) q[0];
rz(0.78293856) q[2];
sx q[2];
rz(-0.73854827) q[2];
sx q[2];
rz(2.1434181) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.60815114) q[1];
sx q[1];
rz(-1.2131629) q[1];
sx q[1];
rz(-2.3945892) q[1];
rz(-pi) q[2];
rz(-2.9070517) q[3];
sx q[3];
rz(-1.4734771) q[3];
sx q[3];
rz(3.0458801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2399981) q[2];
sx q[2];
rz(-2.0549213) q[2];
sx q[2];
rz(1.979801) q[2];
rz(0.53317541) q[3];
sx q[3];
rz(-0.52459255) q[3];
sx q[3];
rz(-2.9840792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48851442) q[0];
sx q[0];
rz(-0.97452679) q[0];
sx q[0];
rz(2.5467806) q[0];
rz(-0.57890233) q[1];
sx q[1];
rz(-1.4756823) q[1];
sx q[1];
rz(2.4822809) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80583094) q[0];
sx q[0];
rz(-1.7044401) q[0];
sx q[0];
rz(-0.62778309) q[0];
x q[1];
rz(-1.1289146) q[2];
sx q[2];
rz(-1.0143447) q[2];
sx q[2];
rz(-0.2054727) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.59567) q[1];
sx q[1];
rz(-2.8575142) q[1];
sx q[1];
rz(-2.0998663) q[1];
rz(-pi) q[2];
rz(-0.26409491) q[3];
sx q[3];
rz(-1.3066829) q[3];
sx q[3];
rz(2.7434289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79002964) q[2];
sx q[2];
rz(-1.1976676) q[2];
sx q[2];
rz(-0.053038049) q[2];
rz(-1.7985581) q[3];
sx q[3];
rz(-1.8648632) q[3];
sx q[3];
rz(3.067335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85528436) q[0];
sx q[0];
rz(-1.7779779) q[0];
sx q[0];
rz(1.5386982) q[0];
rz(3.042649) q[1];
sx q[1];
rz(-1.7715441) q[1];
sx q[1];
rz(-3.0861707) q[1];
rz(0.40545736) q[2];
sx q[2];
rz(-1.5628855) q[2];
sx q[2];
rz(0.93456675) q[2];
rz(-0.74138882) q[3];
sx q[3];
rz(-1.4546053) q[3];
sx q[3];
rz(-1.4447053) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
