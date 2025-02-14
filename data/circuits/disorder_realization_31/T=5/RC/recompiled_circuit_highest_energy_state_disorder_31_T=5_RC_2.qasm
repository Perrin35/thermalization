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
rz(1.2827058) q[0];
sx q[0];
rz(-2.1252706) q[0];
sx q[0];
rz(-1.2555726) q[0];
rz(-1.3610871) q[1];
sx q[1];
rz(-0.95870107) q[1];
sx q[1];
rz(0.60671848) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5131305) q[0];
sx q[0];
rz(-0.89848041) q[0];
sx q[0];
rz(-2.1450514) q[0];
rz(0.032564596) q[2];
sx q[2];
rz(-1.3059907) q[2];
sx q[2];
rz(-0.84003583) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7156034) q[1];
sx q[1];
rz(-1.2066325) q[1];
sx q[1];
rz(0.13873546) q[1];
x q[2];
rz(-1.4503202) q[3];
sx q[3];
rz(-1.7306149) q[3];
sx q[3];
rz(2.8152856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9386193) q[2];
sx q[2];
rz(-1.9846658) q[2];
sx q[2];
rz(2.970676) q[2];
rz(1.1275229) q[3];
sx q[3];
rz(-2.5850962) q[3];
sx q[3];
rz(-3.0881622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.663986) q[0];
sx q[0];
rz(-2.8150616) q[0];
sx q[0];
rz(-2.4531181) q[0];
rz(1.3779878) q[1];
sx q[1];
rz(-2.1207899) q[1];
sx q[1];
rz(0.14150208) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6287891) q[0];
sx q[0];
rz(-2.4806285) q[0];
sx q[0];
rz(2.4557607) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9269518) q[2];
sx q[2];
rz(-2.7032489) q[2];
sx q[2];
rz(0.12210309) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5192109) q[1];
sx q[1];
rz(-1.8528717) q[1];
sx q[1];
rz(-0.36954986) q[1];
rz(1.1859975) q[3];
sx q[3];
rz(-1.4963004) q[3];
sx q[3];
rz(2.7423679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3839533) q[2];
sx q[2];
rz(-2.6944104) q[2];
sx q[2];
rz(-3.1128856) q[2];
rz(2.7590397) q[3];
sx q[3];
rz(-1.5111204) q[3];
sx q[3];
rz(2.9401275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89556995) q[0];
sx q[0];
rz(-1.0423132) q[0];
sx q[0];
rz(0.37164715) q[0];
rz(2.9314575) q[1];
sx q[1];
rz(-2.7919283) q[1];
sx q[1];
rz(1.9336611) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4233404) q[0];
sx q[0];
rz(-0.99159843) q[0];
sx q[0];
rz(2.1787326) q[0];
rz(-pi) q[1];
x q[1];
rz(0.082449989) q[2];
sx q[2];
rz(-1.5497314) q[2];
sx q[2];
rz(0.17848554) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.722586) q[1];
sx q[1];
rz(-1.0204781) q[1];
sx q[1];
rz(0.87700486) q[1];
rz(1.6116834) q[3];
sx q[3];
rz(-1.2256116) q[3];
sx q[3];
rz(0.57905771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.62260425) q[2];
sx q[2];
rz(-0.092978803) q[2];
sx q[2];
rz(-0.63817111) q[2];
rz(-2.7291164) q[3];
sx q[3];
rz(-1.4700438) q[3];
sx q[3];
rz(-2.0392058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0647122) q[0];
sx q[0];
rz(-0.91016722) q[0];
sx q[0];
rz(1.8338715) q[0];
rz(-2.4619596) q[1];
sx q[1];
rz(-0.72703397) q[1];
sx q[1];
rz(1.9576498) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7789202) q[0];
sx q[0];
rz(-2.2244172) q[0];
sx q[0];
rz(1.9833167) q[0];
x q[1];
rz(-0.45071256) q[2];
sx q[2];
rz(-0.26077429) q[2];
sx q[2];
rz(-0.63577543) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7412926) q[1];
sx q[1];
rz(-1.502599) q[1];
sx q[1];
rz(0.3505716) q[1];
rz(0.13069066) q[3];
sx q[3];
rz(-1.4447235) q[3];
sx q[3];
rz(-3.0265615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1188941) q[2];
sx q[2];
rz(-2.9939632) q[2];
sx q[2];
rz(-1.0559319) q[2];
rz(2.8583156) q[3];
sx q[3];
rz(-1.8746459) q[3];
sx q[3];
rz(-1.7578846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.088575514) q[0];
sx q[0];
rz(-0.96240369) q[0];
sx q[0];
rz(1.4085294) q[0];
rz(-0.22964302) q[1];
sx q[1];
rz(-2.2848928) q[1];
sx q[1];
rz(-1.9502669) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9339122) q[0];
sx q[0];
rz(-1.7125074) q[0];
sx q[0];
rz(1.974154) q[0];
rz(-pi) q[1];
rz(1.8764087) q[2];
sx q[2];
rz(-1.6907881) q[2];
sx q[2];
rz(-0.12745276) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7745668) q[1];
sx q[1];
rz(-1.3210249) q[1];
sx q[1];
rz(2.7109409) q[1];
rz(-2.5727083) q[3];
sx q[3];
rz(-2.1422022) q[3];
sx q[3];
rz(2.6517154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4952937) q[2];
sx q[2];
rz(-2.7806492) q[2];
sx q[2];
rz(-1.6634644) q[2];
rz(0.16119257) q[3];
sx q[3];
rz(-1.5321923) q[3];
sx q[3];
rz(-0.78981367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6394871) q[0];
sx q[0];
rz(-0.4244856) q[0];
sx q[0];
rz(-0.91833997) q[0];
rz(0.25686747) q[1];
sx q[1];
rz(-2.145642) q[1];
sx q[1];
rz(-1.1161944) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6300121) q[0];
sx q[0];
rz(-2.0697547) q[0];
sx q[0];
rz(-0.60006882) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6565422) q[2];
sx q[2];
rz(-1.8433136) q[2];
sx q[2];
rz(2.9086824) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4914272) q[1];
sx q[1];
rz(-2.8313447) q[1];
sx q[1];
rz(0.68474557) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5258216) q[3];
sx q[3];
rz(-2.8463461) q[3];
sx q[3];
rz(-2.8930882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9548698) q[2];
sx q[2];
rz(-2.335151) q[2];
sx q[2];
rz(0.40697971) q[2];
rz(2.5278029) q[3];
sx q[3];
rz(-1.5301306) q[3];
sx q[3];
rz(2.6500402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75519049) q[0];
sx q[0];
rz(-0.84119216) q[0];
sx q[0];
rz(0.92591539) q[0];
rz(0.47261247) q[1];
sx q[1];
rz(-2.0874529) q[1];
sx q[1];
rz(-2.6714163) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9191943) q[0];
sx q[0];
rz(-0.83756768) q[0];
sx q[0];
rz(-0.11239692) q[0];
rz(-pi) q[1];
rz(2.0704001) q[2];
sx q[2];
rz(-2.4322369) q[2];
sx q[2];
rz(1.6215404) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9035569) q[1];
sx q[1];
rz(-0.78552526) q[1];
sx q[1];
rz(0.051427917) q[1];
rz(0.84104611) q[3];
sx q[3];
rz(-1.5980835) q[3];
sx q[3];
rz(-0.59852615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4818695) q[2];
sx q[2];
rz(-1.0708151) q[2];
sx q[2];
rz(-2.3107963) q[2];
rz(0.080549084) q[3];
sx q[3];
rz(-1.0815257) q[3];
sx q[3];
rz(-2.1848333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.047121) q[0];
sx q[0];
rz(-1.6166649) q[0];
sx q[0];
rz(2.9834874) q[0];
rz(-0.72599167) q[1];
sx q[1];
rz(-2.7114365) q[1];
sx q[1];
rz(0.37011883) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67664941) q[0];
sx q[0];
rz(-1.2386564) q[0];
sx q[0];
rz(-1.4737096) q[0];
x q[1];
rz(-2.3838777) q[2];
sx q[2];
rz(-1.9279459) q[2];
sx q[2];
rz(2.0983608) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2850879) q[1];
sx q[1];
rz(-2.5743432) q[1];
sx q[1];
rz(0.7817661) q[1];
rz(0.94929578) q[3];
sx q[3];
rz(-2.5007476) q[3];
sx q[3];
rz(2.8261938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.98081723) q[2];
sx q[2];
rz(-1.846401) q[2];
sx q[2];
rz(3.004461) q[2];
rz(-1.6454227) q[3];
sx q[3];
rz(-0.6936332) q[3];
sx q[3];
rz(-2.6579198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5982323) q[0];
sx q[0];
rz(-1.0022663) q[0];
sx q[0];
rz(-0.12635669) q[0];
rz(-0.57451808) q[1];
sx q[1];
rz(-2.2210329) q[1];
sx q[1];
rz(-0.7472907) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80602801) q[0];
sx q[0];
rz(-1.757041) q[0];
sx q[0];
rz(1.1769017) q[0];
x q[1];
rz(-0.66102513) q[2];
sx q[2];
rz(-0.82023652) q[2];
sx q[2];
rz(-0.30414061) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9979228) q[1];
sx q[1];
rz(-1.9410681) q[1];
sx q[1];
rz(-1.0940432) q[1];
x q[2];
rz(-2.8793192) q[3];
sx q[3];
rz(-1.6991391) q[3];
sx q[3];
rz(-0.94576361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1025461) q[2];
sx q[2];
rz(-2.0206082) q[2];
sx q[2];
rz(2.5940564) q[2];
rz(-1.3012137) q[3];
sx q[3];
rz(-1.1373212) q[3];
sx q[3];
rz(3.014452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.309677) q[0];
sx q[0];
rz(-1.3818106) q[0];
sx q[0];
rz(2.5610913) q[0];
rz(2.8864587) q[1];
sx q[1];
rz(-0.83108035) q[1];
sx q[1];
rz(1.9560122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0609293) q[0];
sx q[0];
rz(-0.97968819) q[0];
sx q[0];
rz(0.76158686) q[0];
x q[1];
rz(-1.2078778) q[2];
sx q[2];
rz(-1.5722154) q[2];
sx q[2];
rz(-3.0848173) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.23981389) q[1];
sx q[1];
rz(-2.0396775) q[1];
sx q[1];
rz(-1.6704876) q[1];
rz(1.8099635) q[3];
sx q[3];
rz(-0.73430919) q[3];
sx q[3];
rz(1.7074761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.4127976) q[2];
sx q[2];
rz(-2.2788861) q[2];
sx q[2];
rz(2.1346788) q[2];
rz(-1.5236731) q[3];
sx q[3];
rz(-0.98098743) q[3];
sx q[3];
rz(1.1716918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11378743) q[0];
sx q[0];
rz(-1.8814977) q[0];
sx q[0];
rz(0.30497288) q[0];
rz(-1.2420568) q[1];
sx q[1];
rz(-1.8546974) q[1];
sx q[1];
rz(0.7484662) q[1];
rz(-2.9681846) q[2];
sx q[2];
rz(-1.1447403) q[2];
sx q[2];
rz(-0.8131342) q[2];
rz(0.59769825) q[3];
sx q[3];
rz(-2.2168474) q[3];
sx q[3];
rz(-0.26099152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
