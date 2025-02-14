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
rz(-0.14577785) q[0];
sx q[0];
rz(-1.3480027) q[0];
sx q[0];
rz(0.3682799) q[0];
rz(-0.085973099) q[1];
sx q[1];
rz(-2.2987125) q[1];
sx q[1];
rz(0.24922961) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1220142) q[0];
sx q[0];
rz(-1.459157) q[0];
sx q[0];
rz(-2.6606002) q[0];
rz(-pi) q[1];
rz(2.1051718) q[2];
sx q[2];
rz(-1.3003361) q[2];
sx q[2];
rz(-2.0024744) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6492086) q[1];
sx q[1];
rz(-1.3351775) q[1];
sx q[1];
rz(-2.431303) q[1];
rz(-pi) q[2];
rz(2.8461012) q[3];
sx q[3];
rz(-0.37335941) q[3];
sx q[3];
rz(2.1779324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2117846) q[2];
sx q[2];
rz(-2.4346508) q[2];
sx q[2];
rz(2.9014273) q[2];
rz(0.54720488) q[3];
sx q[3];
rz(-1.5144843) q[3];
sx q[3];
rz(2.1521294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7996247) q[0];
sx q[0];
rz(-2.7662321) q[0];
sx q[0];
rz(2.4976835) q[0];
rz(-2.0945235) q[1];
sx q[1];
rz(-1.1771076) q[1];
sx q[1];
rz(-2.8783096) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1564668) q[0];
sx q[0];
rz(-1.5432285) q[0];
sx q[0];
rz(-3.1113202) q[0];
x q[1];
rz(1.5059169) q[2];
sx q[2];
rz(-1.5047964) q[2];
sx q[2];
rz(-0.22021107) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.66052478) q[1];
sx q[1];
rz(-2.1612289) q[1];
sx q[1];
rz(0.51312889) q[1];
rz(-pi) q[2];
rz(1.0399148) q[3];
sx q[3];
rz(-0.71456075) q[3];
sx q[3];
rz(-1.2604063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6836267) q[2];
sx q[2];
rz(-1.3764952) q[2];
sx q[2];
rz(-0.72570428) q[2];
rz(-0.30722412) q[3];
sx q[3];
rz(-1.8351464) q[3];
sx q[3];
rz(2.0871625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3211408) q[0];
sx q[0];
rz(-1.4844002) q[0];
sx q[0];
rz(-0.76903525) q[0];
rz(-0.10737315) q[1];
sx q[1];
rz(-1.4391876) q[1];
sx q[1];
rz(0.69951397) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0663129) q[0];
sx q[0];
rz(-1.9121721) q[0];
sx q[0];
rz(0.34324788) q[0];
rz(-pi) q[1];
rz(3.1351481) q[2];
sx q[2];
rz(-0.61476189) q[2];
sx q[2];
rz(2.4561735) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.67058439) q[1];
sx q[1];
rz(-2.0194897) q[1];
sx q[1];
rz(3.0506794) q[1];
rz(-pi) q[2];
rz(-2.3885257) q[3];
sx q[3];
rz(-2.1478466) q[3];
sx q[3];
rz(2.6073539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7677782) q[2];
sx q[2];
rz(-0.57743293) q[2];
sx q[2];
rz(2.2334297) q[2];
rz(-2.5670037) q[3];
sx q[3];
rz(-1.8879954) q[3];
sx q[3];
rz(-2.7558034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3278367) q[0];
sx q[0];
rz(-1.8151374) q[0];
sx q[0];
rz(-0.5089708) q[0];
rz(-0.058188997) q[1];
sx q[1];
rz(-1.6845082) q[1];
sx q[1];
rz(-1.0675272) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4315486) q[0];
sx q[0];
rz(-2.2674034) q[0];
sx q[0];
rz(-1.0846653) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3994965) q[2];
sx q[2];
rz(-2.7242085) q[2];
sx q[2];
rz(0.69533747) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0587973) q[1];
sx q[1];
rz(-1.4230295) q[1];
sx q[1];
rz(-3.0331103) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8640609) q[3];
sx q[3];
rz(-1.505027) q[3];
sx q[3];
rz(-1.6848967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.031614583) q[2];
sx q[2];
rz(-1.5974533) q[2];
sx q[2];
rz(0.99008375) q[2];
rz(1.3958942) q[3];
sx q[3];
rz(-3.0315704) q[3];
sx q[3];
rz(-1.8949738) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9083967) q[0];
sx q[0];
rz(-1.7447423) q[0];
sx q[0];
rz(-2.0821849) q[0];
rz(1.1147095) q[1];
sx q[1];
rz(-1.6672986) q[1];
sx q[1];
rz(1.4310736) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16742584) q[0];
sx q[0];
rz(-1.7210809) q[0];
sx q[0];
rz(-1.7420064) q[0];
x q[1];
rz(2.8818508) q[2];
sx q[2];
rz(-2.0613823) q[2];
sx q[2];
rz(0.39133137) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3472673) q[1];
sx q[1];
rz(-1.8121983) q[1];
sx q[1];
rz(2.4893087) q[1];
x q[2];
rz(-2.9274213) q[3];
sx q[3];
rz(-1.7167544) q[3];
sx q[3];
rz(2.0623061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7913671) q[2];
sx q[2];
rz(-0.75608772) q[2];
sx q[2];
rz(1.6737326) q[2];
rz(-1.0427467) q[3];
sx q[3];
rz(-1.0576893) q[3];
sx q[3];
rz(-2.3500672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4333711) q[0];
sx q[0];
rz(-1.2935761) q[0];
sx q[0];
rz(2.9300387) q[0];
rz(-0.64282203) q[1];
sx q[1];
rz(-2.0170409) q[1];
sx q[1];
rz(1.0483673) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7599568) q[0];
sx q[0];
rz(-2.4279729) q[0];
sx q[0];
rz(-0.19605331) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6238113) q[2];
sx q[2];
rz(-1.3919733) q[2];
sx q[2];
rz(-0.12763466) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0901186) q[1];
sx q[1];
rz(-1.5940545) q[1];
sx q[1];
rz(0.6384797) q[1];
x q[2];
rz(-1.6392751) q[3];
sx q[3];
rz(-2.5279495) q[3];
sx q[3];
rz(-0.63921038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.50531498) q[2];
sx q[2];
rz(-1.9376398) q[2];
sx q[2];
rz(-2.8875276) q[2];
rz(0.17357477) q[3];
sx q[3];
rz(-0.55145276) q[3];
sx q[3];
rz(1.180163) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5467065) q[0];
sx q[0];
rz(-2.7017024) q[0];
sx q[0];
rz(-2.34483) q[0];
rz(0.88675371) q[1];
sx q[1];
rz(-1.3462857) q[1];
sx q[1];
rz(0.60060445) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2830505) q[0];
sx q[0];
rz(-1.5158537) q[0];
sx q[0];
rz(2.8514329) q[0];
rz(2.372416) q[2];
sx q[2];
rz(-1.4750136) q[2];
sx q[2];
rz(-2.8623919) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6930042) q[1];
sx q[1];
rz(-1.4088165) q[1];
sx q[1];
rz(0.88108351) q[1];
rz(-pi) q[2];
rz(-2.10403) q[3];
sx q[3];
rz(-1.2640461) q[3];
sx q[3];
rz(2.3045134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.1116921) q[2];
sx q[2];
rz(-1.7039958) q[2];
sx q[2];
rz(2.0491484) q[2];
rz(1.368329) q[3];
sx q[3];
rz(-0.60951257) q[3];
sx q[3];
rz(0.4121367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71548897) q[0];
sx q[0];
rz(-2.0273209) q[0];
sx q[0];
rz(0.45147595) q[0];
rz(0.077797912) q[1];
sx q[1];
rz(-2.1092238) q[1];
sx q[1];
rz(-2.480004) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63115727) q[0];
sx q[0];
rz(-1.1589739) q[0];
sx q[0];
rz(-0.23958448) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0691452) q[2];
sx q[2];
rz(-1.632164) q[2];
sx q[2];
rz(-0.59624404) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.50345031) q[1];
sx q[1];
rz(-1.7026445) q[1];
sx q[1];
rz(1.1080469) q[1];
x q[2];
rz(-2.745146) q[3];
sx q[3];
rz(-1.8155193) q[3];
sx q[3];
rz(0.94371599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0563858) q[2];
sx q[2];
rz(-1.8737917) q[2];
sx q[2];
rz(-2.3351604) q[2];
rz(1.2422397) q[3];
sx q[3];
rz(-2.9759088) q[3];
sx q[3];
rz(0.14036673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(0.30271444) q[0];
sx q[0];
rz(-1.1352204) q[0];
sx q[0];
rz(0.43933991) q[0];
rz(-1.9413403) q[1];
sx q[1];
rz(-2.3019583) q[1];
sx q[1];
rz(-0.95019788) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.882975) q[0];
sx q[0];
rz(-1.1990917) q[0];
sx q[0];
rz(-2.2540493) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0884041) q[2];
sx q[2];
rz(-2.8206499) q[2];
sx q[2];
rz(-2.1118922) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.31017329) q[1];
sx q[1];
rz(-1.2155639) q[1];
sx q[1];
rz(-0.11174496) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1113273) q[3];
sx q[3];
rz(-1.221962) q[3];
sx q[3];
rz(-0.23603786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3085559) q[2];
sx q[2];
rz(-0.50806442) q[2];
sx q[2];
rz(1.1888095) q[2];
rz(0.10558072) q[3];
sx q[3];
rz(-0.9413541) q[3];
sx q[3];
rz(2.1281435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078876615) q[0];
sx q[0];
rz(-0.3083516) q[0];
sx q[0];
rz(2.4421332) q[0];
rz(-1.7309011) q[1];
sx q[1];
rz(-0.78838333) q[1];
sx q[1];
rz(-2.9404822) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58696207) q[0];
sx q[0];
rz(-2.2788958) q[0];
sx q[0];
rz(2.5434341) q[0];
x q[1];
rz(-2.6568065) q[2];
sx q[2];
rz(-2.1109627) q[2];
sx q[2];
rz(2.6933984) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2664908) q[1];
sx q[1];
rz(-1.3492341) q[1];
sx q[1];
rz(0.87163493) q[1];
x q[2];
rz(-2.9617025) q[3];
sx q[3];
rz(-1.4818061) q[3];
sx q[3];
rz(-0.51489553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9451399) q[2];
sx q[2];
rz(-1.8818734) q[2];
sx q[2];
rz(-3.073976) q[2];
rz(1.128528) q[3];
sx q[3];
rz(-2.0566548) q[3];
sx q[3];
rz(1.6245925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45831281) q[0];
sx q[0];
rz(-1.0673609) q[0];
sx q[0];
rz(1.3491032) q[0];
rz(1.5051399) q[1];
sx q[1];
rz(-2.9121193) q[1];
sx q[1];
rz(3.0519003) q[1];
rz(-0.60293052) q[2];
sx q[2];
rz(-1.3984192) q[2];
sx q[2];
rz(2.0655823) q[2];
rz(-2.811089) q[3];
sx q[3];
rz(-1.8310343) q[3];
sx q[3];
rz(-1.9559947) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
