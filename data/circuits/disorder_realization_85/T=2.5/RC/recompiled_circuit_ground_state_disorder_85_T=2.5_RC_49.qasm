OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6326132) q[0];
sx q[0];
rz(-2.296663) q[0];
sx q[0];
rz(2.4281969) q[0];
rz(2.9333935) q[1];
sx q[1];
rz(2.9543076) q[1];
sx q[1];
rz(2.0050144) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5644585) q[0];
sx q[0];
rz(-2.4043879) q[0];
sx q[0];
rz(-1.2080824) q[0];
rz(2.2520653) q[2];
sx q[2];
rz(-1.5233381) q[2];
sx q[2];
rz(-2.0840381) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.40219992) q[1];
sx q[1];
rz(-0.2570411) q[1];
sx q[1];
rz(0.57769026) q[1];
rz(-pi) q[2];
rz(2.7375701) q[3];
sx q[3];
rz(-2.8131313) q[3];
sx q[3];
rz(-0.31620178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.95181695) q[2];
sx q[2];
rz(-1.8686029) q[2];
sx q[2];
rz(0.58958685) q[2];
rz(3.12674) q[3];
sx q[3];
rz(-1.868052) q[3];
sx q[3];
rz(2.5428037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0066765229) q[0];
sx q[0];
rz(-2.0966457) q[0];
sx q[0];
rz(2.6892804) q[0];
rz(0.85330883) q[1];
sx q[1];
rz(-0.49214688) q[1];
sx q[1];
rz(2.7795627) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25681296) q[0];
sx q[0];
rz(-2.1636336) q[0];
sx q[0];
rz(-1.9624016) q[0];
rz(-pi) q[1];
rz(0.30899991) q[2];
sx q[2];
rz(-0.51181343) q[2];
sx q[2];
rz(-1.9644033) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7343062) q[1];
sx q[1];
rz(-0.9883259) q[1];
sx q[1];
rz(-0.72093074) q[1];
rz(0.70254247) q[3];
sx q[3];
rz(-1.4558534) q[3];
sx q[3];
rz(0.045199423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8028458) q[2];
sx q[2];
rz(-0.57463988) q[2];
sx q[2];
rz(2.23526) q[2];
rz(1.686442) q[3];
sx q[3];
rz(-0.95821277) q[3];
sx q[3];
rz(-1.5549392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(2.0543268) q[0];
sx q[0];
rz(-1.9540906) q[0];
sx q[0];
rz(-0.97386709) q[0];
rz(0.57526881) q[1];
sx q[1];
rz(-1.560805) q[1];
sx q[1];
rz(1.5781933) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3546936) q[0];
sx q[0];
rz(-1.5074566) q[0];
sx q[0];
rz(3.1194434) q[0];
rz(-0.84196426) q[2];
sx q[2];
rz(-1.2281111) q[2];
sx q[2];
rz(2.2322643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.089019589) q[1];
sx q[1];
rz(-2.0327489) q[1];
sx q[1];
rz(-2.3506052) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3493912) q[3];
sx q[3];
rz(-0.81238895) q[3];
sx q[3];
rz(-2.098907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.32597184) q[2];
sx q[2];
rz(-1.9591363) q[2];
sx q[2];
rz(1.5284485) q[2];
rz(3.1260887) q[3];
sx q[3];
rz(-1.1752693) q[3];
sx q[3];
rz(-1.4874682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4267321) q[0];
sx q[0];
rz(-2.4437014) q[0];
sx q[0];
rz(-3.0791855) q[0];
rz(-1.9852253) q[1];
sx q[1];
rz(-0.9461177) q[1];
sx q[1];
rz(2.3074493) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7251563) q[0];
sx q[0];
rz(-1.8653989) q[0];
sx q[0];
rz(-1.0480773) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6970728) q[2];
sx q[2];
rz(-1.6202721) q[2];
sx q[2];
rz(0.49903185) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3662009) q[1];
sx q[1];
rz(-1.5171261) q[1];
sx q[1];
rz(0.97241511) q[1];
x q[2];
rz(1.2150202) q[3];
sx q[3];
rz(-2.5952383) q[3];
sx q[3];
rz(-0.90701963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9741029) q[2];
sx q[2];
rz(-0.17593273) q[2];
sx q[2];
rz(1.9195732) q[2];
rz(2.7410298) q[3];
sx q[3];
rz(-2.1192854) q[3];
sx q[3];
rz(-2.9266973) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062531384) q[0];
sx q[0];
rz(-2.489594) q[0];
sx q[0];
rz(0.25667152) q[0];
rz(-1.3447064) q[1];
sx q[1];
rz(-1.2213629) q[1];
sx q[1];
rz(1.7324956) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2715816) q[0];
sx q[0];
rz(-2.3973532) q[0];
sx q[0];
rz(1.1209473) q[0];
x q[1];
rz(-1.8518894) q[2];
sx q[2];
rz(-1.7651193) q[2];
sx q[2];
rz(1.3240726) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0540648) q[1];
sx q[1];
rz(-1.9913186) q[1];
sx q[1];
rz(1.6144714) q[1];
x q[2];
rz(-0.37568521) q[3];
sx q[3];
rz(-0.70313625) q[3];
sx q[3];
rz(0.60763728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1008489) q[2];
sx q[2];
rz(-0.70793968) q[2];
sx q[2];
rz(-2.9948998) q[2];
rz(-1.3692859) q[3];
sx q[3];
rz(-0.36543235) q[3];
sx q[3];
rz(-2.9036314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.9323102) q[0];
sx q[0];
rz(-0.76618659) q[0];
sx q[0];
rz(0.76195088) q[0];
rz(0.85488287) q[1];
sx q[1];
rz(-1.8652752) q[1];
sx q[1];
rz(-0.68861047) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013278339) q[0];
sx q[0];
rz(-1.4869031) q[0];
sx q[0];
rz(1.1021139) q[0];
x q[1];
rz(0.18085955) q[2];
sx q[2];
rz(-0.11520371) q[2];
sx q[2];
rz(1.9242788) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.84160596) q[1];
sx q[1];
rz(-3.0191023) q[1];
sx q[1];
rz(0.56614205) q[1];
x q[2];
rz(1.1375539) q[3];
sx q[3];
rz(-1.8747624) q[3];
sx q[3];
rz(2.7734625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.40053549) q[2];
sx q[2];
rz(-2.1964938) q[2];
sx q[2];
rz(2.1607024) q[2];
rz(2.8716221) q[3];
sx q[3];
rz(-1.2308086) q[3];
sx q[3];
rz(0.009416906) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9735759) q[0];
sx q[0];
rz(-1.6930169) q[0];
sx q[0];
rz(-2.5981405) q[0];
rz(1.5501529) q[1];
sx q[1];
rz(-0.88528577) q[1];
sx q[1];
rz(-2.1563704) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60632432) q[0];
sx q[0];
rz(-2.599975) q[0];
sx q[0];
rz(0.4569502) q[0];
rz(1.0205054) q[2];
sx q[2];
rz(-2.5021363) q[2];
sx q[2];
rz(0.24580641) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4812659) q[1];
sx q[1];
rz(-0.72176343) q[1];
sx q[1];
rz(-0.97246171) q[1];
x q[2];
rz(-0.51237088) q[3];
sx q[3];
rz(-1.2768043) q[3];
sx q[3];
rz(1.9749354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8103509) q[2];
sx q[2];
rz(-2.6369429) q[2];
sx q[2];
rz(1.8028367) q[2];
rz(1.6797558) q[3];
sx q[3];
rz(-1.5509408) q[3];
sx q[3];
rz(2.4707879) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2115693) q[0];
sx q[0];
rz(-1.1299364) q[0];
sx q[0];
rz(1.9071963) q[0];
rz(-1.7878388) q[1];
sx q[1];
rz(-2.6852971) q[1];
sx q[1];
rz(-0.096045883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9864842) q[0];
sx q[0];
rz(-0.86244118) q[0];
sx q[0];
rz(0.95671311) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4951316) q[2];
sx q[2];
rz(-2.1096418) q[2];
sx q[2];
rz(-2.7494753) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2224642) q[1];
sx q[1];
rz(-1.9507196) q[1];
sx q[1];
rz(-1.6187731) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3013294) q[3];
sx q[3];
rz(-1.5397159) q[3];
sx q[3];
rz(-2.4781991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8058406) q[2];
sx q[2];
rz(-0.38613191) q[2];
sx q[2];
rz(1.1248355) q[2];
rz(-2.2624894) q[3];
sx q[3];
rz(-1.9156888) q[3];
sx q[3];
rz(0.44900289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6833471) q[0];
sx q[0];
rz(-1.6793716) q[0];
sx q[0];
rz(2.4264088) q[0];
rz(2.8257651) q[1];
sx q[1];
rz(-0.68604699) q[1];
sx q[1];
rz(-2.979523) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93955646) q[0];
sx q[0];
rz(-1.1456155) q[0];
sx q[0];
rz(-0.74799882) q[0];
rz(-0.58312441) q[2];
sx q[2];
rz(-1.4004502) q[2];
sx q[2];
rz(-1.9161461) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9980364) q[1];
sx q[1];
rz(-0.85696917) q[1];
sx q[1];
rz(-1.803323) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81406419) q[3];
sx q[3];
rz(-0.45032516) q[3];
sx q[3];
rz(-1.7801156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.38966236) q[2];
sx q[2];
rz(-1.1952362) q[2];
sx q[2];
rz(1.8771089) q[2];
rz(-1.8090931) q[3];
sx q[3];
rz(-1.2800565) q[3];
sx q[3];
rz(0.31681028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7952809) q[0];
sx q[0];
rz(-0.77487159) q[0];
sx q[0];
rz(-2.0922022) q[0];
rz(0.28114444) q[1];
sx q[1];
rz(-1.0422336) q[1];
sx q[1];
rz(-1.8195456) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3157177) q[0];
sx q[0];
rz(-2.115504) q[0];
sx q[0];
rz(-1.8648575) q[0];
x q[1];
rz(-2.5192267) q[2];
sx q[2];
rz(-2.4020148) q[2];
sx q[2];
rz(-3.0837631) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7483298) q[1];
sx q[1];
rz(-2.7460881) q[1];
sx q[1];
rz(-1.9935745) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38615245) q[3];
sx q[3];
rz(-1.6623467) q[3];
sx q[3];
rz(1.9492675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.37810024) q[2];
sx q[2];
rz(-1.3214448) q[2];
sx q[2];
rz(-3.0901001) q[2];
rz(1.817305) q[3];
sx q[3];
rz(-1.6591502) q[3];
sx q[3];
rz(1.9342559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1235724) q[0];
sx q[0];
rz(-1.3326895) q[0];
sx q[0];
rz(1.4154758) q[0];
rz(0.82072683) q[1];
sx q[1];
rz(-1.4438933) q[1];
sx q[1];
rz(-1.874598) q[1];
rz(0.61931284) q[2];
sx q[2];
rz(-1.3765341) q[2];
sx q[2];
rz(1.7159009) q[2];
rz(-0.5140082) q[3];
sx q[3];
rz(-1.8655984) q[3];
sx q[3];
rz(1.8023103) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
