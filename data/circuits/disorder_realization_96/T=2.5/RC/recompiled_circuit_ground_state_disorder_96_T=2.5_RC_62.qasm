OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.628217) q[0];
sx q[0];
rz(-1.4613419) q[0];
sx q[0];
rz(-2.2267377) q[0];
rz(2.272361) q[1];
sx q[1];
rz(-0.72328049) q[1];
sx q[1];
rz(-1.7887315) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3769414) q[0];
sx q[0];
rz(-2.1989555) q[0];
sx q[0];
rz(0.92002042) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3981744) q[2];
sx q[2];
rz(-1.2890491) q[2];
sx q[2];
rz(0.27458686) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9761524) q[1];
sx q[1];
rz(-0.70164499) q[1];
sx q[1];
rz(1.8681504) q[1];
rz(-pi) q[2];
rz(-2.9269591) q[3];
sx q[3];
rz(-2.1691536) q[3];
sx q[3];
rz(-0.90902381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73504084) q[2];
sx q[2];
rz(-0.81321365) q[2];
sx q[2];
rz(-1.0068033) q[2];
rz(1.6793647) q[3];
sx q[3];
rz(-1.4454602) q[3];
sx q[3];
rz(1.6026976) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8524858) q[0];
sx q[0];
rz(-2.2302581) q[0];
sx q[0];
rz(-0.12595969) q[0];
rz(-1.1360629) q[1];
sx q[1];
rz(-1.845153) q[1];
sx q[1];
rz(-1.0596641) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9618767) q[0];
sx q[0];
rz(-1.5715412) q[0];
sx q[0];
rz(-0.0019093328) q[0];
rz(-pi) q[1];
rz(-0.95086348) q[2];
sx q[2];
rz(-2.1427832) q[2];
sx q[2];
rz(0.36617719) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6525065) q[1];
sx q[1];
rz(-2.2678657) q[1];
sx q[1];
rz(0.98784312) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7823388) q[3];
sx q[3];
rz(-1.5240961) q[3];
sx q[3];
rz(-0.31622313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58553592) q[2];
sx q[2];
rz(-1.5289331) q[2];
sx q[2];
rz(2.5118828) q[2];
rz(-1.9017879) q[3];
sx q[3];
rz(-2.3972378) q[3];
sx q[3];
rz(2.0201717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.020551) q[0];
sx q[0];
rz(-2.910055) q[0];
sx q[0];
rz(1.2303906) q[0];
rz(-0.70368272) q[1];
sx q[1];
rz(-2.2129702) q[1];
sx q[1];
rz(-1.5603265) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9890922) q[0];
sx q[0];
rz(-1.9643503) q[0];
sx q[0];
rz(-0.016274289) q[0];
x q[1];
rz(-2.4034924) q[2];
sx q[2];
rz(-1.4261275) q[2];
sx q[2];
rz(-1.5371145) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.091895122) q[1];
sx q[1];
rz(-1.1593684) q[1];
sx q[1];
rz(2.0446214) q[1];
rz(-pi) q[2];
rz(2.3558667) q[3];
sx q[3];
rz(-1.8453987) q[3];
sx q[3];
rz(2.6089607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8065717) q[2];
sx q[2];
rz(-2.0609914) q[2];
sx q[2];
rz(3.0873155) q[2];
rz(-1.0862167) q[3];
sx q[3];
rz(-0.79038668) q[3];
sx q[3];
rz(0.49200341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7739968) q[0];
sx q[0];
rz(-0.01570276) q[0];
sx q[0];
rz(2.1570461) q[0];
rz(1.4060075) q[1];
sx q[1];
rz(-1.6008629) q[1];
sx q[1];
rz(-0.23449177) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0011976) q[0];
sx q[0];
rz(-2.2463628) q[0];
sx q[0];
rz(0.35517774) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86799829) q[2];
sx q[2];
rz(-2.7722087) q[2];
sx q[2];
rz(-2.546865) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2114276) q[1];
sx q[1];
rz(-2.549173) q[1];
sx q[1];
rz(0.38229015) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.111895) q[3];
sx q[3];
rz(-2.1468024) q[3];
sx q[3];
rz(2.821049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24388127) q[2];
sx q[2];
rz(-0.66456866) q[2];
sx q[2];
rz(-0.74771869) q[2];
rz(-2.054935) q[3];
sx q[3];
rz(-2.1472411) q[3];
sx q[3];
rz(-2.7896816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20027593) q[0];
sx q[0];
rz(-2.2714697) q[0];
sx q[0];
rz(-1.548832) q[0];
rz(3.0039655) q[1];
sx q[1];
rz(-2.177114) q[1];
sx q[1];
rz(0.67759883) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1573059) q[0];
sx q[0];
rz(-1.56347) q[0];
sx q[0];
rz(0.013989438) q[0];
rz(-pi) q[1];
rz(-0.19385152) q[2];
sx q[2];
rz(-1.7192013) q[2];
sx q[2];
rz(-2.1231164) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4921129) q[1];
sx q[1];
rz(-1.4339916) q[1];
sx q[1];
rz(-0.0059195267) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4047071) q[3];
sx q[3];
rz(-1.3892655) q[3];
sx q[3];
rz(-1.5893861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4981093) q[2];
sx q[2];
rz(-0.52916873) q[2];
sx q[2];
rz(-0.37437487) q[2];
rz(0.74786413) q[3];
sx q[3];
rz(-1.6467843) q[3];
sx q[3];
rz(-1.098746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13136524) q[0];
sx q[0];
rz(-1.7759198) q[0];
sx q[0];
rz(-0.22450547) q[0];
rz(0.35047105) q[1];
sx q[1];
rz(-1.9572565) q[1];
sx q[1];
rz(1.9559466) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3315546) q[0];
sx q[0];
rz(-1.4228978) q[0];
sx q[0];
rz(0.055813544) q[0];
rz(-pi) q[1];
rz(-2.6780532) q[2];
sx q[2];
rz(-1.4847418) q[2];
sx q[2];
rz(1.8870837) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9486432) q[1];
sx q[1];
rz(-1.1426326) q[1];
sx q[1];
rz(2.3938426) q[1];
x q[2];
rz(-0.91573651) q[3];
sx q[3];
rz(-1.011542) q[3];
sx q[3];
rz(-0.19266927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6422358) q[2];
sx q[2];
rz(-0.88358742) q[2];
sx q[2];
rz(0.16481608) q[2];
rz(-0.76588255) q[3];
sx q[3];
rz(-2.3131148) q[3];
sx q[3];
rz(0.59144056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0931452) q[0];
sx q[0];
rz(-3.1146545) q[0];
sx q[0];
rz(-0.41937605) q[0];
rz(-1.5834437) q[1];
sx q[1];
rz(-0.42834586) q[1];
sx q[1];
rz(1.8289808) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87154138) q[0];
sx q[0];
rz(-1.7333374) q[0];
sx q[0];
rz(0.97515653) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75757005) q[2];
sx q[2];
rz(-0.84572863) q[2];
sx q[2];
rz(-2.6427694) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6083642) q[1];
sx q[1];
rz(-1.8177515) q[1];
sx q[1];
rz(-2.9465531) q[1];
rz(-pi) q[2];
rz(2.873654) q[3];
sx q[3];
rz(-1.6959785) q[3];
sx q[3];
rz(2.1942558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7340362) q[2];
sx q[2];
rz(-0.67136705) q[2];
sx q[2];
rz(0.35801312) q[2];
rz(-0.19896209) q[3];
sx q[3];
rz(-2.3301061) q[3];
sx q[3];
rz(3.0437886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1426549) q[0];
sx q[0];
rz(-2.1079347) q[0];
sx q[0];
rz(0.79310966) q[0];
rz(-2.2456887) q[1];
sx q[1];
rz(-1.9205807) q[1];
sx q[1];
rz(-2.6633247) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9288328) q[0];
sx q[0];
rz(-0.46657714) q[0];
sx q[0];
rz(-0.27884941) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4594166) q[2];
sx q[2];
rz(-2.6016005) q[2];
sx q[2];
rz(-3.1201975) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7859808) q[1];
sx q[1];
rz(-2.0352239) q[1];
sx q[1];
rz(0.29504882) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0317814) q[3];
sx q[3];
rz(-1.0581019) q[3];
sx q[3];
rz(-1.5231665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.521495) q[2];
sx q[2];
rz(-1.564881) q[2];
sx q[2];
rz(-0.0077956789) q[2];
rz(-1.9589641) q[3];
sx q[3];
rz(-1.3820442) q[3];
sx q[3];
rz(1.8462605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9919392) q[0];
sx q[0];
rz(-2.6872334) q[0];
sx q[0];
rz(0.40865189) q[0];
rz(2.0463792) q[1];
sx q[1];
rz(-1.6588255) q[1];
sx q[1];
rz(-0.85831395) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8105111) q[0];
sx q[0];
rz(-2.9797922) q[0];
sx q[0];
rz(-1.9952379) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0235262) q[2];
sx q[2];
rz(-1.451734) q[2];
sx q[2];
rz(-0.43242142) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5206362) q[1];
sx q[1];
rz(-2.9163216) q[1];
sx q[1];
rz(-1.0566637) q[1];
x q[2];
rz(0.32118843) q[3];
sx q[3];
rz(-0.62231718) q[3];
sx q[3];
rz(0.27143196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3506763) q[2];
sx q[2];
rz(-0.75296593) q[2];
sx q[2];
rz(-2.47993) q[2];
rz(2.3412797) q[3];
sx q[3];
rz(-1.9789663) q[3];
sx q[3];
rz(-2.3683828) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3043587) q[0];
sx q[0];
rz(-0.14166129) q[0];
sx q[0];
rz(-0.71453553) q[0];
rz(-2.6658658) q[1];
sx q[1];
rz(-0.58554244) q[1];
sx q[1];
rz(0.48019662) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4775608) q[0];
sx q[0];
rz(-0.98958221) q[0];
sx q[0];
rz(-2.1607735) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91861208) q[2];
sx q[2];
rz(-1.3038074) q[2];
sx q[2];
rz(-0.3912386) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.811753) q[1];
sx q[1];
rz(-1.2851213) q[1];
sx q[1];
rz(-0.67606725) q[1];
x q[2];
rz(2.8586724) q[3];
sx q[3];
rz(-1.1503387) q[3];
sx q[3];
rz(2.450301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.51018888) q[2];
sx q[2];
rz(-1.9437342) q[2];
sx q[2];
rz(-0.27604827) q[2];
rz(2.7035642) q[3];
sx q[3];
rz(-1.6249526) q[3];
sx q[3];
rz(0.62694222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.940687) q[0];
sx q[0];
rz(-1.9872266) q[0];
sx q[0];
rz(-3.0845597) q[0];
rz(-3.1055462) q[1];
sx q[1];
rz(-1.5279952) q[1];
sx q[1];
rz(-1.5855018) q[1];
rz(2.0519999) q[2];
sx q[2];
rz(-0.54554987) q[2];
sx q[2];
rz(0.050532374) q[2];
rz(-0.25419079) q[3];
sx q[3];
rz(-1.4387476) q[3];
sx q[3];
rz(0.1968255) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
