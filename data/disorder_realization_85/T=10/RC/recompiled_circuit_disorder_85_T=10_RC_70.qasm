OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20733362) q[0];
sx q[0];
rz(-2.5512295) q[0];
sx q[0];
rz(-0.37101775) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(-0.59950221) q[1];
sx q[1];
rz(-1.7655656) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8986172) q[0];
sx q[0];
rz(-2.3529422) q[0];
sx q[0];
rz(0.9289766) q[0];
rz(1.0797834) q[2];
sx q[2];
rz(-2.2248189) q[2];
sx q[2];
rz(2.430928) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8810597) q[1];
sx q[1];
rz(-0.36306371) q[1];
sx q[1];
rz(-2.0102324) q[1];
rz(0.81935482) q[3];
sx q[3];
rz(-1.6091533) q[3];
sx q[3];
rz(2.0824144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.084289) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(-2.1526745) q[2];
rz(-2.3890498) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(-0.83077103) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9343524) q[0];
sx q[0];
rz(-0.11226421) q[0];
sx q[0];
rz(1.1799312) q[0];
rz(2.143899) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(2.4172799) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1441919) q[0];
sx q[0];
rz(-1.2574982) q[0];
sx q[0];
rz(2.2749167) q[0];
rz(-1.6547336) q[2];
sx q[2];
rz(-1.3720023) q[2];
sx q[2];
rz(-0.73560152) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1867379) q[1];
sx q[1];
rz(-1.9222944) q[1];
sx q[1];
rz(-0.21912205) q[1];
rz(0.23726666) q[3];
sx q[3];
rz(-1.7627343) q[3];
sx q[3];
rz(2.3290079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2362242) q[2];
sx q[2];
rz(-1.8111818) q[2];
sx q[2];
rz(-0.36188564) q[2];
rz(3.0055255) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(-3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9996465) q[0];
sx q[0];
rz(-1.0936341) q[0];
sx q[0];
rz(-1.3954337) q[0];
rz(2.6793001) q[1];
sx q[1];
rz(-0.42458436) q[1];
sx q[1];
rz(-1.9225072) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8273979) q[0];
sx q[0];
rz(-0.44566804) q[0];
sx q[0];
rz(-0.92339869) q[0];
rz(-pi) q[1];
rz(-0.8692603) q[2];
sx q[2];
rz(-1.7995036) q[2];
sx q[2];
rz(1.9321835) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0137274) q[1];
sx q[1];
rz(-1.8121108) q[1];
sx q[1];
rz(-1.024854) q[1];
rz(-pi) q[2];
rz(1.6270646) q[3];
sx q[3];
rz(-2.8731822) q[3];
sx q[3];
rz(2.066582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7200155) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(-1.4455618) q[2];
rz(2.5727663) q[3];
sx q[3];
rz(-2.2918662) q[3];
sx q[3];
rz(3.0310757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9179984) q[0];
sx q[0];
rz(-2.693394) q[0];
sx q[0];
rz(-0.51825994) q[0];
rz(-0.7154243) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(-0.82675654) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90137705) q[0];
sx q[0];
rz(-0.6745406) q[0];
sx q[0];
rz(-2.4504689) q[0];
rz(-pi) q[1];
rz(-2.3739359) q[2];
sx q[2];
rz(-1.6677688) q[2];
sx q[2];
rz(1.8061639) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7552232) q[1];
sx q[1];
rz(-1.1896903) q[1];
sx q[1];
rz(1.2632881) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72439648) q[3];
sx q[3];
rz(-1.3339692) q[3];
sx q[3];
rz(2.608992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9892019) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(-1.7626804) q[2];
rz(-0.072323024) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(-1.4962083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.82692659) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(-2.0157053) q[0];
rz(-0.90244883) q[1];
sx q[1];
rz(-0.97389644) q[1];
sx q[1];
rz(-2.856423) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42442214) q[0];
sx q[0];
rz(-1.6132014) q[0];
sx q[0];
rz(2.0457343) q[0];
rz(-pi) q[1];
rz(1.9589013) q[2];
sx q[2];
rz(-2.1974807) q[2];
sx q[2];
rz(-0.91458048) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.83708977) q[1];
sx q[1];
rz(-1.2577406) q[1];
sx q[1];
rz(-0.00043991107) q[1];
x q[2];
rz(-0.40163715) q[3];
sx q[3];
rz(-1.65997) q[3];
sx q[3];
rz(-3.0107486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8482762) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(-0.53945333) q[2];
rz(-2.8347677) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(-0.45421281) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25813112) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(-2.9845797) q[0];
rz(0.69333386) q[1];
sx q[1];
rz(-2.2608829) q[1];
sx q[1];
rz(-1.3670115) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6002944) q[0];
sx q[0];
rz(-1.6389264) q[0];
sx q[0];
rz(3.1085204) q[0];
x q[1];
rz(0.579367) q[2];
sx q[2];
rz(-0.89951347) q[2];
sx q[2];
rz(0.52057779) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1937618) q[1];
sx q[1];
rz(-1.8596706) q[1];
sx q[1];
rz(0.08430251) q[1];
rz(-2.7886224) q[3];
sx q[3];
rz(-1.6680696) q[3];
sx q[3];
rz(-0.73247611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29629016) q[2];
sx q[2];
rz(-0.85001105) q[2];
sx q[2];
rz(2.7381251) q[2];
rz(-2.6599595) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(-0.51923716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.24213174) q[0];
sx q[0];
rz(-0.88328981) q[0];
sx q[0];
rz(-2.2677299) q[0];
rz(-2.6938687) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(-1.9708995) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6252977) q[0];
sx q[0];
rz(-1.5646311) q[0];
sx q[0];
rz(-0.023029285) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41202338) q[2];
sx q[2];
rz(-0.12672666) q[2];
sx q[2];
rz(0.7574946) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.16329855) q[1];
sx q[1];
rz(-2.1566609) q[1];
sx q[1];
rz(-2.3818124) q[1];
rz(0.4865173) q[3];
sx q[3];
rz(-2.5518637) q[3];
sx q[3];
rz(0.041134838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0968904) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(-2.5308385) q[2];
rz(-0.47510535) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8996745) q[0];
sx q[0];
rz(-3.0245259) q[0];
sx q[0];
rz(-0.29712594) q[0];
rz(-1.3946474) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(0.64613211) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32983366) q[0];
sx q[0];
rz(-0.23010294) q[0];
sx q[0];
rz(-2.0287201) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4722733) q[2];
sx q[2];
rz(-1.7211282) q[2];
sx q[2];
rz(-2.3538102) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1505193) q[1];
sx q[1];
rz(-1.0725613) q[1];
sx q[1];
rz(-1.9626161) q[1];
rz(-0.93692245) q[3];
sx q[3];
rz(-2.2735032) q[3];
sx q[3];
rz(3.1261409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0769161) q[2];
sx q[2];
rz(-2.1885394) q[2];
sx q[2];
rz(-2.6867552) q[2];
rz(-0.70139766) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(-1.1340244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69944537) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(2.3440857) q[0];
rz(0.51756716) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(3.033175) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.425697) q[0];
sx q[0];
rz(-1.6285537) q[0];
sx q[0];
rz(2.2349368) q[0];
x q[1];
rz(1.7494781) q[2];
sx q[2];
rz(-1.6166501) q[2];
sx q[2];
rz(0.97464857) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3280231) q[1];
sx q[1];
rz(-1.7377825) q[1];
sx q[1];
rz(1.3519253) q[1];
rz(-pi) q[2];
rz(-0.18687825) q[3];
sx q[3];
rz(-2.2664321) q[3];
sx q[3];
rz(-0.38284341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0124399) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(-2.8016395) q[2];
rz(-0.41839504) q[3];
sx q[3];
rz(-0.59643006) q[3];
sx q[3];
rz(0.72559124) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6091992) q[0];
sx q[0];
rz(-2.7476855) q[0];
sx q[0];
rz(2.4627731) q[0];
rz(-0.36418307) q[1];
sx q[1];
rz(-1.697425) q[1];
sx q[1];
rz(-3.0864339) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77627742) q[0];
sx q[0];
rz(-0.3826097) q[0];
sx q[0];
rz(2.9676653) q[0];
rz(-0.55969413) q[2];
sx q[2];
rz(-1.7494546) q[2];
sx q[2];
rz(0.44911227) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.058051) q[1];
sx q[1];
rz(-1.0293048) q[1];
sx q[1];
rz(-0.70152775) q[1];
rz(1.6454562) q[3];
sx q[3];
rz(-2.6844822) q[3];
sx q[3];
rz(-2.8198869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.98383343) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(-2.6514163) q[2];
rz(-3.0040719) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(0.93808758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(0.72538439) q[0];
sx q[0];
rz(-1.890247) q[0];
sx q[0];
rz(2.4631137) q[0];
rz(-0.2086808) q[1];
sx q[1];
rz(-1.122767) q[1];
sx q[1];
rz(-1.541419) q[1];
rz(1.6677042) q[2];
sx q[2];
rz(-1.8816392) q[2];
sx q[2];
rz(0.30602602) q[2];
rz(0.18110885) q[3];
sx q[3];
rz(-0.24387471) q[3];
sx q[3];
rz(-2.7004164) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
