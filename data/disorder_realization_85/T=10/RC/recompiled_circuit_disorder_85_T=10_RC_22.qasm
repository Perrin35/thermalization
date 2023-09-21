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
rz(1.3289514) q[0];
sx q[0];
rz(-2.0094123) q[0];
sx q[0];
rz(0.89226512) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0618093) q[2];
sx q[2];
rz(-0.91677374) q[2];
sx q[2];
rz(0.71066463) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7264759) q[1];
sx q[1];
rz(-1.8980025) q[1];
sx q[1];
rz(2.9813558) q[1];
rz(-0.81935482) q[3];
sx q[3];
rz(-1.5324394) q[3];
sx q[3];
rz(2.0824144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0573037) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(-0.98891813) q[2];
rz(0.75254285) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(2.3108216) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9343524) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(1.9616615) q[0];
rz(2.143899) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(-0.72431272) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82942078) q[0];
sx q[0];
rz(-2.2342626) q[0];
sx q[0];
rz(-0.40191606) q[0];
rz(-pi) q[1];
x q[1];
rz(1.486859) q[2];
sx q[2];
rz(-1.7695904) q[2];
sx q[2];
rz(-2.4059911) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.449125) q[1];
sx q[1];
rz(-1.3652703) q[1];
sx q[1];
rz(-1.9301901) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3734666) q[3];
sx q[3];
rz(-1.337968) q[3];
sx q[3];
rz(-0.80430921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.90536845) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(-0.36188564) q[2];
rz(-3.0055255) q[3];
sx q[3];
rz(-2.5858904) q[3];
sx q[3];
rz(-3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1419462) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(-1.746159) q[0];
rz(0.46229258) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(1.2190855) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31419471) q[0];
sx q[0];
rz(-2.6959246) q[0];
sx q[0];
rz(2.218194) q[0];
rz(2.2723324) q[2];
sx q[2];
rz(-1.7995036) q[2];
sx q[2];
rz(1.9321835) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29875007) q[1];
sx q[1];
rz(-2.0992273) q[1];
sx q[1];
rz(-0.28038402) q[1];
rz(-pi) q[2];
rz(-1.3027906) q[3];
sx q[3];
rz(-1.5558814) q[3];
sx q[3];
rz(-0.55004317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7200155) q[2];
sx q[2];
rz(-2.4013459) q[2];
sx q[2];
rz(1.6960309) q[2];
rz(0.56882632) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(3.0310757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22359426) q[0];
sx q[0];
rz(-2.693394) q[0];
sx q[0];
rz(-2.6233327) q[0];
rz(2.4261684) q[1];
sx q[1];
rz(-1.1162858) q[1];
sx q[1];
rz(-2.3148361) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90137705) q[0];
sx q[0];
rz(-0.6745406) q[0];
sx q[0];
rz(-2.4504689) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76765676) q[2];
sx q[2];
rz(-1.6677688) q[2];
sx q[2];
rz(1.3354288) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0485059) q[1];
sx q[1];
rz(-2.6566681) q[1];
sx q[1];
rz(-2.4946458) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8825674) q[3];
sx q[3];
rz(-2.2707553) q[3];
sx q[3];
rz(1.242897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15239079) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(-1.7626804) q[2];
rz(-0.072323024) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(-1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82692659) q[0];
sx q[0];
rz(-0.62012726) q[0];
sx q[0];
rz(-2.0157053) q[0];
rz(-0.90244883) q[1];
sx q[1];
rz(-0.97389644) q[1];
sx q[1];
rz(0.28516969) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0774768) q[0];
sx q[0];
rz(-0.47668326) q[0];
sx q[0];
rz(1.478273) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1826913) q[2];
sx q[2];
rz(-2.1974807) q[2];
sx q[2];
rz(0.91458048) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.8385182) q[1];
sx q[1];
rz(-2.8285366) q[1];
sx q[1];
rz(1.5694373) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4739591) q[3];
sx q[3];
rz(-1.9707465) q[3];
sx q[3];
rz(1.6638343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8482762) q[2];
sx q[2];
rz(-0.50555503) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25813112) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(-2.9845797) q[0];
rz(-2.4482588) q[1];
sx q[1];
rz(-2.2608829) q[1];
sx q[1];
rz(-1.3670115) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0317504) q[0];
sx q[0];
rz(-1.5378008) q[0];
sx q[0];
rz(-1.6389636) q[0];
x q[1];
rz(0.96732803) q[2];
sx q[2];
rz(-0.8562932) q[2];
sx q[2];
rz(2.8514903) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.65946603) q[1];
sx q[1];
rz(-2.8409993) q[1];
sx q[1];
rz(-1.8468922) q[1];
x q[2];
rz(-0.27512392) q[3];
sx q[3];
rz(-2.776006) q[3];
sx q[3];
rz(2.5610353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.29629016) q[2];
sx q[2];
rz(-0.85001105) q[2];
sx q[2];
rz(2.7381251) q[2];
rz(-0.48163313) q[3];
sx q[3];
rz(-2.0694331) q[3];
sx q[3];
rz(2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8994609) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(-2.2677299) q[0];
rz(-0.44772398) q[1];
sx q[1];
rz(-0.73900765) q[1];
sx q[1];
rz(1.1706932) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51629492) q[0];
sx q[0];
rz(-1.5646311) q[0];
sx q[0];
rz(3.1185634) q[0];
rz(-pi) q[1];
rz(2.7295693) q[2];
sx q[2];
rz(-3.014866) q[2];
sx q[2];
rz(0.7574946) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2504187) q[1];
sx q[1];
rz(-0.95953566) q[1];
sx q[1];
rz(-2.3120018) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4865173) q[3];
sx q[3];
rz(-0.58972893) q[3];
sx q[3];
rz(-0.041134838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0447023) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(0.61075413) q[2];
rz(-2.6664873) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8996745) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(-2.8444667) q[0];
rz(-1.3946474) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(-2.4954605) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6883811) q[0];
sx q[0];
rz(-1.4697945) q[0];
sx q[0];
rz(1.7779011) q[0];
x q[1];
rz(2.8201639) q[2];
sx q[2];
rz(-2.6476963) q[2];
sx q[2];
rz(1.0682046) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.38478002) q[1];
sx q[1];
rz(-1.2287178) q[1];
sx q[1];
rz(-0.53201075) q[1];
rz(-pi) q[2];
rz(0.81031462) q[3];
sx q[3];
rz(-2.0397566) q[3];
sx q[3];
rz(-1.9988434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0769161) q[2];
sx q[2];
rz(-2.1885394) q[2];
sx q[2];
rz(-2.6867552) q[2];
rz(-2.440195) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(1.1340244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.69944537) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(-2.3440857) q[0];
rz(0.51756716) q[1];
sx q[1];
rz(-2.3289754) q[1];
sx q[1];
rz(0.10841766) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3603044) q[0];
sx q[0];
rz(-0.66626781) q[0];
sx q[0];
rz(1.4772619) q[0];
rz(-pi) q[1];
rz(1.7494781) q[2];
sx q[2];
rz(-1.5249426) q[2];
sx q[2];
rz(2.1669441) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3280231) q[1];
sx q[1];
rz(-1.4038101) q[1];
sx q[1];
rz(-1.7896673) q[1];
rz(0.18687825) q[3];
sx q[3];
rz(-0.87516057) q[3];
sx q[3];
rz(2.7587492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0124399) q[2];
sx q[2];
rz(-1.3468578) q[2];
sx q[2];
rz(-2.8016395) q[2];
rz(-0.41839504) q[3];
sx q[3];
rz(-0.59643006) q[3];
sx q[3];
rz(-2.4160014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5323935) q[0];
sx q[0];
rz(-2.7476855) q[0];
sx q[0];
rz(-0.6788196) q[0];
rz(0.36418307) q[1];
sx q[1];
rz(-1.697425) q[1];
sx q[1];
rz(-0.055158786) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3653152) q[0];
sx q[0];
rz(-2.758983) q[0];
sx q[0];
rz(-0.17392735) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55969413) q[2];
sx q[2];
rz(-1.7494546) q[2];
sx q[2];
rz(-0.44911227) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0810869) q[1];
sx q[1];
rz(-0.85716893) q[1];
sx q[1];
rz(-2.39141) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1147898) q[3];
sx q[3];
rz(-1.5378693) q[3];
sx q[3];
rz(1.1820716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98383343) q[2];
sx q[2];
rz(-2.1201717) q[2];
sx q[2];
rz(0.49017635) q[2];
rz(-0.13752078) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(-0.93808758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72538439) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(0.2086808) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(1.4738884) q[2];
sx q[2];
rz(-1.2599535) q[2];
sx q[2];
rz(-2.8355666) q[2];
rz(1.5260074) q[3];
sx q[3];
rz(-1.3309892) q[3];
sx q[3];
rz(0.25467024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
