OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3988848) q[0];
sx q[0];
rz(-2.3595915) q[0];
sx q[0];
rz(-1.8703823) q[0];
rz(3.4186163) q[1];
sx q[1];
rz(3.613598) q[1];
sx q[1];
rz(9.4233905) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.815925) q[0];
sx q[0];
rz(-2.7357833) q[0];
sx q[0];
rz(1.9947467) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44714655) q[2];
sx q[2];
rz(-1.6086173) q[2];
sx q[2];
rz(-2.56074) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5056155) q[1];
sx q[1];
rz(-0.50230366) q[1];
sx q[1];
rz(2.99519) q[1];
rz(-1.1997585) q[3];
sx q[3];
rz(-1.8098117) q[3];
sx q[3];
rz(1.0306851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9871621) q[2];
sx q[2];
rz(-0.61750948) q[2];
sx q[2];
rz(2.3922065) q[2];
rz(1.0162214) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(0.40482503) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4352903) q[0];
sx q[0];
rz(-0.82536936) q[0];
sx q[0];
rz(2.170927) q[0];
rz(-1.0372112) q[1];
sx q[1];
rz(-1.7036006) q[1];
sx q[1];
rz(-2.326139) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6204651) q[0];
sx q[0];
rz(-1.6308115) q[0];
sx q[0];
rz(-0.64882664) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30351992) q[2];
sx q[2];
rz(-0.75280658) q[2];
sx q[2];
rz(2.7941861) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7533469) q[1];
sx q[1];
rz(-1.8716295) q[1];
sx q[1];
rz(1.6080329) q[1];
rz(-0.64579441) q[3];
sx q[3];
rz(-1.7495973) q[3];
sx q[3];
rz(-1.7088695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6796391) q[2];
sx q[2];
rz(-1.5660428) q[2];
sx q[2];
rz(2.5088076) q[2];
rz(1.1535545) q[3];
sx q[3];
rz(-2.3735235) q[3];
sx q[3];
rz(-2.8320584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84045029) q[0];
sx q[0];
rz(-1.2121032) q[0];
sx q[0];
rz(0.87483037) q[0];
rz(1.8114999) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(-0.99951807) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32854983) q[0];
sx q[0];
rz(-1.8694436) q[0];
sx q[0];
rz(2.9850328) q[0];
x q[1];
rz(2.3109762) q[2];
sx q[2];
rz(-2.3306371) q[2];
sx q[2];
rz(-2.9142771) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.42400186) q[1];
sx q[1];
rz(-2.058299) q[1];
sx q[1];
rz(-1.8600149) q[1];
rz(-pi) q[2];
rz(0.079654982) q[3];
sx q[3];
rz(-2.0783391) q[3];
sx q[3];
rz(-1.5505276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6040566) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(-2.5615454) q[2];
rz(-2.3245658) q[3];
sx q[3];
rz(-1.3823119) q[3];
sx q[3];
rz(1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76628768) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(2.2312009) q[0];
rz(-0.45122775) q[1];
sx q[1];
rz(-1.5463566) q[1];
sx q[1];
rz(-0.27483637) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9592322) q[0];
sx q[0];
rz(-1.4920456) q[0];
sx q[0];
rz(1.6596646) q[0];
x q[1];
rz(-1.574013) q[2];
sx q[2];
rz(-0.46041691) q[2];
sx q[2];
rz(0.38052961) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.30157629) q[1];
sx q[1];
rz(-0.77743545) q[1];
sx q[1];
rz(-2.2393054) q[1];
rz(-pi) q[2];
rz(0.87423012) q[3];
sx q[3];
rz(-1.4733553) q[3];
sx q[3];
rz(1.3144573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3466907) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(-1.6332731) q[2];
rz(1.1446965) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(-2.9798853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4836924) q[0];
sx q[0];
rz(-1.2150486) q[0];
sx q[0];
rz(-0.99779469) q[0];
rz(0.18355852) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(-1.6246187) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8109587) q[0];
sx q[0];
rz(-1.9348382) q[0];
sx q[0];
rz(-0.90322687) q[0];
x q[1];
rz(-0.5337358) q[2];
sx q[2];
rz(-1.2297451) q[2];
sx q[2];
rz(1.595572) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9285674) q[1];
sx q[1];
rz(-1.9145609) q[1];
sx q[1];
rz(-0.43069559) q[1];
rz(-0.77002854) q[3];
sx q[3];
rz(-0.75776811) q[3];
sx q[3];
rz(-0.023035223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3395485) q[2];
sx q[2];
rz(-0.99207726) q[2];
sx q[2];
rz(0.3240164) q[2];
rz(1.8185395) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(-1.6103305) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76535392) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(-1.8776241) q[0];
rz(2.2309247) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(2.8009159) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9522889) q[0];
sx q[0];
rz(-1.5049107) q[0];
sx q[0];
rz(3.1097263) q[0];
rz(-pi) q[1];
rz(0.31008115) q[2];
sx q[2];
rz(-2.3611464) q[2];
sx q[2];
rz(-2.9030637) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6319879) q[1];
sx q[1];
rz(-1.0972411) q[1];
sx q[1];
rz(1.8297086) q[1];
rz(-0.42231456) q[3];
sx q[3];
rz(-1.1186244) q[3];
sx q[3];
rz(-0.57002588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.95057758) q[2];
sx q[2];
rz(-2.5577929) q[2];
sx q[2];
rz(0.77159709) q[2];
rz(2.5937882) q[3];
sx q[3];
rz(-0.9698202) q[3];
sx q[3];
rz(0.56345338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9724378) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(2.7959438) q[0];
rz(-0.06282839) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(0.46494928) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9879887) q[0];
sx q[0];
rz(-1.1055595) q[0];
sx q[0];
rz(-2.9265755) q[0];
x q[1];
rz(0.99636997) q[2];
sx q[2];
rz(-1.9949706) q[2];
sx q[2];
rz(0.043957274) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9629434) q[1];
sx q[1];
rz(-1.4649676) q[1];
sx q[1];
rz(-0.21957285) q[1];
rz(-2.081359) q[3];
sx q[3];
rz(-3.0887103) q[3];
sx q[3];
rz(-2.3898861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0039625) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(-2.8239992) q[2];
rz(2.5701304) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(0.28731829) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6222318) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(0.28433329) q[0];
rz(-2.590086) q[1];
sx q[1];
rz(-2.9998144) q[1];
sx q[1];
rz(3.0632339) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0067622234) q[0];
sx q[0];
rz(-1.5015331) q[0];
sx q[0];
rz(2.3368821) q[0];
rz(-pi) q[1];
rz(-1.8163082) q[2];
sx q[2];
rz(-1.3405372) q[2];
sx q[2];
rz(-2.2895209) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9662712) q[1];
sx q[1];
rz(-2.2178855) q[1];
sx q[1];
rz(-2.9233169) q[1];
x q[2];
rz(-1.3571635) q[3];
sx q[3];
rz(-2.9812818) q[3];
sx q[3];
rz(-2.8182639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4006965) q[2];
sx q[2];
rz(-2.233278) q[2];
sx q[2];
rz(-0.25137869) q[2];
rz(-2.5583983) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(-1.73197) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79779977) q[0];
sx q[0];
rz(-3.0594337) q[0];
sx q[0];
rz(3.0902241) q[0];
rz(-0.92357606) q[1];
sx q[1];
rz(-0.66134614) q[1];
sx q[1];
rz(-2.267568) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066814518) q[0];
sx q[0];
rz(-1.5794808) q[0];
sx q[0];
rz(-2.3735223) q[0];
rz(-pi) q[1];
rz(2.2888695) q[2];
sx q[2];
rz(-2.6211779) q[2];
sx q[2];
rz(-2.4783217) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6898432) q[1];
sx q[1];
rz(-0.50536957) q[1];
sx q[1];
rz(1.4228574) q[1];
rz(-pi) q[2];
rz(-2.6690528) q[3];
sx q[3];
rz(-0.66415411) q[3];
sx q[3];
rz(3.1411375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.727227) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(0.61974636) q[2];
rz(-1.184458) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(-1.7782036) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0062362) q[0];
sx q[0];
rz(-1.0422491) q[0];
sx q[0];
rz(-2.4172879) q[0];
rz(2.9528217) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(1.1788517) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2154685) q[0];
sx q[0];
rz(-1.3660396) q[0];
sx q[0];
rz(1.8490851) q[0];
x q[1];
rz(-3.0145698) q[2];
sx q[2];
rz(-1.4782895) q[2];
sx q[2];
rz(-0.18735838) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5621592) q[1];
sx q[1];
rz(-2.4907618) q[1];
sx q[1];
rz(1.8821554) q[1];
x q[2];
rz(-0.3801109) q[3];
sx q[3];
rz(-2.5862525) q[3];
sx q[3];
rz(-1.7172608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0835691) q[2];
sx q[2];
rz(-2.100779) q[2];
sx q[2];
rz(-2.2422092) q[2];
rz(2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(-1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5205004) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(-2.3836366) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(1.6451251) q[2];
sx q[2];
rz(-2.7004514) q[2];
sx q[2];
rz(-2.2218291) q[2];
rz(0.27579565) q[3];
sx q[3];
rz(-1.4281359) q[3];
sx q[3];
rz(1.7208163) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];