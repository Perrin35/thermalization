OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4560661) q[0];
sx q[0];
rz(-0.38903061) q[0];
sx q[0];
rz(2.2580137) q[0];
rz(3.1318624) q[1];
sx q[1];
rz(4.598773) q[1];
sx q[1];
rz(7.4814319) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0177901) q[0];
sx q[0];
rz(-1.4746116) q[0];
sx q[0];
rz(-2.9642446) q[0];
x q[1];
rz(-0.3281524) q[2];
sx q[2];
rz(-0.80291623) q[2];
sx q[2];
rz(0.83628718) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6557505) q[1];
sx q[1];
rz(-1.9006923) q[1];
sx q[1];
rz(-0.54539036) q[1];
x q[2];
rz(-1.8331068) q[3];
sx q[3];
rz(-2.9648818) q[3];
sx q[3];
rz(1.7929329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1951695) q[2];
sx q[2];
rz(-2.158458) q[2];
sx q[2];
rz(-0.18134376) q[2];
rz(-2.8803853) q[3];
sx q[3];
rz(-1.8758592) q[3];
sx q[3];
rz(-0.75631022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8392035) q[0];
sx q[0];
rz(-0.2897245) q[0];
sx q[0];
rz(2.7547577) q[0];
rz(0.50239262) q[1];
sx q[1];
rz(-2.1680809) q[1];
sx q[1];
rz(1.5418672) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72774784) q[0];
sx q[0];
rz(-1.1607329) q[0];
sx q[0];
rz(-0.64322612) q[0];
rz(-pi) q[1];
rz(-2.9773211) q[2];
sx q[2];
rz(-1.1195682) q[2];
sx q[2];
rz(2.3300366) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9215556) q[1];
sx q[1];
rz(-1.5350635) q[1];
sx q[1];
rz(1.5807371) q[1];
rz(2.6582791) q[3];
sx q[3];
rz(-1.1477074) q[3];
sx q[3];
rz(-1.9139293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.60275045) q[2];
sx q[2];
rz(-2.2637612) q[2];
sx q[2];
rz(1.4146457) q[2];
rz(0.85033068) q[3];
sx q[3];
rz(-0.43262216) q[3];
sx q[3];
rz(1.6833646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26281115) q[0];
sx q[0];
rz(-1.5314064) q[0];
sx q[0];
rz(-0.61022726) q[0];
rz(-1.2894851) q[1];
sx q[1];
rz(-2.162343) q[1];
sx q[1];
rz(2.1496225) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1410071) q[0];
sx q[0];
rz(-1.3864543) q[0];
sx q[0];
rz(-1.0087476) q[0];
rz(-pi) q[1];
rz(0.66833468) q[2];
sx q[2];
rz(-1.4112345) q[2];
sx q[2];
rz(-2.4192686) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.386258) q[1];
sx q[1];
rz(-1.4812246) q[1];
sx q[1];
rz(0.44134015) q[1];
rz(-pi) q[2];
rz(1.3258341) q[3];
sx q[3];
rz(-2.0162597) q[3];
sx q[3];
rz(2.9374292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.758574) q[2];
sx q[2];
rz(-2.980361) q[2];
sx q[2];
rz(-2.8733011) q[2];
rz(0.39408436) q[3];
sx q[3];
rz(-1.9106617) q[3];
sx q[3];
rz(-2.9530853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2878993) q[0];
sx q[0];
rz(-2.6278966) q[0];
sx q[0];
rz(-0.40507856) q[0];
rz(0.69008094) q[1];
sx q[1];
rz(-1.9837374) q[1];
sx q[1];
rz(-2.4437723) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47397428) q[0];
sx q[0];
rz(-1.6606497) q[0];
sx q[0];
rz(-3.1180179) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5970846) q[2];
sx q[2];
rz(-1.9015357) q[2];
sx q[2];
rz(-2.7059908) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7792369) q[1];
sx q[1];
rz(-0.52473611) q[1];
sx q[1];
rz(-1.9441821) q[1];
rz(1.9911733) q[3];
sx q[3];
rz(-1.8432518) q[3];
sx q[3];
rz(0.099893173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4219389) q[2];
sx q[2];
rz(-1.3976588) q[2];
sx q[2];
rz(-1.5092124) q[2];
rz(0.40431067) q[3];
sx q[3];
rz(-0.68250889) q[3];
sx q[3];
rz(-1.6633165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8835835) q[0];
sx q[0];
rz(-1.785935) q[0];
sx q[0];
rz(1.0255381) q[0];
rz(0.57199663) q[1];
sx q[1];
rz(-1.0943741) q[1];
sx q[1];
rz(0.62932032) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9176661) q[0];
sx q[0];
rz(-1.9804945) q[0];
sx q[0];
rz(1.2439338) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8848558) q[2];
sx q[2];
rz(-0.70054189) q[2];
sx q[2];
rz(-2.8741921) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1197966) q[1];
sx q[1];
rz(-0.52638678) q[1];
sx q[1];
rz(0.22967931) q[1];
rz(-pi) q[2];
rz(2.1719645) q[3];
sx q[3];
rz(-1.3519577) q[3];
sx q[3];
rz(0.18274433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5489674) q[2];
sx q[2];
rz(-0.36965814) q[2];
sx q[2];
rz(0.34234753) q[2];
rz(-1.7957548) q[3];
sx q[3];
rz(-1.4474409) q[3];
sx q[3];
rz(-2.9455744) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.489007) q[0];
sx q[0];
rz(-1.2359897) q[0];
sx q[0];
rz(2.956399) q[0];
rz(-1.7350896) q[1];
sx q[1];
rz(-2.0506737) q[1];
sx q[1];
rz(-1.7746183) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0104116) q[0];
sx q[0];
rz(-0.62958065) q[0];
sx q[0];
rz(-0.38139831) q[0];
rz(-pi) q[1];
x q[1];
rz(2.30079) q[2];
sx q[2];
rz(-1.7992939) q[2];
sx q[2];
rz(-2.218354) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.023932799) q[1];
sx q[1];
rz(-1.6624833) q[1];
sx q[1];
rz(1.7232643) q[1];
x q[2];
rz(-1.1432511) q[3];
sx q[3];
rz(-2.098009) q[3];
sx q[3];
rz(-0.15641071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8967445) q[2];
sx q[2];
rz(-2.5881793) q[2];
sx q[2];
rz(-2.2582167) q[2];
rz(1.4128489) q[3];
sx q[3];
rz(-0.69245517) q[3];
sx q[3];
rz(0.81378716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8900523) q[0];
sx q[0];
rz(-0.10245704) q[0];
sx q[0];
rz(-1.863377) q[0];
rz(0.037840769) q[1];
sx q[1];
rz(-2.3262639) q[1];
sx q[1];
rz(1.7657123) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1139039) q[0];
sx q[0];
rz(-1.2776327) q[0];
sx q[0];
rz(-2.3011074) q[0];
rz(-pi) q[1];
rz(1.0972294) q[2];
sx q[2];
rz(-1.5175022) q[2];
sx q[2];
rz(3.1380944) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3657276) q[1];
sx q[1];
rz(-0.15656808) q[1];
sx q[1];
rz(2.4003029) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7405628) q[3];
sx q[3];
rz(-2.8068636) q[3];
sx q[3];
rz(-1.7498121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4574796) q[2];
sx q[2];
rz(-0.71762466) q[2];
sx q[2];
rz(-0.73105556) q[2];
rz(3.030792) q[3];
sx q[3];
rz(-1.5564857) q[3];
sx q[3];
rz(0.69537648) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59657997) q[0];
sx q[0];
rz(-2.2667363) q[0];
sx q[0];
rz(-0.73356432) q[0];
rz(2.5336174) q[1];
sx q[1];
rz(-1.1939476) q[1];
sx q[1];
rz(0.2342934) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5406571) q[0];
sx q[0];
rz(-2.0458851) q[0];
sx q[0];
rz(-2.8500772) q[0];
x q[1];
rz(0.90227622) q[2];
sx q[2];
rz(-2.1736645) q[2];
sx q[2];
rz(-0.55842802) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5768847) q[1];
sx q[1];
rz(-1.3815834) q[1];
sx q[1];
rz(0.8072168) q[1];
rz(-pi) q[2];
rz(-0.9849302) q[3];
sx q[3];
rz(-1.8551991) q[3];
sx q[3];
rz(-1.5581074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0155448) q[2];
sx q[2];
rz(-1.3092821) q[2];
sx q[2];
rz(-1.7162494) q[2];
rz(-1.6783293) q[3];
sx q[3];
rz(-0.78444702) q[3];
sx q[3];
rz(-1.6459758) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54206806) q[0];
sx q[0];
rz(-2.8064089) q[0];
sx q[0];
rz(-1.19338) q[0];
rz(1.2606196) q[1];
sx q[1];
rz(-1.3648938) q[1];
sx q[1];
rz(-0.9448005) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1856857) q[0];
sx q[0];
rz(-2.3392896) q[0];
sx q[0];
rz(2.5923652) q[0];
rz(-pi) q[1];
rz(-0.94061942) q[2];
sx q[2];
rz(-2.665328) q[2];
sx q[2];
rz(-2.1911052) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24430844) q[1];
sx q[1];
rz(-0.329204) q[1];
sx q[1];
rz(0.55397482) q[1];
rz(-2.0637022) q[3];
sx q[3];
rz(-0.3993881) q[3];
sx q[3];
rz(2.3139017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.21645674) q[2];
sx q[2];
rz(-1.0711121) q[2];
sx q[2];
rz(-0.28820583) q[2];
rz(0.47973412) q[3];
sx q[3];
rz(-2.0917442) q[3];
sx q[3];
rz(-1.5079927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0297246) q[0];
sx q[0];
rz(-2.8653963) q[0];
sx q[0];
rz(-2.2286041) q[0];
rz(0.37462014) q[1];
sx q[1];
rz(-1.4034142) q[1];
sx q[1];
rz(-2.250681) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4961632) q[0];
sx q[0];
rz(-1.9053639) q[0];
sx q[0];
rz(1.3076925) q[0];
x q[1];
rz(-0.28384039) q[2];
sx q[2];
rz(-0.97728697) q[2];
sx q[2];
rz(-0.80988353) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18894698) q[1];
sx q[1];
rz(-1.2332321) q[1];
sx q[1];
rz(1.7811716) q[1];
rz(-pi) q[2];
rz(1.4952412) q[3];
sx q[3];
rz(-1.3831426) q[3];
sx q[3];
rz(2.5726266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9051819) q[2];
sx q[2];
rz(-0.85835251) q[2];
sx q[2];
rz(2.6399844) q[2];
rz(1.8524528) q[3];
sx q[3];
rz(-1.4533549) q[3];
sx q[3];
rz(-0.44617173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63507737) q[0];
sx q[0];
rz(-1.4415393) q[0];
sx q[0];
rz(-2.517979) q[0];
rz(1.1322017) q[1];
sx q[1];
rz(-0.75695801) q[1];
sx q[1];
rz(-3.0523041) q[1];
rz(2.8749309) q[2];
sx q[2];
rz(-1.4229792) q[2];
sx q[2];
rz(-0.40001043) q[2];
rz(-1.0241667) q[3];
sx q[3];
rz(-1.6129941) q[3];
sx q[3];
rz(2.2254406) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
