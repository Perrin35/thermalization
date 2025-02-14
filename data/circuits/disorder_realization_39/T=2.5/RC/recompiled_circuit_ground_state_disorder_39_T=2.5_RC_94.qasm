OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7118536) q[0];
sx q[0];
rz(-1.1450333) q[0];
sx q[0];
rz(1.0166919) q[0];
rz(-0.84438762) q[1];
sx q[1];
rz(-0.8165741) q[1];
sx q[1];
rz(-2.5752697) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1745464) q[0];
sx q[0];
rz(-0.61515795) q[0];
sx q[0];
rz(-1.4021723) q[0];
rz(3.064134) q[2];
sx q[2];
rz(-1.768154) q[2];
sx q[2];
rz(3.0826621) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29222088) q[1];
sx q[1];
rz(-1.9927036) q[1];
sx q[1];
rz(1.6590236) q[1];
rz(-pi) q[2];
rz(-2.8352691) q[3];
sx q[3];
rz(-2.5342083) q[3];
sx q[3];
rz(2.4709783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5113968) q[2];
sx q[2];
rz(-2.2679195) q[2];
sx q[2];
rz(-2.013618) q[2];
rz(1.6389716) q[3];
sx q[3];
rz(-2.2962544) q[3];
sx q[3];
rz(-0.42475167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2750435) q[0];
sx q[0];
rz(-0.4568704) q[0];
sx q[0];
rz(-0.04091111) q[0];
rz(0.54967898) q[1];
sx q[1];
rz(-0.37893852) q[1];
sx q[1];
rz(-0.080370195) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98438063) q[0];
sx q[0];
rz(-1.5665905) q[0];
sx q[0];
rz(2.4192817) q[0];
x q[1];
rz(-0.92795579) q[2];
sx q[2];
rz(-1.11737) q[2];
sx q[2];
rz(-0.7446028) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7521022) q[1];
sx q[1];
rz(-2.6536149) q[1];
sx q[1];
rz(1.4533024) q[1];
x q[2];
rz(-1.188813) q[3];
sx q[3];
rz(-0.63900164) q[3];
sx q[3];
rz(-1.5843713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.72767672) q[2];
sx q[2];
rz(-0.84299403) q[2];
sx q[2];
rz(-2.8576039) q[2];
rz(1.5208288) q[3];
sx q[3];
rz(-1.5903383) q[3];
sx q[3];
rz(-2.3782597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3182227) q[0];
sx q[0];
rz(-0.1544054) q[0];
sx q[0];
rz(2.7247317) q[0];
rz(0.93337026) q[1];
sx q[1];
rz(-0.88637543) q[1];
sx q[1];
rz(-1.0881759) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8899208) q[0];
sx q[0];
rz(-1.3265189) q[0];
sx q[0];
rz(2.1273462) q[0];
rz(-2.5809885) q[2];
sx q[2];
rz(-1.5365155) q[2];
sx q[2];
rz(0.83002582) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0644898) q[1];
sx q[1];
rz(-1.2716115) q[1];
sx q[1];
rz(-1.283518) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.635664) q[3];
sx q[3];
rz(-1.2838351) q[3];
sx q[3];
rz(0.6395517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1209968) q[2];
sx q[2];
rz(-2.4211297) q[2];
sx q[2];
rz(0.32568112) q[2];
rz(-0.39564141) q[3];
sx q[3];
rz(-0.49042693) q[3];
sx q[3];
rz(0.71780786) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91738236) q[0];
sx q[0];
rz(-0.93467394) q[0];
sx q[0];
rz(3.000946) q[0];
rz(0.39528254) q[1];
sx q[1];
rz(-1.5970767) q[1];
sx q[1];
rz(-2.7093754) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2518766) q[0];
sx q[0];
rz(-1.1036627) q[0];
sx q[0];
rz(0.45527311) q[0];
x q[1];
rz(-1.0645116) q[2];
sx q[2];
rz(-1.6461275) q[2];
sx q[2];
rz(0.34092262) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.144369) q[1];
sx q[1];
rz(-2.8565116) q[1];
sx q[1];
rz(2.8530733) q[1];
rz(0.95000259) q[3];
sx q[3];
rz(-2.0696313) q[3];
sx q[3];
rz(-1.0390337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29493368) q[2];
sx q[2];
rz(-2.0294971) q[2];
sx q[2];
rz(0.46889949) q[2];
rz(2.9851959) q[3];
sx q[3];
rz(-2.1726435) q[3];
sx q[3];
rz(1.2676574) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073931996) q[0];
sx q[0];
rz(-2.0079375) q[0];
sx q[0];
rz(2.3611948) q[0];
rz(-0.91267768) q[1];
sx q[1];
rz(-0.78881216) q[1];
sx q[1];
rz(1.3891634) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75994223) q[0];
sx q[0];
rz(-1.3033322) q[0];
sx q[0];
rz(2.777369) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19540968) q[2];
sx q[2];
rz(-0.59048803) q[2];
sx q[2];
rz(-0.54412847) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.011574419) q[1];
sx q[1];
rz(-0.60638529) q[1];
sx q[1];
rz(-2.9467086) q[1];
x q[2];
rz(2.2551304) q[3];
sx q[3];
rz(-1.2008948) q[3];
sx q[3];
rz(0.97742888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3461561) q[2];
sx q[2];
rz(-2.4861591) q[2];
sx q[2];
rz(0.1795086) q[2];
rz(-1.6620212) q[3];
sx q[3];
rz(-2.447465) q[3];
sx q[3];
rz(0.0050541335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805304) q[0];
sx q[0];
rz(-1.5164277) q[0];
sx q[0];
rz(1.6777212) q[0];
rz(3.1278817) q[1];
sx q[1];
rz(-1.6746215) q[1];
sx q[1];
rz(-0.94508583) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6701407) q[0];
sx q[0];
rz(-1.2218804) q[0];
sx q[0];
rz(-0.80100153) q[0];
rz(2.5588972) q[2];
sx q[2];
rz(-0.36755633) q[2];
sx q[2];
rz(-2.3876397) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3774776) q[1];
sx q[1];
rz(-0.65122737) q[1];
sx q[1];
rz(0.73173827) q[1];
rz(-2.2919214) q[3];
sx q[3];
rz(-2.4115218) q[3];
sx q[3];
rz(1.3645614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1549687) q[2];
sx q[2];
rz(-2.0782317) q[2];
sx q[2];
rz(1.0490136) q[2];
rz(1.6074041) q[3];
sx q[3];
rz(-1.0020071) q[3];
sx q[3];
rz(-1.3656507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55117115) q[0];
sx q[0];
rz(-0.48749247) q[0];
sx q[0];
rz(-0.76229873) q[0];
rz(2.6679692) q[1];
sx q[1];
rz(-1.381424) q[1];
sx q[1];
rz(-2.2332938) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0001091) q[0];
sx q[0];
rz(-2.3628652) q[0];
sx q[0];
rz(1.578192) q[0];
rz(-pi) q[1];
rz(-2.2328159) q[2];
sx q[2];
rz(-0.27784608) q[2];
sx q[2];
rz(-1.8709594) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6544853) q[1];
sx q[1];
rz(-1.5562378) q[1];
sx q[1];
rz(1.8076987) q[1];
rz(-pi) q[2];
rz(1.1486111) q[3];
sx q[3];
rz(-3.1109538) q[3];
sx q[3];
rz(2.2922761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.121754) q[2];
sx q[2];
rz(-2.8826931) q[2];
sx q[2];
rz(0.961595) q[2];
rz(-0.52729765) q[3];
sx q[3];
rz(-1.3798102) q[3];
sx q[3];
rz(-2.5808047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90802646) q[0];
sx q[0];
rz(-2.5416424) q[0];
sx q[0];
rz(2.7954234) q[0];
rz(1.2973805) q[1];
sx q[1];
rz(-0.66645122) q[1];
sx q[1];
rz(-0.60874879) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70102126) q[0];
sx q[0];
rz(-2.2445569) q[0];
sx q[0];
rz(3.0340723) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66242156) q[2];
sx q[2];
rz(-2.8067744) q[2];
sx q[2];
rz(-1.3864558) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4989165) q[1];
sx q[1];
rz(-1.3060946) q[1];
sx q[1];
rz(1.4993673) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5777588) q[3];
sx q[3];
rz(-2.2632416) q[3];
sx q[3];
rz(-0.13217029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.79494563) q[2];
sx q[2];
rz(-2.5865159) q[2];
sx q[2];
rz(-2.4897599) q[2];
rz(0.73976222) q[3];
sx q[3];
rz(-1.0659822) q[3];
sx q[3];
rz(2.2447926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5480492) q[0];
sx q[0];
rz(-2.6683922) q[0];
sx q[0];
rz(-2.5439673) q[0];
rz(-1.8400486) q[1];
sx q[1];
rz(-1.5751244) q[1];
sx q[1];
rz(0.32599932) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4833925) q[0];
sx q[0];
rz(-2.3478386) q[0];
sx q[0];
rz(1.9824337) q[0];
rz(-0.32749657) q[2];
sx q[2];
rz(-1.8153037) q[2];
sx q[2];
rz(1.2933019) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.9218775) q[1];
sx q[1];
rz(-1.829521) q[1];
sx q[1];
rz(-2.2398276) q[1];
x q[2];
rz(0.80504412) q[3];
sx q[3];
rz(-0.16318233) q[3];
sx q[3];
rz(1.9022579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7035383) q[2];
sx q[2];
rz(-1.9376829) q[2];
sx q[2];
rz(0.34029141) q[2];
rz(1.641364) q[3];
sx q[3];
rz(-0.54268018) q[3];
sx q[3];
rz(-0.59804183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21094766) q[0];
sx q[0];
rz(-2.0238545) q[0];
sx q[0];
rz(-3.0979544) q[0];
rz(-2.4234407) q[1];
sx q[1];
rz(-1.8225881) q[1];
sx q[1];
rz(1.7125548) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14075101) q[0];
sx q[0];
rz(-1.4006097) q[0];
sx q[0];
rz(2.9508181) q[0];
rz(0.67769717) q[2];
sx q[2];
rz(-1.4645828) q[2];
sx q[2];
rz(-2.8863557) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3188827) q[1];
sx q[1];
rz(-1.4850933) q[1];
sx q[1];
rz(2.6515168) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7765482) q[3];
sx q[3];
rz(-1.2170346) q[3];
sx q[3];
rz(-0.34788528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67937294) q[2];
sx q[2];
rz(-2.2023109) q[2];
sx q[2];
rz(-2.3636554) q[2];
rz(2.098162) q[3];
sx q[3];
rz(-1.611004) q[3];
sx q[3];
rz(-1.2354318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048024561) q[0];
sx q[0];
rz(-2.1401736) q[0];
sx q[0];
rz(-2.370136) q[0];
rz(-1.7780766) q[1];
sx q[1];
rz(-1.9693146) q[1];
sx q[1];
rz(1.6065425) q[1];
rz(-2.4189998) q[2];
sx q[2];
rz(-1.9836139) q[2];
sx q[2];
rz(0.2000533) q[2];
rz(-0.4637152) q[3];
sx q[3];
rz(-2.188893) q[3];
sx q[3];
rz(-2.7617674) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
