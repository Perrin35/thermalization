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
rz(0.39245519) q[0];
sx q[0];
rz(-0.28372228) q[0];
sx q[0];
rz(1.1878045) q[0];
rz(-0.44161931) q[1];
sx q[1];
rz(-1.8752357) q[1];
sx q[1];
rz(-1.9904354) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6080127) q[0];
sx q[0];
rz(-1.5488429) q[0];
sx q[0];
rz(1.4826464) q[0];
rz(1.2032937) q[2];
sx q[2];
rz(-1.1137059) q[2];
sx q[2];
rz(0.45705308) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.44635998) q[1];
sx q[1];
rz(-2.6962013) q[1];
sx q[1];
rz(-2.8421418) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0813484) q[3];
sx q[3];
rz(-0.57751211) q[3];
sx q[3];
rz(-0.58693991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2426131) q[2];
sx q[2];
rz(-0.8967163) q[2];
sx q[2];
rz(0.29956451) q[2];
rz(-2.7859531) q[3];
sx q[3];
rz(-2.1483597) q[3];
sx q[3];
rz(-2.7049844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3886609) q[0];
sx q[0];
rz(-2.5932073) q[0];
sx q[0];
rz(-0.27632982) q[0];
rz(0.07946864) q[1];
sx q[1];
rz(-2.6092968) q[1];
sx q[1];
rz(-2.6452549) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2329335) q[0];
sx q[0];
rz(-2.56224) q[0];
sx q[0];
rz(1.1003554) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9251516) q[2];
sx q[2];
rz(-0.34578824) q[2];
sx q[2];
rz(-1.3160637) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6092094) q[1];
sx q[1];
rz(-0.84621284) q[1];
sx q[1];
rz(2.0128472) q[1];
rz(-pi) q[2];
rz(1.5509505) q[3];
sx q[3];
rz(-2.9796585) q[3];
sx q[3];
rz(-1.9753044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3695099) q[2];
sx q[2];
rz(-1.3630867) q[2];
sx q[2];
rz(-2.962964) q[2];
rz(1.5313088) q[3];
sx q[3];
rz(-0.9261927) q[3];
sx q[3];
rz(-2.3791651) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4302706) q[0];
sx q[0];
rz(-2.0549759) q[0];
sx q[0];
rz(1.7311199) q[0];
rz(1.3871644) q[1];
sx q[1];
rz(-1.4924563) q[1];
sx q[1];
rz(14/(3*pi)) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1324155) q[0];
sx q[0];
rz(-2.2542291) q[0];
sx q[0];
rz(-1.528549) q[0];
rz(-2.335821) q[2];
sx q[2];
rz(-1.2670332) q[2];
sx q[2];
rz(-2.8250776) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.3988603) q[1];
sx q[1];
rz(-1.1452864) q[1];
sx q[1];
rz(2.8102283) q[1];
rz(-pi) q[2];
rz(2.0250506) q[3];
sx q[3];
rz(-0.71703934) q[3];
sx q[3];
rz(1.7771135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2014655) q[2];
sx q[2];
rz(-2.1374233) q[2];
sx q[2];
rz(2.6489769) q[2];
rz(2.7228739) q[3];
sx q[3];
rz(-1.0227572) q[3];
sx q[3];
rz(0.59188265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017460499) q[0];
sx q[0];
rz(-3.0191665) q[0];
sx q[0];
rz(0.19805743) q[0];
rz(-0.65912229) q[1];
sx q[1];
rz(-2.5810869) q[1];
sx q[1];
rz(-2.9445599) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7871088) q[0];
sx q[0];
rz(-2.1012044) q[0];
sx q[0];
rz(1.5317937) q[0];
rz(-pi) q[1];
rz(0.87984271) q[2];
sx q[2];
rz(-2.7759399) q[2];
sx q[2];
rz(-1.8118389) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9024762) q[1];
sx q[1];
rz(-0.65985876) q[1];
sx q[1];
rz(0.59999864) q[1];
rz(-1.0258338) q[3];
sx q[3];
rz(-1.539388) q[3];
sx q[3];
rz(-1.1823013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.78843242) q[2];
sx q[2];
rz(-1.8400695) q[2];
sx q[2];
rz(1.8678467) q[2];
rz(-1.5431917) q[3];
sx q[3];
rz(-1.1929932) q[3];
sx q[3];
rz(2.8852692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0850328) q[0];
sx q[0];
rz(-1.4075449) q[0];
sx q[0];
rz(-3.0388167) q[0];
rz(-0.13725266) q[1];
sx q[1];
rz(-2.6922701) q[1];
sx q[1];
rz(-1.4385673) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.611127) q[0];
sx q[0];
rz(-2.1411863) q[0];
sx q[0];
rz(0.80369759) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0497018) q[2];
sx q[2];
rz(-2.6767807) q[2];
sx q[2];
rz(1.6824695) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2074429) q[1];
sx q[1];
rz(-2.8420919) q[1];
sx q[1];
rz(-1.0494166) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6827312) q[3];
sx q[3];
rz(-1.4482968) q[3];
sx q[3];
rz(-1.4089597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.30456257) q[2];
sx q[2];
rz(-1.4745159) q[2];
sx q[2];
rz(2.4414818) q[2];
rz(-2.2434798) q[3];
sx q[3];
rz(-2.5105295) q[3];
sx q[3];
rz(1.496544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20823088) q[0];
sx q[0];
rz(-2.5550186) q[0];
sx q[0];
rz(-1.5775648) q[0];
rz(-0.6598407) q[1];
sx q[1];
rz(-0.78958646) q[1];
sx q[1];
rz(0.97445625) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81466573) q[0];
sx q[0];
rz(-1.5626994) q[0];
sx q[0];
rz(-0.011122313) q[0];
x q[1];
rz(-2.8235213) q[2];
sx q[2];
rz(-2.632395) q[2];
sx q[2];
rz(-1.4105907) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2657915) q[1];
sx q[1];
rz(-2.6398924) q[1];
sx q[1];
rz(0.022345101) q[1];
rz(-pi) q[2];
rz(0.10147126) q[3];
sx q[3];
rz(-0.54059282) q[3];
sx q[3];
rz(-2.0848047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.66880354) q[2];
sx q[2];
rz(-0.69837022) q[2];
sx q[2];
rz(-0.93639708) q[2];
rz(-2.2514553) q[3];
sx q[3];
rz(-1.7809159) q[3];
sx q[3];
rz(0.20588017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5278006) q[0];
sx q[0];
rz(-2.9449936) q[0];
sx q[0];
rz(0.11372456) q[0];
rz(1.744489) q[1];
sx q[1];
rz(-2.0753658) q[1];
sx q[1];
rz(0.024959175) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67201383) q[0];
sx q[0];
rz(-0.22525283) q[0];
sx q[0];
rz(1.3242871) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3106842) q[2];
sx q[2];
rz(-0.94365135) q[2];
sx q[2];
rz(-2.8994034) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.61103454) q[1];
sx q[1];
rz(-1.636158) q[1];
sx q[1];
rz(-2.2863131) q[1];
rz(-pi) q[2];
rz(-1.5146144) q[3];
sx q[3];
rz(-1.3732519) q[3];
sx q[3];
rz(0.88314013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52594841) q[2];
sx q[2];
rz(-1.7082461) q[2];
sx q[2];
rz(1.0053267) q[2];
rz(0.89056906) q[3];
sx q[3];
rz(-1.7361448) q[3];
sx q[3];
rz(2.4864206) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6152076) q[0];
sx q[0];
rz(-0.40281519) q[0];
sx q[0];
rz(-2.030754) q[0];
rz(3.0530744) q[1];
sx q[1];
rz(-0.42461494) q[1];
sx q[1];
rz(1.7875338) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18830472) q[0];
sx q[0];
rz(-2.3416356) q[0];
sx q[0];
rz(-0.15733457) q[0];
rz(-1.8432003) q[2];
sx q[2];
rz(-1.4227616) q[2];
sx q[2];
rz(1.7065976) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.718487) q[1];
sx q[1];
rz(-0.25782789) q[1];
sx q[1];
rz(-2.8400374) q[1];
rz(-pi) q[2];
rz(-0.60175394) q[3];
sx q[3];
rz(-0.85769282) q[3];
sx q[3];
rz(-2.8181638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84305772) q[2];
sx q[2];
rz(-0.83710805) q[2];
sx q[2];
rz(-0.38541547) q[2];
rz(-2.391583) q[3];
sx q[3];
rz(-2.1842726) q[3];
sx q[3];
rz(0.065940417) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0176004) q[0];
sx q[0];
rz(-1.0753068) q[0];
sx q[0];
rz(-2.3093834) q[0];
rz(-2.0596793) q[1];
sx q[1];
rz(-2.3257906) q[1];
sx q[1];
rz(1.6197416) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2913301) q[0];
sx q[0];
rz(-1.1801774) q[0];
sx q[0];
rz(-2.0807092) q[0];
x q[1];
rz(-2.3004901) q[2];
sx q[2];
rz(-0.37174598) q[2];
sx q[2];
rz(1.8018307) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.950022) q[1];
sx q[1];
rz(-2.2723928) q[1];
sx q[1];
rz(0.8002886) q[1];
rz(-pi) q[2];
rz(2.2275739) q[3];
sx q[3];
rz(-1.7631712) q[3];
sx q[3];
rz(-3.0466444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.99721432) q[2];
sx q[2];
rz(-1.245456) q[2];
sx q[2];
rz(1.2388371) q[2];
rz(-1.8501806) q[3];
sx q[3];
rz(-1.2715128) q[3];
sx q[3];
rz(3.0102357) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0111897) q[0];
sx q[0];
rz(-2.1881396) q[0];
sx q[0];
rz(0.74380547) q[0];
rz(-0.047608308) q[1];
sx q[1];
rz(-2.036939) q[1];
sx q[1];
rz(1.1415175) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6358546) q[0];
sx q[0];
rz(-0.72162745) q[0];
sx q[0];
rz(-0.82765915) q[0];
rz(-pi) q[1];
rz(0.6170527) q[2];
sx q[2];
rz(-1.8931555) q[2];
sx q[2];
rz(0.53215357) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9824578) q[1];
sx q[1];
rz(-1.1679113) q[1];
sx q[1];
rz(2.3680658) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64472736) q[3];
sx q[3];
rz(-1.4253741) q[3];
sx q[3];
rz(-2.4032088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.49012524) q[2];
sx q[2];
rz(-0.23568025) q[2];
sx q[2];
rz(-2.8362595) q[2];
rz(1.4281645) q[3];
sx q[3];
rz(-1.3729264) q[3];
sx q[3];
rz(1.292424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3931428) q[0];
sx q[0];
rz(-0.23101692) q[0];
sx q[0];
rz(1.1363181) q[0];
rz(-2.582386) q[1];
sx q[1];
rz(-1.8103841) q[1];
sx q[1];
rz(0.24787535) q[1];
rz(2.4743773) q[2];
sx q[2];
rz(-1.8493152) q[2];
sx q[2];
rz(2.2753459) q[2];
rz(-0.1189726) q[3];
sx q[3];
rz(-1.705022) q[3];
sx q[3];
rz(-2.6172887) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
