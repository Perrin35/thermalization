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
rz(1.3746102) q[0];
sx q[0];
rz(4.1144028) q[0];
sx q[0];
rz(10.65539) q[0];
rz(-1.542701) q[1];
sx q[1];
rz(-0.35294947) q[1];
sx q[1];
rz(0.43500873) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15947882) q[0];
sx q[0];
rz(-1.6702939) q[0];
sx q[0];
rz(-1.9051108) q[0];
x q[1];
rz(-0.22763197) q[2];
sx q[2];
rz(-1.442558) q[2];
sx q[2];
rz(1.1007835) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5953345) q[1];
sx q[1];
rz(-1.2915242) q[1];
sx q[1];
rz(0.78650773) q[1];
x q[2];
rz(1.5391971) q[3];
sx q[3];
rz(-1.5571619) q[3];
sx q[3];
rz(-2.0054833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.9445442) q[2];
sx q[2];
rz(-0.41434449) q[2];
sx q[2];
rz(-0.91828263) q[2];
rz(2.7783172) q[3];
sx q[3];
rz(-1.2493635) q[3];
sx q[3];
rz(-1.8173119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0488224) q[0];
sx q[0];
rz(-0.73960441) q[0];
sx q[0];
rz(2.0059465) q[0];
rz(-2.8593235) q[1];
sx q[1];
rz(-1.5156563) q[1];
sx q[1];
rz(0.20271066) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36798635) q[0];
sx q[0];
rz(-1.3895036) q[0];
sx q[0];
rz(1.3884991) q[0];
rz(-pi) q[1];
rz(-1.0890929) q[2];
sx q[2];
rz(-1.5823477) q[2];
sx q[2];
rz(0.80977075) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2139288) q[1];
sx q[1];
rz(-1.1176372) q[1];
sx q[1];
rz(1.434668) q[1];
rz(-1.5810896) q[3];
sx q[3];
rz(-1.5655616) q[3];
sx q[3];
rz(0.43983398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89381605) q[2];
sx q[2];
rz(-0.42422366) q[2];
sx q[2];
rz(0.068537863) q[2];
rz(2.2872772) q[3];
sx q[3];
rz(-1.2645384) q[3];
sx q[3];
rz(3.1356649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0900367) q[0];
sx q[0];
rz(-2.5904901) q[0];
sx q[0];
rz(-2.0447482) q[0];
rz(-2.5966273) q[1];
sx q[1];
rz(-0.68383354) q[1];
sx q[1];
rz(2.9473238) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30667728) q[0];
sx q[0];
rz(-0.55473548) q[0];
sx q[0];
rz(-0.86443211) q[0];
rz(1.0164073) q[2];
sx q[2];
rz(-1.8570327) q[2];
sx q[2];
rz(-0.58974671) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3700876) q[1];
sx q[1];
rz(-1.9290826) q[1];
sx q[1];
rz(-0.53046988) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9220334) q[3];
sx q[3];
rz(-2.7673169) q[3];
sx q[3];
rz(0.55340278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7668309) q[2];
sx q[2];
rz(-1.7685879) q[2];
sx q[2];
rz(-2.4135446) q[2];
rz(-2.9211365) q[3];
sx q[3];
rz(-0.15153344) q[3];
sx q[3];
rz(-2.5098586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24882889) q[0];
sx q[0];
rz(-1.4735824) q[0];
sx q[0];
rz(-1.0695176) q[0];
rz(0.88691521) q[1];
sx q[1];
rz(-0.506217) q[1];
sx q[1];
rz(-1.3217529) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89267545) q[0];
sx q[0];
rz(-1.6098611) q[0];
sx q[0];
rz(-0.43401735) q[0];
x q[1];
rz(1.555205) q[2];
sx q[2];
rz(-1.8097477) q[2];
sx q[2];
rz(-2.3794501) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8458183) q[1];
sx q[1];
rz(-1.1330746) q[1];
sx q[1];
rz(-0.079016642) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6245313) q[3];
sx q[3];
rz(-1.0365465) q[3];
sx q[3];
rz(-1.8815966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6679823) q[2];
sx q[2];
rz(-2.0986291) q[2];
sx q[2];
rz(1.854151) q[2];
rz(1.0188262) q[3];
sx q[3];
rz(-2.5835218) q[3];
sx q[3];
rz(-2.4655925) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2126615) q[0];
sx q[0];
rz(-0.23290817) q[0];
sx q[0];
rz(0.88053298) q[0];
rz(-0.13678837) q[1];
sx q[1];
rz(-0.22577481) q[1];
sx q[1];
rz(-2.5313964) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5870415) q[0];
sx q[0];
rz(-0.64233883) q[0];
sx q[0];
rz(-1.5917529) q[0];
x q[1];
rz(1.4177104) q[2];
sx q[2];
rz(-1.1436074) q[2];
sx q[2];
rz(1.2059905) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6366263) q[1];
sx q[1];
rz(-1.3769048) q[1];
sx q[1];
rz(-2.8411179) q[1];
rz(-pi) q[2];
rz(0.28308308) q[3];
sx q[3];
rz(-0.45520458) q[3];
sx q[3];
rz(-0.20008034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.60787624) q[2];
sx q[2];
rz(-2.4166962) q[2];
sx q[2];
rz(1.9113212) q[2];
rz(-2.2599334) q[3];
sx q[3];
rz(-1.6917546) q[3];
sx q[3];
rz(1.3151883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6864258) q[0];
sx q[0];
rz(-1.5979586) q[0];
sx q[0];
rz(2.0879188) q[0];
rz(1.7651438) q[1];
sx q[1];
rz(-2.0339298) q[1];
sx q[1];
rz(-0.68127662) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.412144) q[0];
sx q[0];
rz(-2.2707457) q[0];
sx q[0];
rz(1.5493718) q[0];
rz(1.9176964) q[2];
sx q[2];
rz(-0.85641801) q[2];
sx q[2];
rz(-2.6831085) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6974042) q[1];
sx q[1];
rz(-0.57198524) q[1];
sx q[1];
rz(-2.0818656) q[1];
rz(3.104171) q[3];
sx q[3];
rz(-0.87914124) q[3];
sx q[3];
rz(0.057400364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0690339) q[2];
sx q[2];
rz(-0.37798887) q[2];
sx q[2];
rz(0.47522137) q[2];
rz(1.2230988) q[3];
sx q[3];
rz(-0.94448543) q[3];
sx q[3];
rz(2.9782915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.908919) q[0];
sx q[0];
rz(-0.95218807) q[0];
sx q[0];
rz(-2.9222144) q[0];
rz(-1.0064005) q[1];
sx q[1];
rz(-1.718113) q[1];
sx q[1];
rz(-1.9361608) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4426098) q[0];
sx q[0];
rz(-0.97215334) q[0];
sx q[0];
rz(-3.1151616) q[0];
rz(-pi) q[1];
rz(-1.3036626) q[2];
sx q[2];
rz(-1.1877737) q[2];
sx q[2];
rz(2.2522425) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6483375) q[1];
sx q[1];
rz(-1.770079) q[1];
sx q[1];
rz(1.5882476) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1262283) q[3];
sx q[3];
rz(-1.6939112) q[3];
sx q[3];
rz(-0.11492226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.22126234) q[2];
sx q[2];
rz(-0.39322501) q[2];
sx q[2];
rz(0.66366759) q[2];
rz(2.4237733) q[3];
sx q[3];
rz(-2.4847173) q[3];
sx q[3];
rz(2.9653911) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32888907) q[0];
sx q[0];
rz(-0.95903522) q[0];
sx q[0];
rz(-1.2411728) q[0];
rz(2.7136956) q[1];
sx q[1];
rz(-1.5792081) q[1];
sx q[1];
rz(-1.7180299) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.099402) q[0];
sx q[0];
rz(-1.4960327) q[0];
sx q[0];
rz(-0.21850234) q[0];
x q[1];
rz(-0.15870416) q[2];
sx q[2];
rz(-0.85904143) q[2];
sx q[2];
rz(-1.2325328) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9729303) q[1];
sx q[1];
rz(-0.902938) q[1];
sx q[1];
rz(-3.1144823) q[1];
x q[2];
rz(1.6070494) q[3];
sx q[3];
rz(-1.9055618) q[3];
sx q[3];
rz(0.70543766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1443892) q[2];
sx q[2];
rz(-2.453697) q[2];
sx q[2];
rz(2.8456861) q[2];
rz(-1.0812673) q[3];
sx q[3];
rz(-1.5238949) q[3];
sx q[3];
rz(2.9937939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5645771) q[0];
sx q[0];
rz(-0.58203375) q[0];
sx q[0];
rz(2.4356595) q[0];
rz(2.9544746) q[1];
sx q[1];
rz(-0.83693224) q[1];
sx q[1];
rz(-0.44609889) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4574544) q[0];
sx q[0];
rz(-1.7230526) q[0];
sx q[0];
rz(-1.5689374) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0995277) q[2];
sx q[2];
rz(-1.285811) q[2];
sx q[2];
rz(1.0302271) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.88893585) q[1];
sx q[1];
rz(-1.2351079) q[1];
sx q[1];
rz(-0.97519213) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57658617) q[3];
sx q[3];
rz(-1.4268677) q[3];
sx q[3];
rz(-1.3670539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.17911653) q[2];
sx q[2];
rz(-2.1160782) q[2];
sx q[2];
rz(-1.7402488) q[2];
rz(-2.5380747) q[3];
sx q[3];
rz(-2.9602435) q[3];
sx q[3];
rz(0.13874273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.03183455) q[0];
sx q[0];
rz(-1.4048445) q[0];
sx q[0];
rz(-3.0314714) q[0];
rz(-1.8203863) q[1];
sx q[1];
rz(-0.92480129) q[1];
sx q[1];
rz(-0.22470156) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53407828) q[0];
sx q[0];
rz(-0.67994962) q[0];
sx q[0];
rz(0.91331608) q[0];
rz(-2.9797033) q[2];
sx q[2];
rz(-1.8542552) q[2];
sx q[2];
rz(-3.0155011) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.68606317) q[1];
sx q[1];
rz(-1.7357329) q[1];
sx q[1];
rz(0.73442028) q[1];
rz(-pi) q[2];
rz(1.4175156) q[3];
sx q[3];
rz(-0.17440052) q[3];
sx q[3];
rz(2.0896951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3118185) q[2];
sx q[2];
rz(-0.31121397) q[2];
sx q[2];
rz(-0.52660006) q[2];
rz(-0.38463587) q[3];
sx q[3];
rz(-1.2472109) q[3];
sx q[3];
rz(-0.23469901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93152355) q[0];
sx q[0];
rz(-2.7061404) q[0];
sx q[0];
rz(2.7192116) q[0];
rz(2.3474563) q[1];
sx q[1];
rz(-1.7105449) q[1];
sx q[1];
rz(1.7078043) q[1];
rz(0.83195291) q[2];
sx q[2];
rz(-0.83870287) q[2];
sx q[2];
rz(2.4626682) q[2];
rz(-0.41153367) q[3];
sx q[3];
rz(-1.2306662) q[3];
sx q[3];
rz(-1.6231404) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
