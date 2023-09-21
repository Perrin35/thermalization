OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5947333) q[0];
sx q[0];
rz(-1.5164627) q[0];
sx q[0];
rz(-2.8773142) q[0];
rz(2.1677986) q[1];
sx q[1];
rz(-1.9314613) q[1];
sx q[1];
rz(-0.73524737) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11156946) q[0];
sx q[0];
rz(-0.91696793) q[0];
sx q[0];
rz(-1.3841188) q[0];
rz(-pi) q[1];
rz(-1.0636319) q[2];
sx q[2];
rz(-2.1858366) q[2];
sx q[2];
rz(1.0539436) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4525675) q[1];
sx q[1];
rz(-1.8568294) q[1];
sx q[1];
rz(0.16098117) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5611476) q[3];
sx q[3];
rz(-1.547303) q[3];
sx q[3];
rz(-1.98222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66951093) q[2];
sx q[2];
rz(-1.8410204) q[2];
sx q[2];
rz(2.0377339) q[2];
rz(-1.8707229) q[3];
sx q[3];
rz(-1.2277675) q[3];
sx q[3];
rz(2.8675573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0274149) q[0];
sx q[0];
rz(-0.52863055) q[0];
sx q[0];
rz(2.7052178) q[0];
rz(-2.6787058) q[1];
sx q[1];
rz(-1.0375689) q[1];
sx q[1];
rz(-0.26611051) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14861815) q[0];
sx q[0];
rz(-2.2668112) q[0];
sx q[0];
rz(0.50177411) q[0];
rz(-pi) q[1];
rz(-0.99533178) q[2];
sx q[2];
rz(-1.7082214) q[2];
sx q[2];
rz(0.95057887) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3937711) q[1];
sx q[1];
rz(-2.1854679) q[1];
sx q[1];
rz(0.6786896) q[1];
x q[2];
rz(-0.18272419) q[3];
sx q[3];
rz(-0.97191873) q[3];
sx q[3];
rz(0.11174186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6767072) q[2];
sx q[2];
rz(-1.2767982) q[2];
sx q[2];
rz(-0.51149386) q[2];
rz(-0.809508) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(0.28765837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4354316) q[0];
sx q[0];
rz(-1.5478739) q[0];
sx q[0];
rz(-2.2128552) q[0];
rz(1.4061032) q[1];
sx q[1];
rz(-2.4415253) q[1];
sx q[1];
rz(-1.7944638) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4192393) q[0];
sx q[0];
rz(-1.1142715) q[0];
sx q[0];
rz(1.8629575) q[0];
rz(2.0688829) q[2];
sx q[2];
rz(-0.32215298) q[2];
sx q[2];
rz(1.0108394) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82145065) q[1];
sx q[1];
rz(-2.4648033) q[1];
sx q[1];
rz(-1.0980781) q[1];
rz(-1.3942765) q[3];
sx q[3];
rz(-2.4268097) q[3];
sx q[3];
rz(2.178758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9946263) q[2];
sx q[2];
rz(-1.7549843) q[2];
sx q[2];
rz(1.6195126) q[2];
rz(-0.26432031) q[3];
sx q[3];
rz(-1.0364573) q[3];
sx q[3];
rz(-0.49595293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79214823) q[0];
sx q[0];
rz(-1.148372) q[0];
sx q[0];
rz(-0.96570063) q[0];
rz(2.4194338) q[1];
sx q[1];
rz(-1.5042217) q[1];
sx q[1];
rz(0.55975634) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20630079) q[0];
sx q[0];
rz(-1.1456523) q[0];
sx q[0];
rz(1.6598808) q[0];
rz(-2.431805) q[2];
sx q[2];
rz(-2.1927532) q[2];
sx q[2];
rz(-0.70993916) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8371115) q[1];
sx q[1];
rz(-1.5251535) q[1];
sx q[1];
rz(-0.95506217) q[1];
rz(-pi) q[2];
rz(2.3523931) q[3];
sx q[3];
rz(-1.4025098) q[3];
sx q[3];
rz(0.60889739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7136148) q[2];
sx q[2];
rz(-1.6327991) q[2];
sx q[2];
rz(1.4245865) q[2];
rz(-0.26040855) q[3];
sx q[3];
rz(-1.3924761) q[3];
sx q[3];
rz(-0.4310472) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.150862) q[0];
sx q[0];
rz(-1.6435511) q[0];
sx q[0];
rz(0.36636233) q[0];
rz(1.5461961) q[1];
sx q[1];
rz(-0.55391824) q[1];
sx q[1];
rz(-0.34367925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2491964) q[0];
sx q[0];
rz(-0.93660347) q[0];
sx q[0];
rz(-2.2618494) q[0];
x q[1];
rz(-2.5630066) q[2];
sx q[2];
rz(-0.95539504) q[2];
sx q[2];
rz(1.7631284) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5695565) q[1];
sx q[1];
rz(-2.0680288) q[1];
sx q[1];
rz(-3.0576502) q[1];
rz(-pi) q[2];
x q[2];
rz(3.120317) q[3];
sx q[3];
rz(-2.2145503) q[3];
sx q[3];
rz(-1.4847886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7498103) q[2];
sx q[2];
rz(-0.85931531) q[2];
sx q[2];
rz(-3.0878477) q[2];
rz(-1.7371477) q[3];
sx q[3];
rz(-0.497538) q[3];
sx q[3];
rz(0.29156175) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4869726) q[0];
sx q[0];
rz(-0.64990652) q[0];
sx q[0];
rz(-0.75575954) q[0];
rz(3.1164363) q[1];
sx q[1];
rz(-0.92725602) q[1];
sx q[1];
rz(2.8818534) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.040859) q[0];
sx q[0];
rz(-1.9319527) q[0];
sx q[0];
rz(-0.33613899) q[0];
rz(1.9070712) q[2];
sx q[2];
rz(-1.5894801) q[2];
sx q[2];
rz(2.3980354) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7139587) q[1];
sx q[1];
rz(-0.40363388) q[1];
sx q[1];
rz(-2.1778818) q[1];
rz(-pi) q[2];
rz(-0.017976956) q[3];
sx q[3];
rz(-2.020105) q[3];
sx q[3];
rz(2.9177005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6529237) q[2];
sx q[2];
rz(-0.80604625) q[2];
sx q[2];
rz(0.10061131) q[2];
rz(2.9595024) q[3];
sx q[3];
rz(-0.87825769) q[3];
sx q[3];
rz(-1.7939059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9942193) q[0];
sx q[0];
rz(-1.4202776) q[0];
sx q[0];
rz(2.8314262) q[0];
rz(0.50225964) q[1];
sx q[1];
rz(-0.52400932) q[1];
sx q[1];
rz(-0.60595864) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2180442) q[0];
sx q[0];
rz(-0.96203066) q[0];
sx q[0];
rz(-1.2952842) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29166834) q[2];
sx q[2];
rz(-1.7632315) q[2];
sx q[2];
rz(2.8067436) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.052918606) q[1];
sx q[1];
rz(-2.5175736) q[1];
sx q[1];
rz(-0.96159972) q[1];
rz(2.3378387) q[3];
sx q[3];
rz(-2.440212) q[3];
sx q[3];
rz(-0.758981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9014152) q[2];
sx q[2];
rz(-1.9577953) q[2];
sx q[2];
rz(-0.85285464) q[2];
rz(1.3700221) q[3];
sx q[3];
rz(-1.6829237) q[3];
sx q[3];
rz(0.20496932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9119499) q[0];
sx q[0];
rz(-2.5456173) q[0];
sx q[0];
rz(1.6802616) q[0];
rz(1.4029067) q[1];
sx q[1];
rz(-0.97424126) q[1];
sx q[1];
rz(-0.064037474) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8768351) q[0];
sx q[0];
rz(-1.7581853) q[0];
sx q[0];
rz(-2.1894356) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30631752) q[2];
sx q[2];
rz(-1.0848572) q[2];
sx q[2];
rz(-1.475032) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9248326) q[1];
sx q[1];
rz(-1.9257345) q[1];
sx q[1];
rz(1.6096398) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1021032) q[3];
sx q[3];
rz(-1.8825304) q[3];
sx q[3];
rz(-1.5962275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0503851) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(-1.5554265) q[2];
rz(-2.2533916) q[3];
sx q[3];
rz(-1.2398088) q[3];
sx q[3];
rz(0.92938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.14426194) q[0];
sx q[0];
rz(-2.0781131) q[0];
sx q[0];
rz(-1.7096747) q[0];
rz(-2.572708) q[1];
sx q[1];
rz(-0.535393) q[1];
sx q[1];
rz(2.0137537) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9040684) q[0];
sx q[0];
rz(-1.5547817) q[0];
sx q[0];
rz(2.9928656) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1986507) q[2];
sx q[2];
rz(-1.1254416) q[2];
sx q[2];
rz(-3.1415591) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7833264) q[1];
sx q[1];
rz(-1.3291385) q[1];
sx q[1];
rz(3.1296455) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7752152) q[3];
sx q[3];
rz(-2.0994224) q[3];
sx q[3];
rz(2.5322994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.52788064) q[2];
sx q[2];
rz(-2.2884559) q[2];
sx q[2];
rz(2.9679427) q[2];
rz(0.33637834) q[3];
sx q[3];
rz(-1.9189546) q[3];
sx q[3];
rz(-2.7500847) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4353452) q[0];
sx q[0];
rz(-0.58723891) q[0];
sx q[0];
rz(1.6760814) q[0];
rz(-0.82410518) q[1];
sx q[1];
rz(-1.5287639) q[1];
sx q[1];
rz(0.5724268) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6515598) q[0];
sx q[0];
rz(-0.76793725) q[0];
sx q[0];
rz(0.72816531) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0296713) q[2];
sx q[2];
rz(-1.3238812) q[2];
sx q[2];
rz(1.0690451) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4116284) q[1];
sx q[1];
rz(-2.8169544) q[1];
sx q[1];
rz(-2.4050729) q[1];
x q[2];
rz(1.8626067) q[3];
sx q[3];
rz(-2.2581873) q[3];
sx q[3];
rz(-0.22112267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77999014) q[2];
sx q[2];
rz(-2.6670691) q[2];
sx q[2];
rz(2.6043716) q[2];
rz(-1.0572664) q[3];
sx q[3];
rz(-2.2500762) q[3];
sx q[3];
rz(2.4479772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2873516) q[0];
sx q[0];
rz(-1.9762522) q[0];
sx q[0];
rz(1.5594788) q[0];
rz(2.6782425) q[1];
sx q[1];
rz(-0.87711038) q[1];
sx q[1];
rz(-1.6323485) q[1];
rz(-1.4434545) q[2];
sx q[2];
rz(-0.6586532) q[2];
sx q[2];
rz(-2.3011343) q[2];
rz(1.6871917) q[3];
sx q[3];
rz(-0.4909066) q[3];
sx q[3];
rz(2.7248513) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];