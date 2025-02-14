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
rz(0.74855411) q[0];
sx q[0];
rz(-1.3286123) q[0];
sx q[0];
rz(-2.7207029) q[0];
rz(2.4731877) q[1];
sx q[1];
rz(-1.547812) q[1];
sx q[1];
rz(-1.9903543) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0092888) q[0];
sx q[0];
rz(-2.0175066) q[0];
sx q[0];
rz(1.7942091) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45999476) q[2];
sx q[2];
rz(-1.056571) q[2];
sx q[2];
rz(-2.0005039) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.1681663) q[1];
sx q[1];
rz(-2.3273673) q[1];
sx q[1];
rz(2.3553215) q[1];
rz(-1.2251159) q[3];
sx q[3];
rz(-2.4596523) q[3];
sx q[3];
rz(0.038393858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0737754) q[2];
sx q[2];
rz(-0.99738085) q[2];
sx q[2];
rz(1.0214405) q[2];
rz(1.3218309) q[3];
sx q[3];
rz(-2.6363711) q[3];
sx q[3];
rz(2.5652313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8030871) q[0];
sx q[0];
rz(-1.0271238) q[0];
sx q[0];
rz(0.3013674) q[0];
rz(-1.1400247) q[1];
sx q[1];
rz(-0.81914425) q[1];
sx q[1];
rz(2.0778621) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057010896) q[0];
sx q[0];
rz(-1.8370596) q[0];
sx q[0];
rz(3.0169456) q[0];
rz(-2.930141) q[2];
sx q[2];
rz(-0.7728979) q[2];
sx q[2];
rz(2.7711902) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1769907) q[1];
sx q[1];
rz(-2.1637205) q[1];
sx q[1];
rz(-1.9511392) q[1];
x q[2];
rz(3.0228457) q[3];
sx q[3];
rz(-1.0380961) q[3];
sx q[3];
rz(-1.1632196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.86338824) q[2];
sx q[2];
rz(-0.74625838) q[2];
sx q[2];
rz(-2.9928652) q[2];
rz(-1.9637828) q[3];
sx q[3];
rz(-1.9494467) q[3];
sx q[3];
rz(-1.8671794) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3026368) q[0];
sx q[0];
rz(-1.9479072) q[0];
sx q[0];
rz(-1.0021915) q[0];
rz(-1.3305371) q[1];
sx q[1];
rz(-1.9621153) q[1];
sx q[1];
rz(0.60752121) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1992545) q[0];
sx q[0];
rz(-0.94616854) q[0];
sx q[0];
rz(-2.6074431) q[0];
rz(-pi) q[1];
rz(1.1639062) q[2];
sx q[2];
rz(-1.8259138) q[2];
sx q[2];
rz(1.6734378) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3866736) q[1];
sx q[1];
rz(-0.52820871) q[1];
sx q[1];
rz(0.22294238) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1488299) q[3];
sx q[3];
rz(-2.2053763) q[3];
sx q[3];
rz(1.5082593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9017631) q[2];
sx q[2];
rz(-2.2679057) q[2];
sx q[2];
rz(-2.152781) q[2];
rz(-1.4766988) q[3];
sx q[3];
rz(-1.3380545) q[3];
sx q[3];
rz(-0.57507676) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7779509) q[0];
sx q[0];
rz(-0.87894428) q[0];
sx q[0];
rz(-1.9212035) q[0];
rz(-1.3149423) q[1];
sx q[1];
rz(-1.6362135) q[1];
sx q[1];
rz(3.0016518) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51516384) q[0];
sx q[0];
rz(-2.4068711) q[0];
sx q[0];
rz(2.5236041) q[0];
rz(2.1084505) q[2];
sx q[2];
rz(-2.2058479) q[2];
sx q[2];
rz(-2.6376574) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.23754696) q[1];
sx q[1];
rz(-2.0260467) q[1];
sx q[1];
rz(-0.039467875) q[1];
rz(-pi) q[2];
rz(-1.9544551) q[3];
sx q[3];
rz(-1.6373489) q[3];
sx q[3];
rz(-0.47011791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6949029) q[2];
sx q[2];
rz(-2.0835154) q[2];
sx q[2];
rz(-0.13738446) q[2];
rz(-1.8047699) q[3];
sx q[3];
rz(-1.5627912) q[3];
sx q[3];
rz(0.00060753879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3985652) q[0];
sx q[0];
rz(-1.407857) q[0];
sx q[0];
rz(2.9803357) q[0];
rz(-2.297961) q[1];
sx q[1];
rz(-0.7470986) q[1];
sx q[1];
rz(1.3074494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8358993) q[0];
sx q[0];
rz(-2.3229555) q[0];
sx q[0];
rz(-2.7710415) q[0];
rz(-pi) q[1];
rz(-0.1990541) q[2];
sx q[2];
rz(-1.2110333) q[2];
sx q[2];
rz(0.86866405) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5999509) q[1];
sx q[1];
rz(-0.58990084) q[1];
sx q[1];
rz(-1.7322503) q[1];
x q[2];
rz(-2.3639115) q[3];
sx q[3];
rz(-2.2984004) q[3];
sx q[3];
rz(-0.2524337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1696986) q[2];
sx q[2];
rz(-0.43115386) q[2];
sx q[2];
rz(-0.92740721) q[2];
rz(-1.6124604) q[3];
sx q[3];
rz(-2.0404787) q[3];
sx q[3];
rz(-2.7839938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9277495) q[0];
sx q[0];
rz(-0.96927154) q[0];
sx q[0];
rz(1.7359605) q[0];
rz(-1.7270145) q[1];
sx q[1];
rz(-2.3443293) q[1];
sx q[1];
rz(2.3419211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23223497) q[0];
sx q[0];
rz(-0.063061558) q[0];
sx q[0];
rz(-1.3841649) q[0];
rz(-pi) q[1];
rz(-1.437101) q[2];
sx q[2];
rz(-1.8779199) q[2];
sx q[2];
rz(1.6973059) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.97669125) q[1];
sx q[1];
rz(-1.8215239) q[1];
sx q[1];
rz(-1.9216955) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7148083) q[3];
sx q[3];
rz(-2.7023661) q[3];
sx q[3];
rz(0.78204417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.25524011) q[2];
sx q[2];
rz(-0.73840529) q[2];
sx q[2];
rz(2.3389471) q[2];
rz(2.8211527) q[3];
sx q[3];
rz(-1.3405864) q[3];
sx q[3];
rz(1.6360487) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7771512) q[0];
sx q[0];
rz(-0.027712263) q[0];
sx q[0];
rz(-3.1112772) q[0];
rz(-1.8346571) q[1];
sx q[1];
rz(-2.0590643) q[1];
sx q[1];
rz(0.62987769) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7138853) q[0];
sx q[0];
rz(-2.6763335) q[0];
sx q[0];
rz(-0.42963877) q[0];
rz(-pi) q[1];
rz(0.81923072) q[2];
sx q[2];
rz(-0.60556817) q[2];
sx q[2];
rz(1.8964963) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4737638) q[1];
sx q[1];
rz(-2.7817717) q[1];
sx q[1];
rz(-2.5303165) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9050352) q[3];
sx q[3];
rz(-0.45522296) q[3];
sx q[3];
rz(-0.61456028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1097172) q[2];
sx q[2];
rz(-2.4745291) q[2];
sx q[2];
rz(1.9943705) q[2];
rz(0.68754292) q[3];
sx q[3];
rz(-2.5496428) q[3];
sx q[3];
rz(1.8182925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
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
rz(-2.3654093) q[0];
sx q[0];
rz(-2.9054346) q[0];
sx q[0];
rz(-2.6614406) q[0];
rz(-2.7393553) q[1];
sx q[1];
rz(-1.4996585) q[1];
sx q[1];
rz(-0.38280907) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7744336) q[0];
sx q[0];
rz(-1.2437789) q[0];
sx q[0];
rz(-2.9190382) q[0];
x q[1];
rz(0.11522861) q[2];
sx q[2];
rz(-1.8602716) q[2];
sx q[2];
rz(0.25184271) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20732197) q[1];
sx q[1];
rz(-1.2756961) q[1];
sx q[1];
rz(-1.9701824) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.945369) q[3];
sx q[3];
rz(-1.5999402) q[3];
sx q[3];
rz(-1.9118229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.39168921) q[2];
sx q[2];
rz(-0.73035556) q[2];
sx q[2];
rz(2.7833617) q[2];
rz(-0.41424888) q[3];
sx q[3];
rz(-1.0284938) q[3];
sx q[3];
rz(2.7481368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8527894) q[0];
sx q[0];
rz(-1.8599334) q[0];
sx q[0];
rz(0.417867) q[0];
rz(1.6425543) q[1];
sx q[1];
rz(-0.42048979) q[1];
sx q[1];
rz(-2.5299759) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9238511) q[0];
sx q[0];
rz(-1.3425266) q[0];
sx q[0];
rz(-1.3598674) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5081399) q[2];
sx q[2];
rz(-1.776172) q[2];
sx q[2];
rz(-2.5938671) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.043165) q[1];
sx q[1];
rz(-1.0994776) q[1];
sx q[1];
rz(2.8275675) q[1];
rz(1.4889273) q[3];
sx q[3];
rz(-1.9165526) q[3];
sx q[3];
rz(-0.053178259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9548107) q[2];
sx q[2];
rz(-1.6839226) q[2];
sx q[2];
rz(2.7853454) q[2];
rz(-0.48480836) q[3];
sx q[3];
rz(-1.7907413) q[3];
sx q[3];
rz(-0.029732186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32271069) q[0];
sx q[0];
rz(-1.1047381) q[0];
sx q[0];
rz(1.4622965) q[0];
rz(2.7541584) q[1];
sx q[1];
rz(-0.92715627) q[1];
sx q[1];
rz(2.3364054) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0793251) q[0];
sx q[0];
rz(-1.0073509) q[0];
sx q[0];
rz(-2.2939373) q[0];
x q[1];
rz(2.7323805) q[2];
sx q[2];
rz(-1.9878888) q[2];
sx q[2];
rz(-2.6155748) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.78811121) q[1];
sx q[1];
rz(-0.32057163) q[1];
sx q[1];
rz(-2.660257) q[1];
x q[2];
rz(-1.322454) q[3];
sx q[3];
rz(-1.9232009) q[3];
sx q[3];
rz(3.0737511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.33599535) q[2];
sx q[2];
rz(-2.5184641) q[2];
sx q[2];
rz(1.222329) q[2];
rz(-0.81021106) q[3];
sx q[3];
rz(-0.92456341) q[3];
sx q[3];
rz(-1.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41351086) q[0];
sx q[0];
rz(-1.5039197) q[0];
sx q[0];
rz(1.529827) q[0];
rz(-1.8660846) q[1];
sx q[1];
rz(-2.7559912) q[1];
sx q[1];
rz(-1.4082946) q[1];
rz(3.0048117) q[2];
sx q[2];
rz(-1.9521458) q[2];
sx q[2];
rz(-1.4804179) q[2];
rz(0.90295193) q[3];
sx q[3];
rz(-1.9330238) q[3];
sx q[3];
rz(-0.15180363) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
