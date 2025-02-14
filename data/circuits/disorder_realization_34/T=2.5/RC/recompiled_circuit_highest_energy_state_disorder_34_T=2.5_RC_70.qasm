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
rz(1.1512383) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4822455) q[0];
sx q[0];
rz(-1.3696241) q[0];
sx q[0];
rz(0.45659275) q[0];
x q[1];
rz(-2.6815979) q[2];
sx q[2];
rz(-1.056571) q[2];
sx q[2];
rz(1.1410888) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1382153) q[1];
sx q[1];
rz(-1.0312407) q[1];
sx q[1];
rz(0.92745499) q[1];
x q[2];
rz(-0.26845308) q[3];
sx q[3];
rz(-2.2055948) q[3];
sx q[3];
rz(0.47273794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0737754) q[2];
sx q[2];
rz(-0.99738085) q[2];
sx q[2];
rz(-1.0214405) q[2];
rz(-1.8197618) q[3];
sx q[3];
rz(-0.50522155) q[3];
sx q[3];
rz(-2.5652313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8030871) q[0];
sx q[0];
rz(-1.0271238) q[0];
sx q[0];
rz(-2.8402253) q[0];
rz(1.1400247) q[1];
sx q[1];
rz(-0.81914425) q[1];
sx q[1];
rz(-2.0778621) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6401705) q[0];
sx q[0];
rz(-2.8482264) q[0];
sx q[0];
rz(1.1430995) q[0];
rz(-pi) q[1];
rz(-2.3799494) q[2];
sx q[2];
rz(-1.4237262) q[2];
sx q[2];
rz(-2.0936793) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7980405) q[1];
sx q[1];
rz(-0.69188373) q[1];
sx q[1];
rz(-2.6380098) q[1];
rz(-pi) q[2];
rz(1.7691019) q[3];
sx q[3];
rz(-2.5970659) q[3];
sx q[3];
rz(-2.209112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2782044) q[2];
sx q[2];
rz(-0.74625838) q[2];
sx q[2];
rz(-2.9928652) q[2];
rz(-1.1778098) q[3];
sx q[3];
rz(-1.1921459) q[3];
sx q[3];
rz(1.2744133) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3026368) q[0];
sx q[0];
rz(-1.9479072) q[0];
sx q[0];
rz(-2.1394011) q[0];
rz(1.3305371) q[1];
sx q[1];
rz(-1.9621153) q[1];
sx q[1];
rz(2.5340714) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40815263) q[0];
sx q[0];
rz(-2.3436553) q[0];
sx q[0];
rz(-0.95592461) q[0];
rz(0.27670105) q[2];
sx q[2];
rz(-1.9637798) q[2];
sx q[2];
rz(-0.0056841141) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3866736) q[1];
sx q[1];
rz(-0.52820871) q[1];
sx q[1];
rz(2.9186503) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63849498) q[3];
sx q[3];
rz(-0.83052626) q[3];
sx q[3];
rz(2.3414224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9017631) q[2];
sx q[2];
rz(-0.87368691) q[2];
sx q[2];
rz(-0.98881161) q[2];
rz(-1.4766988) q[3];
sx q[3];
rz(-1.8035382) q[3];
sx q[3];
rz(0.57507676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3636417) q[0];
sx q[0];
rz(-2.2626484) q[0];
sx q[0];
rz(1.2203891) q[0];
rz(-1.3149423) q[1];
sx q[1];
rz(-1.5053791) q[1];
sx q[1];
rz(-3.0016518) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51516384) q[0];
sx q[0];
rz(-0.73472154) q[0];
sx q[0];
rz(2.5236041) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0331421) q[2];
sx q[2];
rz(-2.2058479) q[2];
sx q[2];
rz(2.6376574) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3158886) q[1];
sx q[1];
rz(-1.53535) q[1];
sx q[1];
rz(-1.1152382) q[1];
rz(-pi) q[2];
rz(-1.7470105) q[3];
sx q[3];
rz(-0.38910633) q[3];
sx q[3];
rz(-1.2639623) q[3];
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
rz(3.0042082) q[2];
rz(1.8047699) q[3];
sx q[3];
rz(-1.5788014) q[3];
sx q[3];
rz(-3.1409851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-2.7430275) q[0];
sx q[0];
rz(-1.7337357) q[0];
sx q[0];
rz(0.16125691) q[0];
rz(-2.297961) q[1];
sx q[1];
rz(-0.7470986) q[1];
sx q[1];
rz(1.3074494) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4754514) q[0];
sx q[0];
rz(-1.3031811) q[0];
sx q[0];
rz(0.78351905) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0867811) q[2];
sx q[2];
rz(-2.7325411) q[2];
sx q[2];
rz(1.3889927) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.73518411) q[1];
sx q[1];
rz(-0.98957634) q[1];
sx q[1];
rz(3.0343949) q[1];
rz(-pi) q[2];
rz(-2.3639115) q[3];
sx q[3];
rz(-2.2984004) q[3];
sx q[3];
rz(2.889159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.97189409) q[2];
sx q[2];
rz(-0.43115386) q[2];
sx q[2];
rz(-2.2141854) q[2];
rz(1.5291322) q[3];
sx q[3];
rz(-2.0404787) q[3];
sx q[3];
rz(0.35759887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2138432) q[0];
sx q[0];
rz(-2.1723211) q[0];
sx q[0];
rz(1.4056322) q[0];
rz(1.4145781) q[1];
sx q[1];
rz(-0.79726338) q[1];
sx q[1];
rz(-2.3419211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7223631) q[0];
sx q[0];
rz(-1.6327614) q[0];
sx q[0];
rz(-3.1298766) q[0];
rz(-pi) q[1];
rz(-0.30971618) q[2];
sx q[2];
rz(-1.6982007) q[2];
sx q[2];
rz(-2.9744444) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.0015687167) q[1];
sx q[1];
rz(-0.42823165) q[1];
sx q[1];
rz(0.93044859) q[1];
rz(2.7148083) q[3];
sx q[3];
rz(-0.43922654) q[3];
sx q[3];
rz(2.3595485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8863525) q[2];
sx q[2];
rz(-2.4031874) q[2];
sx q[2];
rz(0.8026455) q[2];
rz(0.32043996) q[3];
sx q[3];
rz(-1.8010062) q[3];
sx q[3];
rz(1.6360487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7771512) q[0];
sx q[0];
rz(-3.1138804) q[0];
sx q[0];
rz(3.1112772) q[0];
rz(-1.8346571) q[1];
sx q[1];
rz(-1.0825284) q[1];
sx q[1];
rz(-0.62987769) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9014466) q[0];
sx q[0];
rz(-1.9909262) q[0];
sx q[0];
rz(1.7769369) q[0];
rz(2.3223619) q[2];
sx q[2];
rz(-0.60556817) q[2];
sx q[2];
rz(1.2450964) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4640208) q[1];
sx q[1];
rz(-1.7742762) q[1];
sx q[1];
rz(-2.8427441) q[1];
x q[2];
rz(0.23655741) q[3];
sx q[3];
rz(-0.45522296) q[3];
sx q[3];
rz(-0.61456028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1097172) q[2];
sx q[2];
rz(-0.66706359) q[2];
sx q[2];
rz(1.9943705) q[2];
rz(-2.4540497) q[3];
sx q[3];
rz(-0.59194982) q[3];
sx q[3];
rz(-1.8182925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3654093) q[0];
sx q[0];
rz(-2.9054346) q[0];
sx q[0];
rz(0.48015204) q[0];
rz(0.40223739) q[1];
sx q[1];
rz(-1.6419342) q[1];
sx q[1];
rz(0.38280907) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3881809) q[0];
sx q[0];
rz(-2.7482902) q[0];
sx q[0];
rz(0.9939145) q[0];
rz(-pi) q[1];
rz(-1.2023964) q[2];
sx q[2];
rz(-0.3109667) q[2];
sx q[2];
rz(0.63705618) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.20732197) q[1];
sx q[1];
rz(-1.2756961) q[1];
sx q[1];
rz(-1.9701824) q[1];
rz(-2.945369) q[3];
sx q[3];
rz(-1.5999402) q[3];
sx q[3];
rz(-1.9118229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39168921) q[2];
sx q[2];
rz(-0.73035556) q[2];
sx q[2];
rz(0.35823092) q[2];
rz(2.7273438) q[3];
sx q[3];
rz(-1.0284938) q[3];
sx q[3];
rz(-0.39345583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2888032) q[0];
sx q[0];
rz(-1.2816592) q[0];
sx q[0];
rz(-2.7237256) q[0];
rz(1.6425543) q[1];
sx q[1];
rz(-0.42048979) q[1];
sx q[1];
rz(0.61161673) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2177416) q[0];
sx q[0];
rz(-1.799066) q[0];
sx q[0];
rz(1.3598674) q[0];
rz(0.20576818) q[2];
sx q[2];
rz(-1.5094583) q[2];
sx q[2];
rz(-2.105728) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.043165) q[1];
sx q[1];
rz(-1.0994776) q[1];
sx q[1];
rz(-0.31402513) q[1];
rz(-1.6526653) q[3];
sx q[3];
rz(-1.22504) q[3];
sx q[3];
rz(0.053178259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.186782) q[2];
sx q[2];
rz(-1.45767) q[2];
sx q[2];
rz(-2.7853454) q[2];
rz(-2.6567843) q[3];
sx q[3];
rz(-1.7907413) q[3];
sx q[3];
rz(-3.1118605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32271069) q[0];
sx q[0];
rz(-2.0368545) q[0];
sx q[0];
rz(1.4622965) q[0];
rz(-0.38743427) q[1];
sx q[1];
rz(-0.92715627) q[1];
sx q[1];
rz(-0.80518728) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0793251) q[0];
sx q[0];
rz(-2.1342417) q[0];
sx q[0];
rz(-0.84765537) q[0];
rz(-pi) q[1];
rz(-0.83909713) q[2];
sx q[2];
rz(-2.5658414) q[2];
sx q[2];
rz(-0.29345278) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2428875) q[1];
sx q[1];
rz(-1.4243898) q[1];
sx q[1];
rz(2.8553748) q[1];
rz(-pi) q[2];
x q[2];
rz(1.322454) q[3];
sx q[3];
rz(-1.2183918) q[3];
sx q[3];
rz(3.0737511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8055973) q[2];
sx q[2];
rz(-0.62312859) q[2];
sx q[2];
rz(-1.222329) q[2];
rz(2.3313816) q[3];
sx q[3];
rz(-0.92456341) q[3];
sx q[3];
rz(1.5626102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41351086) q[0];
sx q[0];
rz(-1.6376729) q[0];
sx q[0];
rz(-1.6117657) q[0];
rz(1.275508) q[1];
sx q[1];
rz(-2.7559912) q[1];
sx q[1];
rz(-1.4082946) q[1];
rz(1.8985844) q[2];
sx q[2];
rz(-2.7375883) q[2];
sx q[2];
rz(1.3069456) q[2];
rz(-2.2386407) q[3];
sx q[3];
rz(-1.9330238) q[3];
sx q[3];
rz(-0.15180363) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
