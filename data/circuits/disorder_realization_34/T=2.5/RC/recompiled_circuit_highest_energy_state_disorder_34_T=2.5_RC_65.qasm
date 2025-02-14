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
rz(1.8129803) q[0];
sx q[0];
rz(12.145481) q[0];
rz(2.4731877) q[1];
sx q[1];
rz(-1.547812) q[1];
sx q[1];
rz(-1.9903543) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4822455) q[0];
sx q[0];
rz(-1.7719685) q[0];
sx q[0];
rz(-2.6849999) q[0];
rz(-pi) q[1];
rz(-1.0082863) q[2];
sx q[2];
rz(-1.1739302) q[2];
sx q[2];
rz(-2.4728554) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0033773) q[1];
sx q[1];
rz(-2.1103519) q[1];
sx q[1];
rz(-0.92745499) q[1];
x q[2];
rz(-2.2231021) q[3];
sx q[3];
rz(-1.3555694) q[3];
sx q[3];
rz(1.8818242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.067817299) q[2];
sx q[2];
rz(-0.99738085) q[2];
sx q[2];
rz(1.0214405) q[2];
rz(1.3218309) q[3];
sx q[3];
rz(-2.6363711) q[3];
sx q[3];
rz(-0.57636133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3385056) q[0];
sx q[0];
rz(-2.1144688) q[0];
sx q[0];
rz(2.8402253) q[0];
rz(2.001568) q[1];
sx q[1];
rz(-2.3224484) q[1];
sx q[1];
rz(1.0637306) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50142215) q[0];
sx q[0];
rz(-2.8482264) q[0];
sx q[0];
rz(-1.1430995) q[0];
rz(0.21145169) q[2];
sx q[2];
rz(-2.3686948) q[2];
sx q[2];
rz(0.37040243) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1769907) q[1];
sx q[1];
rz(-2.1637205) q[1];
sx q[1];
rz(-1.9511392) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11874695) q[3];
sx q[3];
rz(-1.0380961) q[3];
sx q[3];
rz(1.1632196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.86338824) q[2];
sx q[2];
rz(-0.74625838) q[2];
sx q[2];
rz(0.14872742) q[2];
rz(1.9637828) q[3];
sx q[3];
rz(-1.9494467) q[3];
sx q[3];
rz(1.8671794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3026368) q[0];
sx q[0];
rz(-1.1936854) q[0];
sx q[0];
rz(2.1394011) q[0];
rz(-1.3305371) q[1];
sx q[1];
rz(-1.1794773) q[1];
sx q[1];
rz(-0.60752121) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70456767) q[0];
sx q[0];
rz(-1.9965197) q[0];
sx q[0];
rz(-2.2680437) q[0];
rz(-1.9776865) q[2];
sx q[2];
rz(-1.8259138) q[2];
sx q[2];
rz(-1.4681548) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62251112) q[1];
sx q[1];
rz(-1.4591328) q[1];
sx q[1];
rz(-0.51736781) q[1];
rz(0.99276279) q[3];
sx q[3];
rz(-2.2053763) q[3];
sx q[3];
rz(-1.5082593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2398296) q[2];
sx q[2];
rz(-2.2679057) q[2];
sx q[2];
rz(0.98881161) q[2];
rz(1.4766988) q[3];
sx q[3];
rz(-1.8035382) q[3];
sx q[3];
rz(-0.57507676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3636417) q[0];
sx q[0];
rz(-0.87894428) q[0];
sx q[0];
rz(1.2203891) q[0];
rz(-1.3149423) q[1];
sx q[1];
rz(-1.5053791) q[1];
sx q[1];
rz(0.13994089) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8927887) q[0];
sx q[0];
rz(-2.148845) q[0];
sx q[0];
rz(1.0885574) q[0];
rz(-0.60735382) q[2];
sx q[2];
rz(-0.80728856) q[2];
sx q[2];
rz(-1.2918778) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3158886) q[1];
sx q[1];
rz(-1.53535) q[1];
sx q[1];
rz(-2.0263544) q[1];
x q[2];
rz(-0.071752944) q[3];
sx q[3];
rz(-1.9535616) q[3];
sx q[3];
rz(-1.07384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6949029) q[2];
sx q[2];
rz(-2.0835154) q[2];
sx q[2];
rz(-3.0042082) q[2];
rz(-1.3368227) q[3];
sx q[3];
rz(-1.5788014) q[3];
sx q[3];
rz(-3.1409851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7430275) q[0];
sx q[0];
rz(-1.7337357) q[0];
sx q[0];
rz(0.16125691) q[0];
rz(-0.84363168) q[1];
sx q[1];
rz(-2.3944941) q[1];
sx q[1];
rz(-1.8341433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3529328) q[0];
sx q[0];
rz(-0.82214117) q[0];
sx q[0];
rz(-1.9400807) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9425386) q[2];
sx q[2];
rz(-1.2110333) q[2];
sx q[2];
rz(2.2729286) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4064085) q[1];
sx q[1];
rz(-0.98957634) q[1];
sx q[1];
rz(-3.0343949) q[1];
rz(0.90353538) q[3];
sx q[3];
rz(-2.1318815) q[3];
sx q[3];
rz(0.72431952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97189409) q[2];
sx q[2];
rz(-0.43115386) q[2];
sx q[2];
rz(2.2141854) q[2];
rz(1.5291322) q[3];
sx q[3];
rz(-1.1011139) q[3];
sx q[3];
rz(2.7839938) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9277495) q[0];
sx q[0];
rz(-0.96927154) q[0];
sx q[0];
rz(1.7359605) q[0];
rz(1.4145781) q[1];
sx q[1];
rz(-0.79726338) q[1];
sx q[1];
rz(0.79967156) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23223497) q[0];
sx q[0];
rz(-3.0785311) q[0];
sx q[0];
rz(-1.3841649) q[0];
rz(-1.437101) q[2];
sx q[2];
rz(-1.8779199) q[2];
sx q[2];
rz(-1.4442867) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1400239) q[1];
sx q[1];
rz(-0.42823165) q[1];
sx q[1];
rz(-2.2111441) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42678435) q[3];
sx q[3];
rz(-0.43922654) q[3];
sx q[3];
rz(2.3595485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.25524011) q[2];
sx q[2];
rz(-0.73840529) q[2];
sx q[2];
rz(0.8026455) q[2];
rz(0.32043996) q[3];
sx q[3];
rz(-1.3405864) q[3];
sx q[3];
rz(1.505544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7771512) q[0];
sx q[0];
rz(-0.027712263) q[0];
sx q[0];
rz(-3.1112772) q[0];
rz(-1.3069356) q[1];
sx q[1];
rz(-1.0825284) q[1];
sx q[1];
rz(0.62987769) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4277074) q[0];
sx q[0];
rz(-0.46525912) q[0];
sx q[0];
rz(-2.7119539) q[0];
x q[1];
rz(0.44158641) q[2];
sx q[2];
rz(-1.1418742) q[2];
sx q[2];
rz(-0.39583229) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4640208) q[1];
sx q[1];
rz(-1.7742762) q[1];
sx q[1];
rz(-0.29884853) q[1];
rz(-pi) q[2];
rz(-2.9050352) q[3];
sx q[3];
rz(-2.6863697) q[3];
sx q[3];
rz(-2.5270324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1097172) q[2];
sx q[2];
rz(-2.4745291) q[2];
sx q[2];
rz(-1.1472222) q[2];
rz(-0.68754292) q[3];
sx q[3];
rz(-2.5496428) q[3];
sx q[3];
rz(1.3233002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77618337) q[0];
sx q[0];
rz(-2.9054346) q[0];
sx q[0];
rz(0.48015204) q[0];
rz(2.7393553) q[1];
sx q[1];
rz(-1.4996585) q[1];
sx q[1];
rz(-2.7587836) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8653899) q[0];
sx q[0];
rz(-1.360219) q[0];
sx q[0];
rz(-1.2361069) q[0];
rz(-pi) q[1];
rz(1.2794957) q[2];
sx q[2];
rz(-1.460382) q[2];
sx q[2];
rz(1.2859273) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20732197) q[1];
sx q[1];
rz(-1.8658966) q[1];
sx q[1];
rz(1.1714102) q[1];
x q[2];
rz(-0.14842378) q[3];
sx q[3];
rz(-0.19834861) q[3];
sx q[3];
rz(2.655011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7499034) q[2];
sx q[2];
rz(-0.73035556) q[2];
sx q[2];
rz(0.35823092) q[2];
rz(-0.41424888) q[3];
sx q[3];
rz(-2.1130989) q[3];
sx q[3];
rz(-2.7481368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8527894) q[0];
sx q[0];
rz(-1.2816592) q[0];
sx q[0];
rz(-2.7237256) q[0];
rz(-1.4990384) q[1];
sx q[1];
rz(-2.7211029) q[1];
sx q[1];
rz(2.5299759) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8369524) q[0];
sx q[0];
rz(-1.3654183) q[0];
sx q[0];
rz(-2.9083328) q[0];
rz(2.8496004) q[2];
sx q[2];
rz(-0.21459178) q[2];
sx q[2];
rz(-2.8923182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4222718) q[1];
sx q[1];
rz(-0.55972717) q[1];
sx q[1];
rz(-1.0258963) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4889273) q[3];
sx q[3];
rz(-1.9165526) q[3];
sx q[3];
rz(-3.0884144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.186782) q[2];
sx q[2];
rz(-1.45767) q[2];
sx q[2];
rz(-2.7853454) q[2];
rz(0.48480836) q[3];
sx q[3];
rz(-1.3508513) q[3];
sx q[3];
rz(-0.029732186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.32271069) q[0];
sx q[0];
rz(-1.1047381) q[0];
sx q[0];
rz(1.4622965) q[0];
rz(0.38743427) q[1];
sx q[1];
rz(-0.92715627) q[1];
sx q[1];
rz(-2.3364054) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1772386) q[0];
sx q[0];
rz(-0.88429175) q[0];
sx q[0];
rz(-0.80857386) q[0];
rz(-pi) q[1];
rz(2.7323805) q[2];
sx q[2];
rz(-1.9878888) q[2];
sx q[2];
rz(0.52601782) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8565909) q[1];
sx q[1];
rz(-1.2877255) q[1];
sx q[1];
rz(1.4182752) q[1];
rz(-1.322454) q[3];
sx q[3];
rz(-1.9232009) q[3];
sx q[3];
rz(3.0737511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8055973) q[2];
sx q[2];
rz(-0.62312859) q[2];
sx q[2];
rz(1.222329) q[2];
rz(2.3313816) q[3];
sx q[3];
rz(-0.92456341) q[3];
sx q[3];
rz(-1.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(3.0048117) q[2];
sx q[2];
rz(-1.9521458) q[2];
sx q[2];
rz(-1.4804179) q[2];
rz(0.44966251) q[3];
sx q[3];
rz(-2.1884314) q[3];
sx q[3];
rz(-1.9951453) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
