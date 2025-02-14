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
rz(0.42088977) q[0];
rz(-0.66840494) q[1];
sx q[1];
rz(-1.5937807) q[1];
sx q[1];
rz(-1.1512383) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4822455) q[0];
sx q[0];
rz(-1.3696241) q[0];
sx q[0];
rz(2.6849999) q[0];
rz(-pi) q[1];
rz(-2.2368512) q[2];
sx q[2];
rz(-2.4657365) q[2];
sx q[2];
rz(-0.35206282) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1681663) q[1];
sx q[1];
rz(-0.81422537) q[1];
sx q[1];
rz(-0.78627118) q[1];
rz(1.2251159) q[3];
sx q[3];
rz(-0.68194032) q[3];
sx q[3];
rz(0.038393858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.067817299) q[2];
sx q[2];
rz(-0.99738085) q[2];
sx q[2];
rz(2.1201521) q[2];
rz(1.3218309) q[3];
sx q[3];
rz(-0.50522155) q[3];
sx q[3];
rz(-2.5652313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8030871) q[0];
sx q[0];
rz(-2.1144688) q[0];
sx q[0];
rz(0.3013674) q[0];
rz(-2.001568) q[1];
sx q[1];
rz(-2.3224484) q[1];
sx q[1];
rz(2.0778621) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057010896) q[0];
sx q[0];
rz(-1.304533) q[0];
sx q[0];
rz(0.12464704) q[0];
rz(-0.76164328) q[2];
sx q[2];
rz(-1.4237262) q[2];
sx q[2];
rz(-1.0479133) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96460198) q[1];
sx q[1];
rz(-2.1637205) q[1];
sx q[1];
rz(-1.1904535) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0349991) q[3];
sx q[3];
rz(-1.6730273) q[3];
sx q[3];
rz(0.4680948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2782044) q[2];
sx q[2];
rz(-0.74625838) q[2];
sx q[2];
rz(2.9928652) q[2];
rz(1.1778098) q[3];
sx q[3];
rz(-1.9494467) q[3];
sx q[3];
rz(-1.8671794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8389559) q[0];
sx q[0];
rz(-1.1936854) q[0];
sx q[0];
rz(-1.0021915) q[0];
rz(1.3305371) q[1];
sx q[1];
rz(-1.1794773) q[1];
sx q[1];
rz(-2.5340714) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.73344) q[0];
sx q[0];
rz(-2.3436553) q[0];
sx q[0];
rz(2.185668) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27670105) q[2];
sx q[2];
rz(-1.1778129) q[2];
sx q[2];
rz(0.0056841141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.62251112) q[1];
sx q[1];
rz(-1.6824598) q[1];
sx q[1];
rz(-0.51736781) q[1];
rz(0.72106169) q[3];
sx q[3];
rz(-2.0264056) q[3];
sx q[3];
rz(2.8351015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2398296) q[2];
sx q[2];
rz(-2.2679057) q[2];
sx q[2];
rz(2.152781) q[2];
rz(-1.6648939) q[3];
sx q[3];
rz(-1.8035382) q[3];
sx q[3];
rz(2.5665159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7779509) q[0];
sx q[0];
rz(-2.2626484) q[0];
sx q[0];
rz(1.9212035) q[0];
rz(1.3149423) q[1];
sx q[1];
rz(-1.5053791) q[1];
sx q[1];
rz(3.0016518) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8927887) q[0];
sx q[0];
rz(-0.99274764) q[0];
sx q[0];
rz(-1.0885574) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5342388) q[2];
sx q[2];
rz(-0.80728856) q[2];
sx q[2];
rz(-1.8497149) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3158886) q[1];
sx q[1];
rz(-1.6062427) q[1];
sx q[1];
rz(-2.0263544) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9544551) q[3];
sx q[3];
rz(-1.5042437) q[3];
sx q[3];
rz(2.6714747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4466897) q[2];
sx q[2];
rz(-2.0835154) q[2];
sx q[2];
rz(0.13738446) q[2];
rz(-1.8047699) q[3];
sx q[3];
rz(-1.5627912) q[3];
sx q[3];
rz(-3.1409851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3985652) q[0];
sx q[0];
rz(-1.7337357) q[0];
sx q[0];
rz(-2.9803357) q[0];
rz(2.297961) q[1];
sx q[1];
rz(-0.7470986) q[1];
sx q[1];
rz(1.8341433) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7886598) q[0];
sx q[0];
rz(-2.3194515) q[0];
sx q[0];
rz(-1.9400807) q[0];
rz(1.0867811) q[2];
sx q[2];
rz(-2.7325411) q[2];
sx q[2];
rz(-1.7526) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4064085) q[1];
sx q[1];
rz(-2.1520163) q[1];
sx q[1];
rz(-3.0343949) q[1];
rz(-pi) q[2];
rz(-0.67477711) q[3];
sx q[3];
rz(-2.1222847) q[3];
sx q[3];
rz(1.2434208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.97189409) q[2];
sx q[2];
rz(-0.43115386) q[2];
sx q[2];
rz(-0.92740721) q[2];
rz(-1.5291322) q[3];
sx q[3];
rz(-1.1011139) q[3];
sx q[3];
rz(0.35759887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.41922955) q[0];
sx q[0];
rz(-1.5088313) q[0];
sx q[0];
rz(3.1298766) q[0];
rz(-pi) q[1];
x q[1];
rz(1.437101) q[2];
sx q[2];
rz(-1.2636728) q[2];
sx q[2];
rz(-1.4442867) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4569163) q[1];
sx q[1];
rz(-1.9102736) q[1];
sx q[1];
rz(0.26626069) q[1];
rz(-pi) q[2];
rz(2.7148083) q[3];
sx q[3];
rz(-2.7023661) q[3];
sx q[3];
rz(0.78204417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25524011) q[2];
sx q[2];
rz(-0.73840529) q[2];
sx q[2];
rz(0.8026455) q[2];
rz(-0.32043996) q[3];
sx q[3];
rz(-1.3405864) q[3];
sx q[3];
rz(-1.505544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7771512) q[0];
sx q[0];
rz(-0.027712263) q[0];
sx q[0];
rz(3.1112772) q[0];
rz(1.3069356) q[1];
sx q[1];
rz(-2.0590643) q[1];
sx q[1];
rz(0.62987769) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8960285) q[0];
sx q[0];
rz(-1.3828074) q[0];
sx q[0];
rz(2.7134368) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44158641) q[2];
sx q[2];
rz(-1.9997184) q[2];
sx q[2];
rz(2.7457604) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4640208) q[1];
sx q[1];
rz(-1.3673165) q[1];
sx q[1];
rz(0.29884853) q[1];
x q[2];
rz(-0.44416645) q[3];
sx q[3];
rz(-1.6740177) q[3];
sx q[3];
rz(0.74300649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.031875413) q[2];
sx q[2];
rz(-0.66706359) q[2];
sx q[2];
rz(1.1472222) q[2];
rz(0.68754292) q[3];
sx q[3];
rz(-2.5496428) q[3];
sx q[3];
rz(1.8182925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77618337) q[0];
sx q[0];
rz(-0.23615806) q[0];
sx q[0];
rz(2.6614406) q[0];
rz(-2.7393553) q[1];
sx q[1];
rz(-1.6419342) q[1];
sx q[1];
rz(-2.7587836) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8653899) q[0];
sx q[0];
rz(-1.360219) q[0];
sx q[0];
rz(1.2361069) q[0];
rz(-3.026364) q[2];
sx q[2];
rz(-1.8602716) q[2];
sx q[2];
rz(-2.8897499) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.20732197) q[1];
sx q[1];
rz(-1.8658966) q[1];
sx q[1];
rz(-1.1714102) q[1];
rz(1.6005101) q[3];
sx q[3];
rz(-1.374657) q[3];
sx q[3];
rz(2.8063584) q[3];
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
rz(2.7481368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8527894) q[0];
sx q[0];
rz(-1.2816592) q[0];
sx q[0];
rz(-0.417867) q[0];
rz(-1.4990384) q[1];
sx q[1];
rz(-0.42048979) q[1];
sx q[1];
rz(-2.5299759) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9238511) q[0];
sx q[0];
rz(-1.799066) q[0];
sx q[0];
rz(1.7817253) q[0];
rz(-pi) q[1];
rz(-0.29199227) q[2];
sx q[2];
rz(-2.9270009) q[2];
sx q[2];
rz(-0.24927441) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.043165) q[1];
sx q[1];
rz(-1.0994776) q[1];
sx q[1];
rz(2.8275675) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2232333) q[3];
sx q[3];
rz(-2.7866552) q[3];
sx q[3];
rz(-0.29071537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9548107) q[2];
sx q[2];
rz(-1.6839226) q[2];
sx q[2];
rz(2.7853454) q[2];
rz(2.6567843) q[3];
sx q[3];
rz(-1.7907413) q[3];
sx q[3];
rz(3.1118605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32271069) q[0];
sx q[0];
rz(-1.1047381) q[0];
sx q[0];
rz(1.4622965) q[0];
rz(-2.7541584) q[1];
sx q[1];
rz(-0.92715627) q[1];
sx q[1];
rz(-2.3364054) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1925404) q[0];
sx q[0];
rz(-2.1645191) q[0];
sx q[0];
rz(-2.441382) q[0];
rz(-0.83909713) q[2];
sx q[2];
rz(-2.5658414) q[2];
sx q[2];
rz(-0.29345278) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.78811121) q[1];
sx q[1];
rz(-2.821021) q[1];
sx q[1];
rz(-0.48133565) q[1];
rz(-pi) q[2];
rz(2.7789742) q[3];
sx q[3];
rz(-1.3380074) q[3];
sx q[3];
rz(-1.5902558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.33599535) q[2];
sx q[2];
rz(-2.5184641) q[2];
sx q[2];
rz(-1.9192637) q[2];
rz(0.81021106) q[3];
sx q[3];
rz(-0.92456341) q[3];
sx q[3];
rz(1.5789825) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41351086) q[0];
sx q[0];
rz(-1.5039197) q[0];
sx q[0];
rz(1.529827) q[0];
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
