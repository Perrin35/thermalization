OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(1.8338058) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(4.0772822) q[1];
sx q[1];
rz(4.7128591) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21068621) q[0];
sx q[0];
rz(-1.2331729) q[0];
sx q[0];
rz(0.36436413) q[0];
rz(-pi) q[1];
rz(-0.70648944) q[2];
sx q[2];
rz(-2.2405365) q[2];
sx q[2];
rz(2.0073839) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5174487) q[1];
sx q[1];
rz(-1.9074719) q[1];
sx q[1];
rz(1.8271853) q[1];
rz(-2.9458463) q[3];
sx q[3];
rz(-1.2139075) q[3];
sx q[3];
rz(1.2008047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2661665) q[2];
sx q[2];
rz(-0.29310075) q[2];
sx q[2];
rz(1.1323294) q[2];
rz(1.6752361) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(-1.0124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19673008) q[0];
sx q[0];
rz(-2.9319627) q[0];
sx q[0];
rz(0.18584132) q[0];
rz(0.56022412) q[1];
sx q[1];
rz(-1.2954243) q[1];
sx q[1];
rz(2.9247608) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8502055) q[0];
sx q[0];
rz(-2.4225525) q[0];
sx q[0];
rz(2.015381) q[0];
x q[1];
rz(2.6685643) q[2];
sx q[2];
rz(-1.066726) q[2];
sx q[2];
rz(0.31993983) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5237907) q[1];
sx q[1];
rz(-2.4317867) q[1];
sx q[1];
rz(1.8989423) q[1];
rz(-0.87644491) q[3];
sx q[3];
rz(-2.090824) q[3];
sx q[3];
rz(2.2976573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.310114) q[2];
sx q[2];
rz(-2.3159413) q[2];
sx q[2];
rz(-1.8537834) q[2];
rz(0.76256049) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(-2.8365703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.6771616) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(-0.60423869) q[0];
rz(1.8151981) q[1];
sx q[1];
rz(-1.7809968) q[1];
sx q[1];
rz(0.93260971) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55721012) q[0];
sx q[0];
rz(-0.32214468) q[0];
sx q[0];
rz(-1.7412709) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7980174) q[2];
sx q[2];
rz(-1.7482687) q[2];
sx q[2];
rz(-0.15572671) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.31308094) q[1];
sx q[1];
rz(-0.57144895) q[1];
sx q[1];
rz(2.9213195) q[1];
rz(0.91498418) q[3];
sx q[3];
rz(-0.6647771) q[3];
sx q[3];
rz(1.2802326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.147826) q[2];
sx q[2];
rz(-1.0819165) q[2];
sx q[2];
rz(1.0926584) q[2];
rz(0.5422194) q[3];
sx q[3];
rz(-1.0850302) q[3];
sx q[3];
rz(0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3595235) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(-0.50022593) q[0];
rz(-2.3362828) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(-1.6436228) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7786176) q[0];
sx q[0];
rz(-0.59016363) q[0];
sx q[0];
rz(-2.1182563) q[0];
rz(1.285032) q[2];
sx q[2];
rz(-0.19837241) q[2];
sx q[2];
rz(2.6464268) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8151617) q[1];
sx q[1];
rz(-1.3109428) q[1];
sx q[1];
rz(-1.3304779) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46521503) q[3];
sx q[3];
rz(-1.1606996) q[3];
sx q[3];
rz(1.0614392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.74636373) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(-0.70181075) q[2];
rz(0.83135215) q[3];
sx q[3];
rz(-0.96389198) q[3];
sx q[3];
rz(-2.519616) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2410626) q[0];
sx q[0];
rz(-2.5456972) q[0];
sx q[0];
rz(-0.81533122) q[0];
rz(-1.5218081) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(-2.0933847) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5532903) q[0];
sx q[0];
rz(-2.7748845) q[0];
sx q[0];
rz(-2.328863) q[0];
rz(-pi) q[1];
x q[1];
rz(0.016096073) q[2];
sx q[2];
rz(-2.490009) q[2];
sx q[2];
rz(-0.96166699) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8909) q[1];
sx q[1];
rz(-2.7437468) q[1];
sx q[1];
rz(2.489151) q[1];
rz(-pi) q[2];
rz(-1.0190373) q[3];
sx q[3];
rz(-2.4487547) q[3];
sx q[3];
rz(2.2321731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6158225) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(1.099951) q[2];
rz(2.3163017) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(2.2560789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0734171) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(-2.2391879) q[0];
rz(1.0166608) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(-0.12983233) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3467305) q[0];
sx q[0];
rz(-0.71160337) q[0];
sx q[0];
rz(-2.5559588) q[0];
rz(-pi) q[1];
rz(-0.36826276) q[2];
sx q[2];
rz(-1.0266745) q[2];
sx q[2];
rz(-2.1899109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9784769) q[1];
sx q[1];
rz(-0.69677959) q[1];
sx q[1];
rz(0.58660581) q[1];
rz(2.5927605) q[3];
sx q[3];
rz(-1.4026814) q[3];
sx q[3];
rz(2.7041534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3123902) q[2];
sx q[2];
rz(-2.1924993) q[2];
sx q[2];
rz(-2.9373346) q[2];
rz(1.2060818) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(-2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7234574) q[0];
sx q[0];
rz(-1.3292987) q[0];
sx q[0];
rz(-1.6947421) q[0];
rz(1.2591259) q[1];
sx q[1];
rz(-2.1513758) q[1];
sx q[1];
rz(2.4553305) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26353729) q[0];
sx q[0];
rz(-1.1362846) q[0];
sx q[0];
rz(2.6146019) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1967728) q[2];
sx q[2];
rz(-0.84569028) q[2];
sx q[2];
rz(0.041989728) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4916617) q[1];
sx q[1];
rz(-1.4964536) q[1];
sx q[1];
rz(1.3851623) q[1];
rz(-pi) q[2];
rz(-1.8920184) q[3];
sx q[3];
rz(-1.1676844) q[3];
sx q[3];
rz(1.5954799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69616047) q[2];
sx q[2];
rz(-1.3779209) q[2];
sx q[2];
rz(0.0017722842) q[2];
rz(2.5799675) q[3];
sx q[3];
rz(-2.2300945) q[3];
sx q[3];
rz(1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5381662) q[0];
sx q[0];
rz(-2.4551233) q[0];
sx q[0];
rz(1.4461393) q[0];
rz(-2.360545) q[1];
sx q[1];
rz(-1.8361517) q[1];
sx q[1];
rz(-1.5015645) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7088889) q[0];
sx q[0];
rz(-0.49312691) q[0];
sx q[0];
rz(1.6119484) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5090452) q[2];
sx q[2];
rz(-1.5442863) q[2];
sx q[2];
rz(-0.74552958) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0280684) q[1];
sx q[1];
rz(-1.4816195) q[1];
sx q[1];
rz(-0.32797565) q[1];
rz(-pi) q[2];
rz(0.99174188) q[3];
sx q[3];
rz(-2.1175044) q[3];
sx q[3];
rz(0.21608298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.35187307) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(-1.3191351) q[2];
rz(1.9296648) q[3];
sx q[3];
rz(-1.8550248) q[3];
sx q[3];
rz(-2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8050352) q[0];
sx q[0];
rz(-2.5890077) q[0];
sx q[0];
rz(-1.9375027) q[0];
rz(-0.38326344) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(-2.7899172) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29031819) q[0];
sx q[0];
rz(-1.1997249) q[0];
sx q[0];
rz(-1.0822269) q[0];
rz(-pi) q[1];
rz(-2.4497689) q[2];
sx q[2];
rz(-1.8080538) q[2];
sx q[2];
rz(2.0123864) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55191509) q[1];
sx q[1];
rz(-2.3304686) q[1];
sx q[1];
rz(0.07304904) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6230691) q[3];
sx q[3];
rz(-1.4613323) q[3];
sx q[3];
rz(1.048552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.3433156) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(1.8593672) q[2];
rz(1.6451689) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(-2.0857247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6431817) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(0.19432755) q[0];
rz(-1.0378029) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(1.0338354) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68034222) q[0];
sx q[0];
rz(-2.470812) q[0];
sx q[0];
rz(1.781342) q[0];
rz(0.50844426) q[2];
sx q[2];
rz(-2.2061081) q[2];
sx q[2];
rz(-2.9490162) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0018113) q[1];
sx q[1];
rz(-1.0011295) q[1];
sx q[1];
rz(0.12057481) q[1];
x q[2];
rz(-1.7887573) q[3];
sx q[3];
rz(-2.2721604) q[3];
sx q[3];
rz(-2.4234114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0795435) q[2];
sx q[2];
rz(-2.1958308) q[2];
sx q[2];
rz(0.6357843) q[2];
rz(-2.87129) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(-1.5283782) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6939659) q[0];
sx q[0];
rz(-1.3128558) q[0];
sx q[0];
rz(-2.0679612) q[0];
rz(-1.4355961) q[1];
sx q[1];
rz(-1.5626848) q[1];
sx q[1];
rz(-2.3609153) q[1];
rz(-2.1736017) q[2];
sx q[2];
rz(-1.5445166) q[2];
sx q[2];
rz(-1.9894285) q[2];
rz(3.071143) q[3];
sx q[3];
rz(-1.2415213) q[3];
sx q[3];
rz(0.51125676) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];