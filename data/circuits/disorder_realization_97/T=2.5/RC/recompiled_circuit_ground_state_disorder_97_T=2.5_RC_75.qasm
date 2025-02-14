OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.41807732) q[0];
sx q[0];
rz(-0.014208566) q[0];
sx q[0];
rz(3.1007822) q[0];
rz(-0.040925097) q[1];
sx q[1];
rz(-1.0558145) q[1];
sx q[1];
rz(2.5850886) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7303527) q[0];
sx q[0];
rz(-1.7539193) q[0];
sx q[0];
rz(1.2191811) q[0];
rz(2.7950613) q[2];
sx q[2];
rz(-0.60014562) q[2];
sx q[2];
rz(0.66020667) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5750016) q[1];
sx q[1];
rz(-2.0739158) q[1];
sx q[1];
rz(1.9921705) q[1];
x q[2];
rz(-1.1993042) q[3];
sx q[3];
rz(-0.25361922) q[3];
sx q[3];
rz(-0.41244464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9875662) q[2];
sx q[2];
rz(-1.3904927) q[2];
sx q[2];
rz(1.7517368) q[2];
rz(-0.050358199) q[3];
sx q[3];
rz(-1.7213089) q[3];
sx q[3];
rz(-1.095298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20767009) q[0];
sx q[0];
rz(-1.727513) q[0];
sx q[0];
rz(3.0811144) q[0];
rz(-2.1531847) q[1];
sx q[1];
rz(-1.6840839) q[1];
sx q[1];
rz(-2.9284533) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6211301) q[0];
sx q[0];
rz(-2.222837) q[0];
sx q[0];
rz(-1.0494285) q[0];
rz(-pi) q[1];
rz(-0.95560851) q[2];
sx q[2];
rz(-2.9434574) q[2];
sx q[2];
rz(-0.83189925) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4276398) q[1];
sx q[1];
rz(-2.8866077) q[1];
sx q[1];
rz(2.990635) q[1];
x q[2];
rz(1.2422823) q[3];
sx q[3];
rz(-2.6217048) q[3];
sx q[3];
rz(-1.1240791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.1468143) q[2];
sx q[2];
rz(-1.549336) q[2];
sx q[2];
rz(3.1352622) q[2];
rz(-1.4396884) q[3];
sx q[3];
rz(-0.78841811) q[3];
sx q[3];
rz(-2.7009713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.45563662) q[0];
sx q[0];
rz(-0.39304471) q[0];
sx q[0];
rz(1.3156369) q[0];
rz(-1.7542959) q[1];
sx q[1];
rz(-2.5376814) q[1];
sx q[1];
rz(-0.046262892) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14213698) q[0];
sx q[0];
rz(-1.6007413) q[0];
sx q[0];
rz(1.5697451) q[0];
rz(-pi) q[1];
rz(-0.46628433) q[2];
sx q[2];
rz(-1.2215677) q[2];
sx q[2];
rz(0.081307129) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8475593) q[1];
sx q[1];
rz(-0.98638505) q[1];
sx q[1];
rz(-1.9895893) q[1];
x q[2];
rz(0.1667695) q[3];
sx q[3];
rz(-2.178949) q[3];
sx q[3];
rz(1.8812572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5837273) q[2];
sx q[2];
rz(-1.9670308) q[2];
sx q[2];
rz(-2.2306856) q[2];
rz(-1.2298498) q[3];
sx q[3];
rz(-0.76084796) q[3];
sx q[3];
rz(0.41874829) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2428899) q[0];
sx q[0];
rz(-1.7430584) q[0];
sx q[0];
rz(-2.718495) q[0];
rz(-2.2078919) q[1];
sx q[1];
rz(-1.9405148) q[1];
sx q[1];
rz(-0.43412128) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4431209) q[0];
sx q[0];
rz(-1.6054732) q[0];
sx q[0];
rz(0.29129819) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7291315) q[2];
sx q[2];
rz(-2.2248817) q[2];
sx q[2];
rz(-2.9323062) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0435183) q[1];
sx q[1];
rz(-1.6988108) q[1];
sx q[1];
rz(-0.56358595) q[1];
x q[2];
rz(-1.7755327) q[3];
sx q[3];
rz(-1.6833441) q[3];
sx q[3];
rz(-1.1074668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9551129) q[2];
sx q[2];
rz(-0.13804144) q[2];
sx q[2];
rz(2.1407927) q[2];
rz(-0.75438386) q[3];
sx q[3];
rz(-0.84027925) q[3];
sx q[3];
rz(1.3365356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56483018) q[0];
sx q[0];
rz(-0.63693988) q[0];
sx q[0];
rz(1.379396) q[0];
rz(-0.55343902) q[1];
sx q[1];
rz(-1.484) q[1];
sx q[1];
rz(2.4310506) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0934389) q[0];
sx q[0];
rz(-1.7293674) q[0];
sx q[0];
rz(2.3331654) q[0];
rz(-2.8997786) q[2];
sx q[2];
rz(-1.6266357) q[2];
sx q[2];
rz(-0.18130359) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2031999) q[1];
sx q[1];
rz(-1.5364931) q[1];
sx q[1];
rz(-1.8015979) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4031317) q[3];
sx q[3];
rz(-2.1049771) q[3];
sx q[3];
rz(-3.0287865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1826627) q[2];
sx q[2];
rz(-1.6561597) q[2];
sx q[2];
rz(0.043969285) q[2];
rz(0.26642695) q[3];
sx q[3];
rz(-3.0156342) q[3];
sx q[3];
rz(0.21336475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.63785315) q[0];
sx q[0];
rz(-1.6197236) q[0];
sx q[0];
rz(0.36852401) q[0];
rz(2.7875994) q[1];
sx q[1];
rz(-1.1292255) q[1];
sx q[1];
rz(-1.0634208) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2475488) q[0];
sx q[0];
rz(-3.0008033) q[0];
sx q[0];
rz(-1.6611499) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8485799) q[2];
sx q[2];
rz(-1.4627224) q[2];
sx q[2];
rz(1.9982823) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.271558) q[1];
sx q[1];
rz(-1.1599301) q[1];
sx q[1];
rz(-2.690587) q[1];
x q[2];
rz(1.3838718) q[3];
sx q[3];
rz(-0.83079544) q[3];
sx q[3];
rz(-0.79953098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4949558) q[2];
sx q[2];
rz(-1.8821303) q[2];
sx q[2];
rz(-1.4259526) q[2];
rz(-0.87092733) q[3];
sx q[3];
rz(-2.0723074) q[3];
sx q[3];
rz(-0.23437414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9694685) q[0];
sx q[0];
rz(-0.5640465) q[0];
sx q[0];
rz(-0.010757541) q[0];
rz(-2.5116008) q[1];
sx q[1];
rz(-2.0959581) q[1];
sx q[1];
rz(2.2363372) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9173316) q[0];
sx q[0];
rz(-1.0380099) q[0];
sx q[0];
rz(-3.0661746) q[0];
x q[1];
rz(0.33335061) q[2];
sx q[2];
rz(-1.0819737) q[2];
sx q[2];
rz(2.6437063) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.9295302) q[1];
sx q[1];
rz(-1.7133949) q[1];
sx q[1];
rz(2.6349956) q[1];
x q[2];
rz(2.1077044) q[3];
sx q[3];
rz(-2.4126588) q[3];
sx q[3];
rz(-1.1400853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4524727) q[2];
sx q[2];
rz(-1.8931188) q[2];
sx q[2];
rz(3.0242237) q[2];
rz(-0.72702414) q[3];
sx q[3];
rz(-2.2456462) q[3];
sx q[3];
rz(-2.0015543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82379502) q[0];
sx q[0];
rz(-1.057484) q[0];
sx q[0];
rz(-0.29528883) q[0];
rz(1.3914039) q[1];
sx q[1];
rz(-0.65381217) q[1];
sx q[1];
rz(0.48694912) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27010879) q[0];
sx q[0];
rz(-2.919533) q[0];
sx q[0];
rz(-0.47501774) q[0];
x q[1];
rz(-1.0435264) q[2];
sx q[2];
rz(-0.57948105) q[2];
sx q[2];
rz(-0.5652437) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.036186342) q[1];
sx q[1];
rz(-2.3008022) q[1];
sx q[1];
rz(-2.1814815) q[1];
x q[2];
rz(2.8744633) q[3];
sx q[3];
rz(-1.7448269) q[3];
sx q[3];
rz(1.4166299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9140909) q[2];
sx q[2];
rz(-2.5469683) q[2];
sx q[2];
rz(-1.5957069) q[2];
rz(0.60901904) q[3];
sx q[3];
rz(-2.0273429) q[3];
sx q[3];
rz(2.473089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85536999) q[0];
sx q[0];
rz(-0.90317059) q[0];
sx q[0];
rz(-1.5089996) q[0];
rz(1.1434932) q[1];
sx q[1];
rz(-1.9805464) q[1];
sx q[1];
rz(0.11018363) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71188739) q[0];
sx q[0];
rz(-1.4674057) q[0];
sx q[0];
rz(-2.5508419) q[0];
rz(2.0333336) q[2];
sx q[2];
rz(-1.8018844) q[2];
sx q[2];
rz(-3.0840735) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5563123) q[1];
sx q[1];
rz(-1.9686799) q[1];
sx q[1];
rz(2.3817987) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5694782) q[3];
sx q[3];
rz(-2.0681046) q[3];
sx q[3];
rz(-2.8921102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9785568) q[2];
sx q[2];
rz(-2.0608993) q[2];
sx q[2];
rz(2.2829096) q[2];
rz(-1.5052172) q[3];
sx q[3];
rz(-1.7965763) q[3];
sx q[3];
rz(1.4428562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7929512) q[0];
sx q[0];
rz(-2.0703147) q[0];
sx q[0];
rz(-0.68565482) q[0];
rz(0.81949702) q[1];
sx q[1];
rz(-1.4792496) q[1];
sx q[1];
rz(-2.562838) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9070918) q[0];
sx q[0];
rz(-1.0847738) q[0];
sx q[0];
rz(2.4810243) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0039601) q[2];
sx q[2];
rz(-1.3610148) q[2];
sx q[2];
rz(-2.6875241) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.47992771) q[1];
sx q[1];
rz(-2.6746422) q[1];
sx q[1];
rz(-0.60162094) q[1];
rz(-pi) q[2];
rz(-1.5220422) q[3];
sx q[3];
rz(-1.4617059) q[3];
sx q[3];
rz(-1.8439807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8389575) q[2];
sx q[2];
rz(-0.18582782) q[2];
sx q[2];
rz(-2.9616984) q[2];
rz(-2.5707865) q[3];
sx q[3];
rz(-0.72206086) q[3];
sx q[3];
rz(-0.65176898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8410692) q[0];
sx q[0];
rz(-2.4916334) q[0];
sx q[0];
rz(-1.4148225) q[0];
rz(-0.48814804) q[1];
sx q[1];
rz(-0.5562677) q[1];
sx q[1];
rz(1.5604896) q[1];
rz(1.0875779) q[2];
sx q[2];
rz(-0.69003006) q[2];
sx q[2];
rz(-0.48984169) q[2];
rz(0.53849738) q[3];
sx q[3];
rz(-1.8168728) q[3];
sx q[3];
rz(-1.3088165) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
