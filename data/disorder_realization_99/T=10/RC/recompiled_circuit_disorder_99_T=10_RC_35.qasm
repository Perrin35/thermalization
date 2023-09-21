OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1113623) q[0];
sx q[0];
rz(-2.4863939) q[0];
sx q[0];
rz(2.1652048) q[0];
rz(2.2250277) q[1];
sx q[1];
rz(-2.38148) q[1];
sx q[1];
rz(-2.3488933) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7615258) q[0];
sx q[0];
rz(-2.4702284) q[0];
sx q[0];
rz(1.5321561) q[0];
x q[1];
rz(1.8007457) q[2];
sx q[2];
rz(-2.1030428) q[2];
sx q[2];
rz(-1.4872273) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.76268643) q[1];
sx q[1];
rz(-1.1211809) q[1];
sx q[1];
rz(3.0607037) q[1];
rz(-pi) q[2];
rz(-0.71346618) q[3];
sx q[3];
rz(-0.93868512) q[3];
sx q[3];
rz(-1.6085094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67316002) q[2];
sx q[2];
rz(-0.94753733) q[2];
sx q[2];
rz(-2.6426278) q[2];
rz(-2.7178102) q[3];
sx q[3];
rz(-1.9215923) q[3];
sx q[3];
rz(-3.1288778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.38133165) q[0];
sx q[0];
rz(-1.4562162) q[0];
sx q[0];
rz(0.99200845) q[0];
rz(-0.17653067) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(-2.7761249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7468837) q[0];
sx q[0];
rz(-2.611126) q[0];
sx q[0];
rz(-2.6283426) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2720049) q[2];
sx q[2];
rz(-0.80539942) q[2];
sx q[2];
rz(2.373901) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5181527) q[1];
sx q[1];
rz(-1.232051) q[1];
sx q[1];
rz(2.5046196) q[1];
x q[2];
rz(-1.69181) q[3];
sx q[3];
rz(-1.3711208) q[3];
sx q[3];
rz(1.6691085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3084597) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(-2.0642521) q[2];
rz(-0.16263738) q[3];
sx q[3];
rz(-1.8789623) q[3];
sx q[3];
rz(1.8796273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66115528) q[0];
sx q[0];
rz(-2.2876331) q[0];
sx q[0];
rz(2.5307632) q[0];
rz(-1.4340596) q[1];
sx q[1];
rz(-1.3620946) q[1];
sx q[1];
rz(-0.45062137) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075899344) q[0];
sx q[0];
rz(-1.6487299) q[0];
sx q[0];
rz(0.033228544) q[0];
rz(-pi) q[1];
rz(1.1969823) q[2];
sx q[2];
rz(-1.8573559) q[2];
sx q[2];
rz(-2.8754004) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0011117) q[1];
sx q[1];
rz(-1.8794267) q[1];
sx q[1];
rz(2.4806697) q[1];
rz(-pi) q[2];
rz(2.8353325) q[3];
sx q[3];
rz(-2.2420886) q[3];
sx q[3];
rz(1.1989532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0314363) q[2];
sx q[2];
rz(-1.0855731) q[2];
sx q[2];
rz(1.9499367) q[2];
rz(0.84447652) q[3];
sx q[3];
rz(-0.28544393) q[3];
sx q[3];
rz(2.5604131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0464756) q[0];
sx q[0];
rz(-1.7549544) q[0];
sx q[0];
rz(1.1412096) q[0];
rz(-1.8449239) q[1];
sx q[1];
rz(-0.43162391) q[1];
sx q[1];
rz(-0.30532125) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13320623) q[0];
sx q[0];
rz(-1.8828694) q[0];
sx q[0];
rz(-2.5532789) q[0];
rz(-pi) q[1];
rz(-1.9919473) q[2];
sx q[2];
rz(-1.1425758) q[2];
sx q[2];
rz(-0.8840094) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.42974737) q[1];
sx q[1];
rz(-1.6763121) q[1];
sx q[1];
rz(-0.79515102) q[1];
x q[2];
rz(-1.0563072) q[3];
sx q[3];
rz(-2.0293651) q[3];
sx q[3];
rz(-1.4460627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9332463) q[2];
sx q[2];
rz(-1.5758005) q[2];
sx q[2];
rz(0.34710458) q[2];
rz(0.38337213) q[3];
sx q[3];
rz(-2.255286) q[3];
sx q[3];
rz(-2.0736407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(0.40889302) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(2.6026759) q[0];
rz(1.7319038) q[1];
sx q[1];
rz(-0.46202818) q[1];
sx q[1];
rz(0.40571037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1203128) q[0];
sx q[0];
rz(-0.73182523) q[0];
sx q[0];
rz(-2.7243607) q[0];
x q[1];
rz(-2.876725) q[2];
sx q[2];
rz(-2.1200074) q[2];
sx q[2];
rz(2.7153646) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0438784) q[1];
sx q[1];
rz(-1.3188625) q[1];
sx q[1];
rz(-1.9240727) q[1];
x q[2];
rz(2.898786) q[3];
sx q[3];
rz(-2.3620776) q[3];
sx q[3];
rz(-1.3047583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47200176) q[2];
sx q[2];
rz(-0.80299032) q[2];
sx q[2];
rz(1.6873138) q[2];
rz(-0.83868319) q[3];
sx q[3];
rz(-1.8554153) q[3];
sx q[3];
rz(1.3735166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31081653) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(-2.7344761) q[0];
rz(-0.85917568) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(-1.8242594) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8562397) q[0];
sx q[0];
rz(-1.0278773) q[0];
sx q[0];
rz(0.026144233) q[0];
rz(-2.4469354) q[2];
sx q[2];
rz(-0.79840556) q[2];
sx q[2];
rz(2.3677804) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4827022) q[1];
sx q[1];
rz(-1.5735978) q[1];
sx q[1];
rz(-2.3059694) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96885724) q[3];
sx q[3];
rz(-0.67776206) q[3];
sx q[3];
rz(2.4399151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.09981) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(-2.2463634) q[2];
rz(1.9334531) q[3];
sx q[3];
rz(-0.67966214) q[3];
sx q[3];
rz(-2.4334548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5773425) q[0];
sx q[0];
rz(-2.0944493) q[0];
sx q[0];
rz(-1.1279001) q[0];
rz(0.212542) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(1.9099265) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8679778) q[0];
sx q[0];
rz(-0.6379188) q[0];
sx q[0];
rz(-0.62423737) q[0];
rz(-pi) q[1];
rz(0.77881323) q[2];
sx q[2];
rz(-2.3328569) q[2];
sx q[2];
rz(2.8445809) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1612138) q[1];
sx q[1];
rz(-1.3413789) q[1];
sx q[1];
rz(1.6769665) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6024186) q[3];
sx q[3];
rz(-2.5368689) q[3];
sx q[3];
rz(-1.5550169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.82905519) q[2];
sx q[2];
rz(-2.2237491) q[2];
sx q[2];
rz(0.60116872) q[2];
rz(2.6230295) q[3];
sx q[3];
rz(-2.8278567) q[3];
sx q[3];
rz(-1.7140088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9940051) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(0.073154733) q[0];
rz(-2.1144497) q[1];
sx q[1];
rz(-1.606753) q[1];
sx q[1];
rz(-0.59246078) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97544599) q[0];
sx q[0];
rz(-1.5243634) q[0];
sx q[0];
rz(-2.5509044) q[0];
x q[1];
rz(-2.5657797) q[2];
sx q[2];
rz(-1.1278707) q[2];
sx q[2];
rz(-1.5695614) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2295099) q[1];
sx q[1];
rz(-0.76369897) q[1];
sx q[1];
rz(2.0565226) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.153272) q[3];
sx q[3];
rz(-0.47401014) q[3];
sx q[3];
rz(1.43309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5804194) q[2];
sx q[2];
rz(-1.4658804) q[2];
sx q[2];
rz(-2.6994761) q[2];
rz(-0.90028611) q[3];
sx q[3];
rz(-1.0628275) q[3];
sx q[3];
rz(3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.028037926) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(1.3046718) q[0];
rz(-1.6944983) q[1];
sx q[1];
rz(-1.1728975) q[1];
sx q[1];
rz(-2.5794199) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5484555) q[0];
sx q[0];
rz(-2.2459185) q[0];
sx q[0];
rz(2.0329342) q[0];
rz(-pi) q[1];
rz(1.7299037) q[2];
sx q[2];
rz(-2.6052887) q[2];
sx q[2];
rz(1.2250587) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8852946) q[1];
sx q[1];
rz(-2.2513299) q[1];
sx q[1];
rz(0.93974944) q[1];
rz(-pi) q[2];
rz(-2.3638944) q[3];
sx q[3];
rz(-2.5986528) q[3];
sx q[3];
rz(-2.267595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8582981) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(2.1600058) q[2];
rz(3.1372519) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(0.59345746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.8269862) q[0];
sx q[0];
rz(-1.2207323) q[0];
sx q[0];
rz(3.0944968) q[0];
rz(-1.3866407) q[1];
sx q[1];
rz(-1.1562693) q[1];
sx q[1];
rz(3.093739) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78794107) q[0];
sx q[0];
rz(-1.5418058) q[0];
sx q[0];
rz(2.2257462) q[0];
x q[1];
rz(1.1282975) q[2];
sx q[2];
rz(-1.6167165) q[2];
sx q[2];
rz(2.0890582) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0117482) q[1];
sx q[1];
rz(-0.99600345) q[1];
sx q[1];
rz(0.5188491) q[1];
x q[2];
rz(0.1286653) q[3];
sx q[3];
rz(-0.4056969) q[3];
sx q[3];
rz(2.2588244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4973267) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(-0.93635526) q[2];
rz(-0.80091536) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(-0.67805964) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8235648) q[0];
sx q[0];
rz(-2.2021273) q[0];
sx q[0];
rz(0.20903023) q[0];
rz(-2.6387852) q[1];
sx q[1];
rz(-1.8223169) q[1];
sx q[1];
rz(1.0261789) q[1];
rz(3.1059601) q[2];
sx q[2];
rz(-0.55152668) q[2];
sx q[2];
rz(1.3338911) q[2];
rz(1.1669284) q[3];
sx q[3];
rz(-1.231791) q[3];
sx q[3];
rz(0.024699208) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
