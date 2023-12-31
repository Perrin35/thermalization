OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0302304) q[0];
sx q[0];
rz(-0.65519873) q[0];
sx q[0];
rz(0.97638786) q[0];
rz(-0.916565) q[1];
sx q[1];
rz(2.38148) q[1];
sx q[1];
rz(10.217477) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9206031) q[0];
sx q[0];
rz(-1.594829) q[0];
sx q[0];
rz(2.241797) q[0];
x q[1];
rz(2.5976074) q[2];
sx q[2];
rz(-1.3731125) q[2];
sx q[2];
rz(-0.20180861) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1945038) q[1];
sx q[1];
rz(-2.6852486) q[1];
sx q[1];
rz(-1.4049005) q[1];
x q[2];
rz(-0.71346618) q[3];
sx q[3];
rz(-0.93868512) q[3];
sx q[3];
rz(1.5330832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.67316002) q[2];
sx q[2];
rz(-2.1940553) q[2];
sx q[2];
rz(-0.49896487) q[2];
rz(-2.7178102) q[3];
sx q[3];
rz(-1.9215923) q[3];
sx q[3];
rz(0.012714816) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38133165) q[0];
sx q[0];
rz(-1.6853764) q[0];
sx q[0];
rz(0.99200845) q[0];
rz(2.965062) q[1];
sx q[1];
rz(-0.20543988) q[1];
sx q[1];
rz(2.7761249) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9734946) q[0];
sx q[0];
rz(-2.0272278) q[0];
sx q[0];
rz(1.8512076) q[0];
rz(-pi) q[1];
rz(2.2720049) q[2];
sx q[2];
rz(-0.80539942) q[2];
sx q[2];
rz(-2.373901) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9531627) q[1];
sx q[1];
rz(-0.97524446) q[1];
sx q[1];
rz(-1.9838536) q[1];
rz(2.6037381) q[3];
sx q[3];
rz(-2.9085277) q[3];
sx q[3];
rz(2.2190998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3084597) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(-2.0642521) q[2];
rz(2.9789553) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(-1.8796273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4804374) q[0];
sx q[0];
rz(-2.2876331) q[0];
sx q[0];
rz(-2.5307632) q[0];
rz(-1.4340596) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(-2.6909713) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6621646) q[0];
sx q[0];
rz(-0.084708609) q[0];
sx q[0];
rz(-1.1685632) q[0];
x q[1];
rz(0.30655105) q[2];
sx q[2];
rz(-1.2129285) q[2];
sx q[2];
rz(-1.7265665) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.140481) q[1];
sx q[1];
rz(-1.2621659) q[1];
sx q[1];
rz(-0.66092296) q[1];
rz(-pi) q[2];
rz(-0.30626014) q[3];
sx q[3];
rz(-2.2420886) q[3];
sx q[3];
rz(-1.9426395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0314363) q[2];
sx q[2];
rz(-1.0855731) q[2];
sx q[2];
rz(-1.9499367) q[2];
rz(-2.2971161) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(0.58117956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095117005) q[0];
sx q[0];
rz(-1.3866383) q[0];
sx q[0];
rz(-2.0003831) q[0];
rz(1.8449239) q[1];
sx q[1];
rz(-2.7099687) q[1];
sx q[1];
rz(2.8362714) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0083864) q[0];
sx q[0];
rz(-1.2587233) q[0];
sx q[0];
rz(-0.58831373) q[0];
rz(-1.1496454) q[2];
sx q[2];
rz(-1.9990168) q[2];
sx q[2];
rz(-0.8840094) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0381283) q[1];
sx q[1];
rz(-2.3410019) q[1];
sx q[1];
rz(-0.14726463) q[1];
rz(-pi) q[2];
rz(-0.51586113) q[3];
sx q[3];
rz(-2.0277884) q[3];
sx q[3];
rz(0.36992321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2083464) q[2];
sx q[2];
rz(-1.5657921) q[2];
sx q[2];
rz(0.34710458) q[2];
rz(0.38337213) q[3];
sx q[3];
rz(-0.88630668) q[3];
sx q[3];
rz(-1.0679519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.7326996) q[0];
sx q[0];
rz(-0.86298958) q[0];
sx q[0];
rz(-2.6026759) q[0];
rz(-1.4096889) q[1];
sx q[1];
rz(-0.46202818) q[1];
sx q[1];
rz(-2.7358823) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021279871) q[0];
sx q[0];
rz(-0.73182523) q[0];
sx q[0];
rz(-0.41723199) q[0];
x q[1];
rz(1.97498) q[2];
sx q[2];
rz(-0.6037854) q[2];
sx q[2];
rz(-3.0886138) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.09771422) q[1];
sx q[1];
rz(-1.3188625) q[1];
sx q[1];
rz(1.21752) q[1];
rz(-pi) q[2];
rz(0.24280663) q[3];
sx q[3];
rz(-0.77951509) q[3];
sx q[3];
rz(-1.3047583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.47200176) q[2];
sx q[2];
rz(-2.3386023) q[2];
sx q[2];
rz(1.6873138) q[2];
rz(2.3029095) q[3];
sx q[3];
rz(-1.8554153) q[3];
sx q[3];
rz(-1.7680761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.8307761) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(-0.40711656) q[0];
rz(2.282417) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(1.3173332) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2989527) q[0];
sx q[0];
rz(-1.5484122) q[0];
sx q[0];
rz(2.1138666) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69465722) q[2];
sx q[2];
rz(-2.3431871) q[2];
sx q[2];
rz(2.3677804) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4827022) q[1];
sx q[1];
rz(-1.5679949) q[1];
sx q[1];
rz(-2.3059694) q[1];
rz(-pi) q[2];
x q[2];
rz(2.156593) q[3];
sx q[3];
rz(-1.9337774) q[3];
sx q[3];
rz(1.3604878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.09981) q[2];
sx q[2];
rz(-2.2618497) q[2];
sx q[2];
rz(2.2463634) q[2];
rz(-1.2081395) q[3];
sx q[3];
rz(-0.67966214) q[3];
sx q[3];
rz(-2.4334548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5773425) q[0];
sx q[0];
rz(-2.0944493) q[0];
sx q[0];
rz(2.0136925) q[0];
rz(-0.212542) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(1.2316661) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54264442) q[0];
sx q[0];
rz(-1.066474) q[0];
sx q[0];
rz(1.1619316) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2052843) q[2];
sx q[2];
rz(-1.0299183) q[2];
sx q[2];
rz(0.66334954) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5269446) q[1];
sx q[1];
rz(-1.6741747) q[1];
sx q[1];
rz(-2.9109216) q[1];
rz(-pi) q[2];
rz(0.96630649) q[3];
sx q[3];
rz(-1.5887727) q[3];
sx q[3];
rz(3.0997961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3125375) q[2];
sx q[2];
rz(-0.91784358) q[2];
sx q[2];
rz(-0.60116872) q[2];
rz(0.51856315) q[3];
sx q[3];
rz(-2.8278567) q[3];
sx q[3];
rz(-1.4275838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9940051) q[0];
sx q[0];
rz(-0.00031358263) q[0];
sx q[0];
rz(0.073154733) q[0];
rz(-2.1144497) q[1];
sx q[1];
rz(-1.5348397) q[1];
sx q[1];
rz(-2.5491319) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5151278) q[0];
sx q[0];
rz(-2.160762) q[0];
sx q[0];
rz(1.5149087) q[0];
x q[1];
rz(-2.4249435) q[2];
sx q[2];
rz(-0.71084329) q[2];
sx q[2];
rz(0.58472842) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8471624) q[1];
sx q[1];
rz(-1.2420328) q[1];
sx q[1];
rz(-2.2734103) q[1];
x q[2];
rz(1.9755649) q[3];
sx q[3];
rz(-1.3169857) q[3];
sx q[3];
rz(-0.39241957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56117326) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(0.44211659) q[2];
rz(-2.2413065) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1135547) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(1.8369209) q[0];
rz(-1.4470944) q[1];
sx q[1];
rz(-1.9686952) q[1];
sx q[1];
rz(-2.5794199) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5484555) q[0];
sx q[0];
rz(-0.89567417) q[0];
sx q[0];
rz(2.0329342) q[0];
rz(-pi) q[1];
rz(-3.047692) q[2];
sx q[2];
rz(-1.0419838) q[2];
sx q[2];
rz(1.4096066) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1162029) q[1];
sx q[1];
rz(-0.89239489) q[1];
sx q[1];
rz(-0.62979001) q[1];
rz(-pi) q[2];
rz(0.4060679) q[3];
sx q[3];
rz(-1.9417524) q[3];
sx q[3];
rz(0.0036811034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.28329453) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(2.1600058) q[2];
rz(-0.0043407241) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(0.59345746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3146064) q[0];
sx q[0];
rz(-1.9208603) q[0];
sx q[0];
rz(0.04709588) q[0];
rz(1.754952) q[1];
sx q[1];
rz(-1.1562693) q[1];
sx q[1];
rz(3.093739) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76059607) q[0];
sx q[0];
rz(-2.2254235) q[0];
sx q[0];
rz(3.1050443) q[0];
rz(-1.1282975) q[2];
sx q[2];
rz(-1.6167165) q[2];
sx q[2];
rz(-2.0890582) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.88356) q[1];
sx q[1];
rz(-2.0000534) q[1];
sx q[1];
rz(2.2116823) q[1];
rz(1.6258532) q[3];
sx q[3];
rz(-1.9729455) q[3];
sx q[3];
rz(2.1189342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4973267) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(-0.93635526) q[2];
rz(-0.80091536) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(2.463533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8235648) q[0];
sx q[0];
rz(-2.2021273) q[0];
sx q[0];
rz(0.20903023) q[0];
rz(2.6387852) q[1];
sx q[1];
rz(-1.3192758) q[1];
sx q[1];
rz(-2.1154138) q[1];
rz(1.5488831) q[2];
sx q[2];
rz(-2.1219325) q[2];
sx q[2];
rz(1.3757201) q[2];
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
