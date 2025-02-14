OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.06935057) q[0];
sx q[0];
rz(-2.0877593) q[0];
sx q[0];
rz(1.8928438) q[0];
rz(-1.2381923) q[1];
sx q[1];
rz(3.5817322) q[1];
sx q[1];
rz(11.323827) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1093501) q[0];
sx q[0];
rz(-2.1324131) q[0];
sx q[0];
rz(-1.1121145) q[0];
x q[1];
rz(-2.3435107) q[2];
sx q[2];
rz(-1.5028364) q[2];
sx q[2];
rz(-0.77347212) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.203071) q[1];
sx q[1];
rz(-1.8932098) q[1];
sx q[1];
rz(-1.5702269) q[1];
rz(-1.9162769) q[3];
sx q[3];
rz(-2.1519147) q[3];
sx q[3];
rz(-1.2863359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.093546346) q[2];
sx q[2];
rz(-2.8600433) q[2];
sx q[2];
rz(-1.2151037) q[2];
rz(1.9563227) q[3];
sx q[3];
rz(-0.90648854) q[3];
sx q[3];
rz(-1.5989446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1854061) q[0];
sx q[0];
rz(-2.7423999) q[0];
sx q[0];
rz(-1.6831552) q[0];
rz(-0.1419119) q[1];
sx q[1];
rz(-1.9344354) q[1];
sx q[1];
rz(1.9292319) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6120554) q[0];
sx q[0];
rz(-2.3891695) q[0];
sx q[0];
rz(1.0616395) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1375863) q[2];
sx q[2];
rz(-2.5494583) q[2];
sx q[2];
rz(-0.59321813) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.77245678) q[1];
sx q[1];
rz(-1.4765655) q[1];
sx q[1];
rz(2.2460031) q[1];
rz(-pi) q[2];
rz(-2.66414) q[3];
sx q[3];
rz(-2.4702854) q[3];
sx q[3];
rz(2.6539913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.096574664) q[2];
sx q[2];
rz(-0.18288945) q[2];
sx q[2];
rz(2.1796687) q[2];
rz(2.9436881) q[3];
sx q[3];
rz(-1.8353381) q[3];
sx q[3];
rz(-1.4977247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69979954) q[0];
sx q[0];
rz(-0.72023359) q[0];
sx q[0];
rz(1.066347) q[0];
rz(2.9329246) q[1];
sx q[1];
rz(-2.6228948) q[1];
sx q[1];
rz(0.64812237) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84635669) q[0];
sx q[0];
rz(-1.5180249) q[0];
sx q[0];
rz(-2.1480302) q[0];
rz(2.6510973) q[2];
sx q[2];
rz(-1.5083665) q[2];
sx q[2];
rz(1.708781) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3274182) q[1];
sx q[1];
rz(-1.2035554) q[1];
sx q[1];
rz(0.7805853) q[1];
rz(2.4345458) q[3];
sx q[3];
rz(-1.505369) q[3];
sx q[3];
rz(2.9312627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2842747) q[2];
sx q[2];
rz(-1.6782574) q[2];
sx q[2];
rz(0.54764444) q[2];
rz(-1.0037496) q[3];
sx q[3];
rz(-0.32310969) q[3];
sx q[3];
rz(0.16166648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6169287) q[0];
sx q[0];
rz(-1.1268317) q[0];
sx q[0];
rz(-2.774985) q[0];
rz(1.8560575) q[1];
sx q[1];
rz(-1.6211685) q[1];
sx q[1];
rz(2.3191648) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14299195) q[0];
sx q[0];
rz(-1.4944153) q[0];
sx q[0];
rz(0.19799302) q[0];
rz(2.9738178) q[2];
sx q[2];
rz(-1.3047403) q[2];
sx q[2];
rz(-0.38420907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.82717035) q[1];
sx q[1];
rz(-2.2232375) q[1];
sx q[1];
rz(-0.62364044) q[1];
rz(-pi) q[2];
rz(0.27536686) q[3];
sx q[3];
rz(-1.5726701) q[3];
sx q[3];
rz(-0.13260674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.503868) q[2];
sx q[2];
rz(-2.80105) q[2];
sx q[2];
rz(1.5558422) q[2];
rz(1.2268892) q[3];
sx q[3];
rz(-0.63819686) q[3];
sx q[3];
rz(-1.816412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1894839) q[0];
sx q[0];
rz(-1.1186849) q[0];
sx q[0];
rz(2.191191) q[0];
rz(2.9474126) q[1];
sx q[1];
rz(-1.8870528) q[1];
sx q[1];
rz(0.28087428) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0786983) q[0];
sx q[0];
rz(-0.67606976) q[0];
sx q[0];
rz(1.990114) q[0];
rz(2.9921586) q[2];
sx q[2];
rz(-1.583287) q[2];
sx q[2];
rz(-2.6246659) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.54714006) q[1];
sx q[1];
rz(-3.0277589) q[1];
sx q[1];
rz(0.74233858) q[1];
rz(-1.0626385) q[3];
sx q[3];
rz(-1.302037) q[3];
sx q[3];
rz(2.8063959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.48996553) q[2];
sx q[2];
rz(-1.4474892) q[2];
sx q[2];
rz(-2.4408686) q[2];
rz(1.0939595) q[3];
sx q[3];
rz(-2.14812) q[3];
sx q[3];
rz(0.82205621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99264282) q[0];
sx q[0];
rz(-2.094291) q[0];
sx q[0];
rz(0.5067504) q[0];
rz(-1.1243593) q[1];
sx q[1];
rz(-1.1686814) q[1];
sx q[1];
rz(-0.36875025) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71218894) q[0];
sx q[0];
rz(-0.87596873) q[0];
sx q[0];
rz(0.24526986) q[0];
x q[1];
rz(2.7733621) q[2];
sx q[2];
rz(-0.86129649) q[2];
sx q[2];
rz(1.8031098) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.74417392) q[1];
sx q[1];
rz(-1.4793921) q[1];
sx q[1];
rz(-1.869702) q[1];
rz(-2.7140824) q[3];
sx q[3];
rz(-1.9906133) q[3];
sx q[3];
rz(-1.9577718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90699482) q[2];
sx q[2];
rz(-0.58595053) q[2];
sx q[2];
rz(2.9242945) q[2];
rz(-1.6654738) q[3];
sx q[3];
rz(-0.81513351) q[3];
sx q[3];
rz(-0.34172094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9076964) q[0];
sx q[0];
rz(-2.8192769) q[0];
sx q[0];
rz(-2.3436558) q[0];
rz(-1.7012874) q[1];
sx q[1];
rz(-2.1013575) q[1];
sx q[1];
rz(0.61680102) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42617048) q[0];
sx q[0];
rz(-2.2035193) q[0];
sx q[0];
rz(1.6646845) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5011214) q[2];
sx q[2];
rz(-2.5129112) q[2];
sx q[2];
rz(1.6887661) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4492256) q[1];
sx q[1];
rz(-1.4265851) q[1];
sx q[1];
rz(-1.3647791) q[1];
x q[2];
rz(-2.439626) q[3];
sx q[3];
rz(-2.1118374) q[3];
sx q[3];
rz(-0.64149414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1282318) q[2];
sx q[2];
rz(-1.9888473) q[2];
sx q[2];
rz(-0.70179233) q[2];
rz(-0.93891406) q[3];
sx q[3];
rz(-2.2434668) q[3];
sx q[3];
rz(1.2452589) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.668648) q[0];
sx q[0];
rz(-0.15707792) q[0];
sx q[0];
rz(-3.0182086) q[0];
rz(-0.75048796) q[1];
sx q[1];
rz(-1.6672983) q[1];
sx q[1];
rz(-1.0184681) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3833419) q[0];
sx q[0];
rz(-1.4660379) q[0];
sx q[0];
rz(1.3488171) q[0];
rz(-2.0872714) q[2];
sx q[2];
rz(-1.9811842) q[2];
sx q[2];
rz(0.058319969) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2606517) q[1];
sx q[1];
rz(-0.48704942) q[1];
sx q[1];
rz(-3.0609291) q[1];
rz(-1.6945778) q[3];
sx q[3];
rz(-0.63237337) q[3];
sx q[3];
rz(0.06452175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5156775) q[2];
sx q[2];
rz(-1.6720142) q[2];
sx q[2];
rz(-0.53543004) q[2];
rz(1.911602) q[3];
sx q[3];
rz(-1.001469) q[3];
sx q[3];
rz(3.1297704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54315058) q[0];
sx q[0];
rz(-2.4813528) q[0];
sx q[0];
rz(-1.862233) q[0];
rz(-1.2641501) q[1];
sx q[1];
rz(-1.2707571) q[1];
sx q[1];
rz(0.15170161) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3705419) q[0];
sx q[0];
rz(-0.40239516) q[0];
sx q[0];
rz(-2.9950525) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13387605) q[2];
sx q[2];
rz(-1.5907484) q[2];
sx q[2];
rz(-2.6559115) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7599185) q[1];
sx q[1];
rz(-1.8823138) q[1];
sx q[1];
rz(0.98656922) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0284729) q[3];
sx q[3];
rz(-2.2185225) q[3];
sx q[3];
rz(-0.59779378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.58574) q[2];
sx q[2];
rz(-0.93583217) q[2];
sx q[2];
rz(-1.9110511) q[2];
rz(-1.3452283) q[3];
sx q[3];
rz(-1.3627005) q[3];
sx q[3];
rz(-1.2983373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1718488) q[0];
sx q[0];
rz(-1.8010704) q[0];
sx q[0];
rz(0.44542435) q[0];
rz(2.716966) q[1];
sx q[1];
rz(-1.5716962) q[1];
sx q[1];
rz(-0.62579036) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15881702) q[0];
sx q[0];
rz(-0.80251575) q[0];
sx q[0];
rz(-0.83440749) q[0];
rz(-1.884787) q[2];
sx q[2];
rz(-2.6382448) q[2];
sx q[2];
rz(-1.7679917) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6546331) q[1];
sx q[1];
rz(-0.59096293) q[1];
sx q[1];
rz(-2.0396359) q[1];
rz(-pi) q[2];
rz(1.6091691) q[3];
sx q[3];
rz(-1.4202446) q[3];
sx q[3];
rz(0.90177075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.17973913) q[2];
sx q[2];
rz(-2.1135795) q[2];
sx q[2];
rz(-1.4614159) q[2];
rz(0.05750582) q[3];
sx q[3];
rz(-1.688136) q[3];
sx q[3];
rz(0.97829515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4089324) q[0];
sx q[0];
rz(-1.4807777) q[0];
sx q[0];
rz(-2.5430191) q[0];
rz(-0.21845017) q[1];
sx q[1];
rz(-1.5427867) q[1];
sx q[1];
rz(-1.3839518) q[1];
rz(-1.6545532) q[2];
sx q[2];
rz(-1.0261921) q[2];
sx q[2];
rz(1.5941317) q[2];
rz(-2.3575847) q[3];
sx q[3];
rz(-1.5386469) q[3];
sx q[3];
rz(0.28924573) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
