OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0722421) q[0];
sx q[0];
rz(-1.0538333) q[0];
sx q[0];
rz(-1.8928438) q[0];
rz(-1.2381923) q[1];
sx q[1];
rz(-2.7014531) q[1];
sx q[1];
rz(1.8990489) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36168594) q[0];
sx q[0];
rz(-2.4324327) q[0];
sx q[0];
rz(0.61320029) q[0];
rz(-0.79808198) q[2];
sx q[2];
rz(-1.6387562) q[2];
sx q[2];
rz(-0.77347212) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.93852167) q[1];
sx q[1];
rz(-1.2483828) q[1];
sx q[1];
rz(-1.5702269) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60939021) q[3];
sx q[3];
rz(-1.2838138) q[3];
sx q[3];
rz(-0.089394102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0480463) q[2];
sx q[2];
rz(-2.8600433) q[2];
sx q[2];
rz(-1.2151037) q[2];
rz(-1.9563227) q[3];
sx q[3];
rz(-2.2351041) q[3];
sx q[3];
rz(-1.5989446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.1854061) q[0];
sx q[0];
rz(-2.7423999) q[0];
sx q[0];
rz(-1.4584374) q[0];
rz(2.9996808) q[1];
sx q[1];
rz(-1.2071573) q[1];
sx q[1];
rz(1.2123607) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1233009) q[0];
sx q[0];
rz(-2.2102076) q[0];
sx q[0];
rz(0.42808867) q[0];
x q[1];
rz(-1.5734912) q[2];
sx q[2];
rz(-0.97866733) q[2];
sx q[2];
rz(0.58838974) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4184561) q[1];
sx q[1];
rz(-0.89913705) q[1];
sx q[1];
rz(3.0211012) q[1];
x q[2];
rz(-0.47745269) q[3];
sx q[3];
rz(-2.4702854) q[3];
sx q[3];
rz(0.48760133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.045018) q[2];
sx q[2];
rz(-0.18288945) q[2];
sx q[2];
rz(-2.1796687) q[2];
rz(-2.9436881) q[3];
sx q[3];
rz(-1.3062545) q[3];
sx q[3];
rz(1.643868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4417931) q[0];
sx q[0];
rz(-2.4213591) q[0];
sx q[0];
rz(1.066347) q[0];
rz(-0.20866808) q[1];
sx q[1];
rz(-2.6228948) q[1];
sx q[1];
rz(0.64812237) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3828165) q[0];
sx q[0];
rz(-0.99446873) q[0];
sx q[0];
rz(-3.0786425) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0096613) q[2];
sx q[2];
rz(-2.6474617) q[2];
sx q[2];
rz(-3.1199093) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.40884205) q[1];
sx q[1];
rz(-2.2957845) q[1];
sx q[1];
rz(0.50028657) q[1];
rz(-pi) q[2];
rz(1.6567635) q[3];
sx q[3];
rz(-2.2760165) q[3];
sx q[3];
rz(-1.41627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2842747) q[2];
sx q[2];
rz(-1.6782574) q[2];
sx q[2];
rz(2.5939482) q[2];
rz(1.0037496) q[3];
sx q[3];
rz(-2.818483) q[3];
sx q[3];
rz(-2.9799262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6169287) q[0];
sx q[0];
rz(-2.014761) q[0];
sx q[0];
rz(-0.3666077) q[0];
rz(-1.2855351) q[1];
sx q[1];
rz(-1.5204241) q[1];
sx q[1];
rz(-2.3191648) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4124968) q[0];
sx q[0];
rz(-1.3733882) q[0];
sx q[0];
rz(1.4928994) q[0];
x q[1];
rz(-0.16777487) q[2];
sx q[2];
rz(-1.3047403) q[2];
sx q[2];
rz(2.7573836) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.3318204) q[1];
sx q[1];
rz(-2.0533516) q[1];
sx q[1];
rz(2.3259374) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0068914135) q[3];
sx q[3];
rz(-0.27537307) q[3];
sx q[3];
rz(1.7100348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.6377247) q[2];
sx q[2];
rz(-0.34054264) q[2];
sx q[2];
rz(1.5857504) q[2];
rz(-1.2268892) q[3];
sx q[3];
rz(-2.5033958) q[3];
sx q[3];
rz(1.3251806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1894839) q[0];
sx q[0];
rz(-2.0229078) q[0];
sx q[0];
rz(0.9504016) q[0];
rz(0.19418007) q[1];
sx q[1];
rz(-1.8870528) q[1];
sx q[1];
rz(2.8607184) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1732728) q[0];
sx q[0];
rz(-1.8283947) q[0];
sx q[0];
rz(-2.2031242) q[0];
x q[1];
rz(0.083706972) q[2];
sx q[2];
rz(-0.14995126) q[2];
sx q[2];
rz(2.170495) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7627613) q[1];
sx q[1];
rz(-1.4939337) q[1];
sx q[1];
rz(-0.084048653) q[1];
x q[2];
rz(-2.836197) q[3];
sx q[3];
rz(-1.0825233) q[3];
sx q[3];
rz(2.0528169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.48996553) q[2];
sx q[2];
rz(-1.6941035) q[2];
sx q[2];
rz(-0.70072407) q[2];
rz(-1.0939595) q[3];
sx q[3];
rz(-2.14812) q[3];
sx q[3];
rz(2.3195364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99264282) q[0];
sx q[0];
rz(-1.0473017) q[0];
sx q[0];
rz(0.5067504) q[0];
rz(2.0172334) q[1];
sx q[1];
rz(-1.1686814) q[1];
sx q[1];
rz(2.7728424) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0848727) q[0];
sx q[0];
rz(-0.73000008) q[0];
sx q[0];
rz(-1.8542791) q[0];
rz(-pi) q[1];
rz(0.36823057) q[2];
sx q[2];
rz(-0.86129649) q[2];
sx q[2];
rz(-1.8031098) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3430887) q[1];
sx q[1];
rz(-1.868416) q[1];
sx q[1];
rz(-0.095620015) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7140824) q[3];
sx q[3];
rz(-1.9906133) q[3];
sx q[3];
rz(1.9577718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2345978) q[2];
sx q[2];
rz(-2.5556421) q[2];
sx q[2];
rz(2.9242945) q[2];
rz(-1.4761188) q[3];
sx q[3];
rz(-0.81513351) q[3];
sx q[3];
rz(0.34172094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9076964) q[0];
sx q[0];
rz(-2.8192769) q[0];
sx q[0];
rz(0.79793683) q[0];
rz(1.7012874) q[1];
sx q[1];
rz(-2.1013575) q[1];
sx q[1];
rz(2.5247916) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26825702) q[0];
sx q[0];
rz(-0.63870478) q[0];
sx q[0];
rz(-0.12715841) q[0];
x q[1];
rz(1.2347414) q[2];
sx q[2];
rz(-1.0288887) q[2];
sx q[2];
rz(1.0934747) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4492256) q[1];
sx q[1];
rz(-1.7150075) q[1];
sx q[1];
rz(1.7768136) q[1];
rz(-pi) q[2];
rz(0.74941483) q[3];
sx q[3];
rz(-2.2843369) q[3];
sx q[3];
rz(-2.7593222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1282318) q[2];
sx q[2];
rz(-1.9888473) q[2];
sx q[2];
rz(0.70179233) q[2];
rz(2.2026786) q[3];
sx q[3];
rz(-2.2434668) q[3];
sx q[3];
rz(-1.8963337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.668648) q[0];
sx q[0];
rz(-2.9845147) q[0];
sx q[0];
rz(0.12338403) q[0];
rz(-0.75048796) q[1];
sx q[1];
rz(-1.6672983) q[1];
sx q[1];
rz(2.1231245) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3833419) q[0];
sx q[0];
rz(-1.4660379) q[0];
sx q[0];
rz(-1.7927756) q[0];
rz(-2.0872714) q[2];
sx q[2];
rz(-1.9811842) q[2];
sx q[2];
rz(-3.0832727) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97216304) q[1];
sx q[1];
rz(-2.0561245) q[1];
sx q[1];
rz(1.6134439) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94208053) q[3];
sx q[3];
rz(-1.6438369) q[3];
sx q[3];
rz(-1.4062509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5156775) q[2];
sx q[2];
rz(-1.6720142) q[2];
sx q[2];
rz(0.53543004) q[2];
rz(-1.2299906) q[3];
sx q[3];
rz(-1.001469) q[3];
sx q[3];
rz(-0.011822239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5984421) q[0];
sx q[0];
rz(-0.66023985) q[0];
sx q[0];
rz(-1.2793596) q[0];
rz(1.8774425) q[1];
sx q[1];
rz(-1.2707571) q[1];
sx q[1];
rz(0.15170161) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3705419) q[0];
sx q[0];
rz(-0.40239516) q[0];
sx q[0];
rz(2.9950525) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5506641) q[2];
sx q[2];
rz(-1.4369471) q[2];
sx q[2];
rz(1.0878022) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.38167414) q[1];
sx q[1];
rz(-1.2592788) q[1];
sx q[1];
rz(-2.1550234) q[1];
rz(-pi) q[2];
rz(-0.72355481) q[3];
sx q[3];
rz(-1.1465985) q[3];
sx q[3];
rz(1.8198609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5558527) q[2];
sx q[2];
rz(-0.93583217) q[2];
sx q[2];
rz(-1.2305416) q[2];
rz(1.3452283) q[3];
sx q[3];
rz(-1.3627005) q[3];
sx q[3];
rz(1.2983373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96974385) q[0];
sx q[0];
rz(-1.3405223) q[0];
sx q[0];
rz(-2.6961683) q[0];
rz(0.42462665) q[1];
sx q[1];
rz(-1.5716962) q[1];
sx q[1];
rz(-2.5158023) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75794894) q[0];
sx q[0];
rz(-1.008908) q[0];
sx q[0];
rz(0.60737078) q[0];
rz(-pi) q[1];
rz(2.0532512) q[2];
sx q[2];
rz(-1.7203334) q[2];
sx q[2];
rz(-2.66726) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6596131) q[1];
sx q[1];
rz(-1.8252884) q[1];
sx q[1];
rz(-1.031443) q[1];
x q[2];
rz(0.24769737) q[3];
sx q[3];
rz(-2.9862635) q[3];
sx q[3];
rz(0.65117902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9618535) q[2];
sx q[2];
rz(-2.1135795) q[2];
sx q[2];
rz(-1.4614159) q[2];
rz(0.05750582) q[3];
sx q[3];
rz(-1.4534566) q[3];
sx q[3];
rz(2.1632975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4089324) q[0];
sx q[0];
rz(-1.6608149) q[0];
sx q[0];
rz(0.59857359) q[0];
rz(-0.21845017) q[1];
sx q[1];
rz(-1.5427867) q[1];
sx q[1];
rz(-1.3839518) q[1];
rz(-0.54616164) q[2];
sx q[2];
rz(-1.6424137) q[2];
sx q[2];
rz(0.066802468) q[2];
rz(0.045513734) q[3];
sx q[3];
rz(-0.78452605) q[3];
sx q[3];
rz(1.8922643) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
