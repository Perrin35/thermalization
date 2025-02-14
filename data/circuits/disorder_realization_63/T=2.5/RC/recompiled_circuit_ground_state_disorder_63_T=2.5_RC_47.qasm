OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11290057) q[0];
sx q[0];
rz(-0.24554645) q[0];
sx q[0];
rz(-0.052659642) q[0];
rz(-2.0242937) q[1];
sx q[1];
rz(2.8007562) q[1];
sx q[1];
rz(10.513289) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77086335) q[0];
sx q[0];
rz(-1.5077871) q[0];
sx q[0];
rz(1.9669238) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4289189) q[2];
sx q[2];
rz(-1.4617702) q[2];
sx q[2];
rz(3.1224868) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6793078) q[1];
sx q[1];
rz(-2.4801697) q[1];
sx q[1];
rz(0.2703305) q[1];
rz(-pi) q[2];
rz(1.151122) q[3];
sx q[3];
rz(-1.5670781) q[3];
sx q[3];
rz(-2.9633455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5071621) q[2];
sx q[2];
rz(-1.3619225) q[2];
sx q[2];
rz(0.6081028) q[2];
rz(-1.1242584) q[3];
sx q[3];
rz(-1.9367155) q[3];
sx q[3];
rz(-0.89157909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41481498) q[0];
sx q[0];
rz(-1.5730653) q[0];
sx q[0];
rz(-2.538105) q[0];
rz(-2.6314349) q[1];
sx q[1];
rz(-2.3699103) q[1];
sx q[1];
rz(2.1147494) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.802216) q[0];
sx q[0];
rz(-2.908704) q[0];
sx q[0];
rz(2.4075137) q[0];
rz(-pi) q[1];
rz(-0.81342621) q[2];
sx q[2];
rz(-2.1132838) q[2];
sx q[2];
rz(0.39219609) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.75637792) q[1];
sx q[1];
rz(-0.88895386) q[1];
sx q[1];
rz(0.67519501) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92225109) q[3];
sx q[3];
rz(-2.3664858) q[3];
sx q[3];
rz(2.4668281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5828731) q[2];
sx q[2];
rz(-1.9931953) q[2];
sx q[2];
rz(0.64644512) q[2];
rz(-1.2785814) q[3];
sx q[3];
rz(-2.1918112) q[3];
sx q[3];
rz(1.4518552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.932514) q[0];
sx q[0];
rz(-0.13020733) q[0];
sx q[0];
rz(1.5721488) q[0];
rz(0.34700829) q[1];
sx q[1];
rz(-1.6667112) q[1];
sx q[1];
rz(0.86403799) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4915159) q[0];
sx q[0];
rz(-1.361761) q[0];
sx q[0];
rz(2.5670426) q[0];
rz(2.6367497) q[2];
sx q[2];
rz(-2.6528203) q[2];
sx q[2];
rz(-2.8504643) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5123295) q[1];
sx q[1];
rz(-0.95712582) q[1];
sx q[1];
rz(0.83028173) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7517463) q[3];
sx q[3];
rz(-0.53314185) q[3];
sx q[3];
rz(-1.0764493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0916834) q[2];
sx q[2];
rz(-2.9361762) q[2];
sx q[2];
rz(-1.6845711) q[2];
rz(2.125804) q[3];
sx q[3];
rz(-0.53272811) q[3];
sx q[3];
rz(-0.31639019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24093534) q[0];
sx q[0];
rz(-0.37739402) q[0];
sx q[0];
rz(1.578791) q[0];
rz(-2.0511625) q[1];
sx q[1];
rz(-2.5999887) q[1];
sx q[1];
rz(-1.9900367) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8354412) q[0];
sx q[0];
rz(-0.62850289) q[0];
sx q[0];
rz(-2.7449904) q[0];
rz(-pi) q[1];
rz(2.3811023) q[2];
sx q[2];
rz(-1.2142688) q[2];
sx q[2];
rz(0.52887756) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5051401) q[1];
sx q[1];
rz(-0.52524127) q[1];
sx q[1];
rz(1.6401827) q[1];
rz(-pi) q[2];
rz(-2.1776803) q[3];
sx q[3];
rz(-1.5480032) q[3];
sx q[3];
rz(-0.43302872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84732032) q[2];
sx q[2];
rz(-2.5966094) q[2];
sx q[2];
rz(1.4997743) q[2];
rz(0.13449399) q[3];
sx q[3];
rz(-1.8650863) q[3];
sx q[3];
rz(-0.78181481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0624369) q[0];
sx q[0];
rz(-2.2298614) q[0];
sx q[0];
rz(0.4656747) q[0];
rz(-1.2767208) q[1];
sx q[1];
rz(-1.0036422) q[1];
sx q[1];
rz(2.1086955) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55556923) q[0];
sx q[0];
rz(-0.15842552) q[0];
sx q[0];
rz(1.5044295) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0014702) q[2];
sx q[2];
rz(-1.6701785) q[2];
sx q[2];
rz(-1.6665719) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.031449854) q[1];
sx q[1];
rz(-1.8217161) q[1];
sx q[1];
rz(-0.67278426) q[1];
rz(-pi) q[2];
rz(-2.752653) q[3];
sx q[3];
rz(-2.0802214) q[3];
sx q[3];
rz(-3.096599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.59139645) q[2];
sx q[2];
rz(-1.1893716) q[2];
sx q[2];
rz(2.2097394) q[2];
rz(-0.094680928) q[3];
sx q[3];
rz(-2.2414312) q[3];
sx q[3];
rz(-1.8315106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69022995) q[0];
sx q[0];
rz(-0.19915038) q[0];
sx q[0];
rz(-3.0354101) q[0];
rz(0.55446082) q[1];
sx q[1];
rz(-2.3531871) q[1];
sx q[1];
rz(-1.5415446) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6351955) q[0];
sx q[0];
rz(-1.9612086) q[0];
sx q[0];
rz(-1.5968538) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0511398) q[2];
sx q[2];
rz(-0.51562947) q[2];
sx q[2];
rz(1.7898498) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20075837) q[1];
sx q[1];
rz(-1.6929827) q[1];
sx q[1];
rz(-0.69735511) q[1];
rz(-pi) q[2];
rz(-0.21884905) q[3];
sx q[3];
rz(-2.8893724) q[3];
sx q[3];
rz(-0.24009934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9083378) q[2];
sx q[2];
rz(-2.4007863) q[2];
sx q[2];
rz(1.4976658) q[2];
rz(0.46105591) q[3];
sx q[3];
rz(-0.65270972) q[3];
sx q[3];
rz(1.3155931) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6719565) q[0];
sx q[0];
rz(-0.75103432) q[0];
sx q[0];
rz(1.0787971) q[0];
rz(2.5476593) q[1];
sx q[1];
rz(-2.1682231) q[1];
sx q[1];
rz(2.914391) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0086454) q[0];
sx q[0];
rz(-1.7194178) q[0];
sx q[0];
rz(-0.11283837) q[0];
x q[1];
rz(1.960177) q[2];
sx q[2];
rz(-0.35786483) q[2];
sx q[2];
rz(0.030454446) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1006443) q[1];
sx q[1];
rz(-2.4923263) q[1];
sx q[1];
rz(-1.0430286) q[1];
x q[2];
rz(1.9428857) q[3];
sx q[3];
rz(-1.6642906) q[3];
sx q[3];
rz(2.1402511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5220773) q[2];
sx q[2];
rz(-2.7232813) q[2];
sx q[2];
rz(0.72959161) q[2];
rz(1.6884165) q[3];
sx q[3];
rz(-1.4540693) q[3];
sx q[3];
rz(-0.61354536) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5486117) q[0];
sx q[0];
rz(-2.225003) q[0];
sx q[0];
rz(0.35162893) q[0];
rz(0.31373203) q[1];
sx q[1];
rz(-1.2064563) q[1];
sx q[1];
rz(1.7650013) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54395318) q[0];
sx q[0];
rz(-2.7386129) q[0];
sx q[0];
rz(3.0209994) q[0];
rz(-0.96504546) q[2];
sx q[2];
rz(-2.9403983) q[2];
sx q[2];
rz(2.4810042) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.87159705) q[1];
sx q[1];
rz(-1.0078127) q[1];
sx q[1];
rz(1.3857122) q[1];
rz(-pi) q[2];
rz(1.2010716) q[3];
sx q[3];
rz(-1.4784357) q[3];
sx q[3];
rz(-2.8522688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5967963) q[2];
sx q[2];
rz(-0.72303855) q[2];
sx q[2];
rz(-1.6205622) q[2];
rz(-2.9936786) q[3];
sx q[3];
rz(-1.7971797) q[3];
sx q[3];
rz(-1.773905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3677597) q[0];
sx q[0];
rz(-1.2667043) q[0];
sx q[0];
rz(1.897478) q[0];
rz(-1.7077839) q[1];
sx q[1];
rz(-1.560874) q[1];
sx q[1];
rz(0.22020766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51802902) q[0];
sx q[0];
rz(-1.2014148) q[0];
sx q[0];
rz(-0.9250888) q[0];
rz(-pi) q[1];
rz(-2.037165) q[2];
sx q[2];
rz(-1.4254071) q[2];
sx q[2];
rz(-0.25823247) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0427093) q[1];
sx q[1];
rz(-2.115887) q[1];
sx q[1];
rz(0.49867757) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4768019) q[3];
sx q[3];
rz(-1.4729744) q[3];
sx q[3];
rz(1.0536915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1844909) q[2];
sx q[2];
rz(-1.9944921) q[2];
sx q[2];
rz(1.6415049) q[2];
rz(-3.0575276) q[3];
sx q[3];
rz(-1.2596687) q[3];
sx q[3];
rz(-0.070092289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4641007) q[0];
sx q[0];
rz(-0.75834948) q[0];
sx q[0];
rz(-0.41859928) q[0];
rz(-2.1981926) q[1];
sx q[1];
rz(-1.489233) q[1];
sx q[1];
rz(0.30280608) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26696268) q[0];
sx q[0];
rz(-1.4041252) q[0];
sx q[0];
rz(3.002949) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9507016) q[2];
sx q[2];
rz(-1.7451236) q[2];
sx q[2];
rz(-0.088100351) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8157685) q[1];
sx q[1];
rz(-0.093397141) q[1];
sx q[1];
rz(0.48431046) q[1];
rz(-pi) q[2];
rz(-1.2556158) q[3];
sx q[3];
rz(-2.7806296) q[3];
sx q[3];
rz(-0.40004594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6454978) q[2];
sx q[2];
rz(-0.87254137) q[2];
sx q[2];
rz(-2.3910451) q[2];
rz(1.6802855) q[3];
sx q[3];
rz(-0.76992005) q[3];
sx q[3];
rz(-2.6245978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1207598) q[0];
sx q[0];
rz(-1.5626386) q[0];
sx q[0];
rz(1.563969) q[0];
rz(-1.2278521) q[1];
sx q[1];
rz(-0.49929437) q[1];
sx q[1];
rz(1.1566537) q[1];
rz(1.7534134) q[2];
sx q[2];
rz(-1.0116521) q[2];
sx q[2];
rz(0.73965363) q[2];
rz(1.8932992) q[3];
sx q[3];
rz(-1.396198) q[3];
sx q[3];
rz(-0.085343501) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
