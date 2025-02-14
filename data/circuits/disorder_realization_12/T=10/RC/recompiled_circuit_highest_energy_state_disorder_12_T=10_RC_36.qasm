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
rz(2.501261) q[0];
sx q[0];
rz(-2.2202272) q[0];
sx q[0];
rz(1.5322354) q[0];
rz(-2.9345203) q[1];
sx q[1];
rz(-1.1595885) q[1];
sx q[1];
rz(-1.6699189) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4764413) q[0];
sx q[0];
rz(-1.2560606) q[0];
sx q[0];
rz(-1.2405618) q[0];
x q[1];
rz(-0.70090909) q[2];
sx q[2];
rz(-2.1953088) q[2];
sx q[2];
rz(-2.8025742) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.5122657) q[1];
sx q[1];
rz(-1.6696603) q[1];
sx q[1];
rz(1.3310286) q[1];
rz(-0.59572409) q[3];
sx q[3];
rz(-1.7215773) q[3];
sx q[3];
rz(0.25546701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.49429911) q[2];
sx q[2];
rz(-2.1493201) q[2];
sx q[2];
rz(2.5956019) q[2];
rz(2.2030988) q[3];
sx q[3];
rz(-2.6280845) q[3];
sx q[3];
rz(0.45309666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0168125) q[0];
sx q[0];
rz(-2.0240968) q[0];
sx q[0];
rz(-0.6413396) q[0];
rz(1.9524139) q[1];
sx q[1];
rz(-2.7180505) q[1];
sx q[1];
rz(2.0372527) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76726145) q[0];
sx q[0];
rz(-0.64909726) q[0];
sx q[0];
rz(0.24477203) q[0];
rz(-pi) q[1];
rz(2.8139427) q[2];
sx q[2];
rz(-2.0830685) q[2];
sx q[2];
rz(1.22359) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1321261) q[1];
sx q[1];
rz(-2.3915572) q[1];
sx q[1];
rz(-2.2723531) q[1];
rz(-pi) q[2];
rz(-1.6112174) q[3];
sx q[3];
rz(-1.4835525) q[3];
sx q[3];
rz(-0.62007346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7426593) q[2];
sx q[2];
rz(-1.089596) q[2];
sx q[2];
rz(0.99417865) q[2];
rz(1.582877) q[3];
sx q[3];
rz(-2.4623058) q[3];
sx q[3];
rz(0.5510295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-3.1347374) q[0];
sx q[0];
rz(-1.0900494) q[0];
sx q[0];
rz(-2.5110974) q[0];
rz(0.14671239) q[1];
sx q[1];
rz(-1.2598597) q[1];
sx q[1];
rz(2.2781118) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3849626) q[0];
sx q[0];
rz(-2.4061086) q[0];
sx q[0];
rz(1.1138659) q[0];
rz(-pi) q[1];
rz(0.21196694) q[2];
sx q[2];
rz(-1.2951998) q[2];
sx q[2];
rz(2.8805062) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.59592225) q[1];
sx q[1];
rz(-1.1356259) q[1];
sx q[1];
rz(0.11637139) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9894838) q[3];
sx q[3];
rz(-0.31732355) q[3];
sx q[3];
rz(-3.139024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8063987) q[2];
sx q[2];
rz(-1.6890182) q[2];
sx q[2];
rz(-2.2550968) q[2];
rz(0.42996201) q[3];
sx q[3];
rz(-2.6468247) q[3];
sx q[3];
rz(0.89216843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48915136) q[0];
sx q[0];
rz(-2.8471071) q[0];
sx q[0];
rz(-2.6926706) q[0];
rz(2.4849675) q[1];
sx q[1];
rz(-1.4337599) q[1];
sx q[1];
rz(2.8112559) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2097004) q[0];
sx q[0];
rz(-0.39474129) q[0];
sx q[0];
rz(-2.1861211) q[0];
x q[1];
rz(2.9135111) q[2];
sx q[2];
rz(-1.9794399) q[2];
sx q[2];
rz(0.064106031) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8395665) q[1];
sx q[1];
rz(-2.6077787) q[1];
sx q[1];
rz(-2.1697976) q[1];
rz(-pi) q[2];
rz(1.6449241) q[3];
sx q[3];
rz(-0.12136689) q[3];
sx q[3];
rz(1.6182155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0535447) q[2];
sx q[2];
rz(-1.6036754) q[2];
sx q[2];
rz(-0.21928445) q[2];
rz(-2.3413279) q[3];
sx q[3];
rz(-2.181668) q[3];
sx q[3];
rz(2.3805591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6680172) q[0];
sx q[0];
rz(-1.1290843) q[0];
sx q[0];
rz(-1.4833204) q[0];
rz(2.5738916) q[1];
sx q[1];
rz(-1.9300108) q[1];
sx q[1];
rz(1.6426881) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53500861) q[0];
sx q[0];
rz(-0.55130115) q[0];
sx q[0];
rz(2.8583741) q[0];
rz(-pi) q[1];
rz(-0.097857006) q[2];
sx q[2];
rz(-1.5598394) q[2];
sx q[2];
rz(-1.9734073) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0755456) q[1];
sx q[1];
rz(-2.4986364) q[1];
sx q[1];
rz(2.8075904) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12171774) q[3];
sx q[3];
rz(-1.7608402) q[3];
sx q[3];
rz(2.7167883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6296926) q[2];
sx q[2];
rz(-0.7581768) q[2];
sx q[2];
rz(-0.32621113) q[2];
rz(-3.015226) q[3];
sx q[3];
rz(-2.5765403) q[3];
sx q[3];
rz(0.27157426) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3417974) q[0];
sx q[0];
rz(-1.3038776) q[0];
sx q[0];
rz(-0.52249348) q[0];
rz(-0.74784589) q[1];
sx q[1];
rz(-1.6035085) q[1];
sx q[1];
rz(-1.984367) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.889288) q[0];
sx q[0];
rz(-1.2784908) q[0];
sx q[0];
rz(2.0415123) q[0];
rz(-pi) q[1];
rz(-0.87027486) q[2];
sx q[2];
rz(-0.57862568) q[2];
sx q[2];
rz(3.0491997) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.044837601) q[1];
sx q[1];
rz(-1.4787349) q[1];
sx q[1];
rz(1.4274867) q[1];
rz(-pi) q[2];
rz(-2.8145267) q[3];
sx q[3];
rz(-1.6432135) q[3];
sx q[3];
rz(-1.4625957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.130827) q[2];
sx q[2];
rz(-0.5126493) q[2];
sx q[2];
rz(-1.1831247) q[2];
rz(0.042923953) q[3];
sx q[3];
rz(-1.639651) q[3];
sx q[3];
rz(0.24421282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6965028) q[0];
sx q[0];
rz(-0.88560167) q[0];
sx q[0];
rz(2.4687299) q[0];
rz(0.54975763) q[1];
sx q[1];
rz(-1.8472698) q[1];
sx q[1];
rz(-1.4901644) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5888118) q[0];
sx q[0];
rz(-0.59156617) q[0];
sx q[0];
rz(2.8179396) q[0];
x q[1];
rz(-2.0668903) q[2];
sx q[2];
rz(-2.0286244) q[2];
sx q[2];
rz(1.7448448) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6770391) q[1];
sx q[1];
rz(-2.0050075) q[1];
sx q[1];
rz(-1.6766816) q[1];
rz(-pi) q[2];
rz(2.1792322) q[3];
sx q[3];
rz(-0.058789805) q[3];
sx q[3];
rz(2.1421681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0979536) q[2];
sx q[2];
rz(-2.592228) q[2];
sx q[2];
rz(-1.1489541) q[2];
rz(2.7502381) q[3];
sx q[3];
rz(-1.2201744) q[3];
sx q[3];
rz(0.84827387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.5963762) q[0];
sx q[0];
rz(-1.2059809) q[0];
sx q[0];
rz(2.4429831) q[0];
rz(0.85083234) q[1];
sx q[1];
rz(-0.60092503) q[1];
sx q[1];
rz(-0.94815475) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7018775) q[0];
sx q[0];
rz(-2.2902459) q[0];
sx q[0];
rz(1.6802358) q[0];
x q[1];
rz(-3.0879232) q[2];
sx q[2];
rz(-1.2415774) q[2];
sx q[2];
rz(-1.366113) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.588447) q[1];
sx q[1];
rz(-1.1678034) q[1];
sx q[1];
rz(-0.95882596) q[1];
rz(-pi) q[2];
rz(2.3513673) q[3];
sx q[3];
rz(-2.2314615) q[3];
sx q[3];
rz(3.1243711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66905388) q[2];
sx q[2];
rz(-1.6927745) q[2];
sx q[2];
rz(-1.4581468) q[2];
rz(1.091188) q[3];
sx q[3];
rz(-2.2444221) q[3];
sx q[3];
rz(-0.18479656) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3882196) q[0];
sx q[0];
rz(-2.9599157) q[0];
sx q[0];
rz(-1.3004119) q[0];
rz(2.6822207) q[1];
sx q[1];
rz(-1.4540902) q[1];
sx q[1];
rz(-2.557911) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65410173) q[0];
sx q[0];
rz(-2.0493453) q[0];
sx q[0];
rz(2.3951981) q[0];
x q[1];
rz(-1.6064209) q[2];
sx q[2];
rz(-1.5322663) q[2];
sx q[2];
rz(0.57254475) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3724842) q[1];
sx q[1];
rz(-2.5021837) q[1];
sx q[1];
rz(0.95692371) q[1];
x q[2];
rz(-2.0784057) q[3];
sx q[3];
rz(-0.48590966) q[3];
sx q[3];
rz(2.8775281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5251081) q[2];
sx q[2];
rz(-2.8664092) q[2];
sx q[2];
rz(0.52213651) q[2];
rz(-2.3790242) q[3];
sx q[3];
rz(-1.6430166) q[3];
sx q[3];
rz(-0.34408072) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46227118) q[0];
sx q[0];
rz(-1.9336047) q[0];
sx q[0];
rz(0.59360498) q[0];
rz(-2.4484334) q[1];
sx q[1];
rz(-2.3489372) q[1];
sx q[1];
rz(2.3902182) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40461883) q[0];
sx q[0];
rz(-1.3122137) q[0];
sx q[0];
rz(0.41249852) q[0];
rz(2.2518806) q[2];
sx q[2];
rz(-1.0719704) q[2];
sx q[2];
rz(2.6798583) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.33022949) q[1];
sx q[1];
rz(-0.6682446) q[1];
sx q[1];
rz(3.038637) q[1];
rz(-1.4937711) q[3];
sx q[3];
rz(-1.4024156) q[3];
sx q[3];
rz(-1.2428969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0005325) q[2];
sx q[2];
rz(-0.70912051) q[2];
sx q[2];
rz(3.1179023) q[2];
rz(2.0247816) q[3];
sx q[3];
rz(-1.4477372) q[3];
sx q[3];
rz(1.1344596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5524207) q[0];
sx q[0];
rz(-1.9235274) q[0];
sx q[0];
rz(-2.0148475) q[0];
rz(1.8213656) q[1];
sx q[1];
rz(-1.2757433) q[1];
sx q[1];
rz(-2.129) q[1];
rz(-2.7350551) q[2];
sx q[2];
rz(-2.5589522) q[2];
sx q[2];
rz(2.027579) q[2];
rz(-1.3397401) q[3];
sx q[3];
rz(-1.7920296) q[3];
sx q[3];
rz(-1.5248004) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
