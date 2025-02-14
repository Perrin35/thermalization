OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8074789) q[0];
sx q[0];
rz(2.4456094) q[0];
sx q[0];
rz(13.453475) q[0];
rz(2.4298985) q[1];
sx q[1];
rz(-2.0790172) q[1];
sx q[1];
rz(0.33023155) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4667128) q[0];
sx q[0];
rz(-1.1322563) q[0];
sx q[0];
rz(-0.86535485) q[0];
rz(2.5895533) q[2];
sx q[2];
rz(-0.62969724) q[2];
sx q[2];
rz(-0.80314595) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6373693) q[1];
sx q[1];
rz(-2.0614701) q[1];
sx q[1];
rz(-1.4549903) q[1];
x q[2];
rz(1.4012778) q[3];
sx q[3];
rz(-2.7054686) q[3];
sx q[3];
rz(-0.90986246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3868788) q[2];
sx q[2];
rz(-2.7140996) q[2];
sx q[2];
rz(1.0270366) q[2];
rz(-0.88296452) q[3];
sx q[3];
rz(-1.797978) q[3];
sx q[3];
rz(2.5880255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7850194) q[0];
sx q[0];
rz(-0.078332575) q[0];
sx q[0];
rz(2.8739492) q[0];
rz(-1.74125) q[1];
sx q[1];
rz(-2.5282271) q[1];
sx q[1];
rz(-2.5852481) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2927383) q[0];
sx q[0];
rz(-1.8195099) q[0];
sx q[0];
rz(1.066904) q[0];
rz(-pi) q[1];
x q[1];
rz(1.914996) q[2];
sx q[2];
rz(-0.89620334) q[2];
sx q[2];
rz(-0.062002484) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2399064) q[1];
sx q[1];
rz(-2.1237897) q[1];
sx q[1];
rz(2.6917975) q[1];
x q[2];
rz(-0.90640599) q[3];
sx q[3];
rz(-1.0944546) q[3];
sx q[3];
rz(1.2211459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9050425) q[2];
sx q[2];
rz(-1.6908129) q[2];
sx q[2];
rz(1.2518008) q[2];
rz(1.8768138) q[3];
sx q[3];
rz(-0.71139657) q[3];
sx q[3];
rz(2.0701764) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65997893) q[0];
sx q[0];
rz(-2.5857506) q[0];
sx q[0];
rz(0.83928338) q[0];
rz(1.5045769) q[1];
sx q[1];
rz(-1.4312276) q[1];
sx q[1];
rz(-1.379871) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31034841) q[0];
sx q[0];
rz(-1.6524685) q[0];
sx q[0];
rz(-1.6947338) q[0];
rz(-pi) q[1];
rz(2.1160385) q[2];
sx q[2];
rz(-1.7548867) q[2];
sx q[2];
rz(1.7005526) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.38244694) q[1];
sx q[1];
rz(-2.5576572) q[1];
sx q[1];
rz(0.96453519) q[1];
x q[2];
rz(2.1228509) q[3];
sx q[3];
rz(-2.0716487) q[3];
sx q[3];
rz(-0.93397442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69924361) q[2];
sx q[2];
rz(-1.2171429) q[2];
sx q[2];
rz(0.71845636) q[2];
rz(0.62026223) q[3];
sx q[3];
rz(-2.2060427) q[3];
sx q[3];
rz(0.81234318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6830867) q[0];
sx q[0];
rz(-1.2309256) q[0];
sx q[0];
rz(1.7339535) q[0];
rz(-2.5900224) q[1];
sx q[1];
rz(-3.0307814) q[1];
sx q[1];
rz(-1.3161906) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16770076) q[0];
sx q[0];
rz(-0.97588723) q[0];
sx q[0];
rz(-0.3428726) q[0];
x q[1];
rz(1.614093) q[2];
sx q[2];
rz(-0.27972066) q[2];
sx q[2];
rz(-2.2808275) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11034695) q[1];
sx q[1];
rz(-1.7025885) q[1];
sx q[1];
rz(-2.736997) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67263453) q[3];
sx q[3];
rz(-0.97390538) q[3];
sx q[3];
rz(-1.8565282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59336415) q[2];
sx q[2];
rz(-2.7061988) q[2];
sx q[2];
rz(1.0080053) q[2];
rz(-2.8625782) q[3];
sx q[3];
rz(-0.81441003) q[3];
sx q[3];
rz(2.6830955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27610436) q[0];
sx q[0];
rz(-1.3300329) q[0];
sx q[0];
rz(-0.27648595) q[0];
rz(2.8522988) q[1];
sx q[1];
rz(-1.9464867) q[1];
sx q[1];
rz(0.2624661) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.58231) q[0];
sx q[0];
rz(-1.773343) q[0];
sx q[0];
rz(1.3900533) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97945531) q[2];
sx q[2];
rz(-1.5618513) q[2];
sx q[2];
rz(1.1696377) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.20690675) q[1];
sx q[1];
rz(-2.1175724) q[1];
sx q[1];
rz(1.9711167) q[1];
rz(-pi) q[2];
rz(2.9093303) q[3];
sx q[3];
rz(-2.5506936) q[3];
sx q[3];
rz(-0.10400203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.7485973) q[2];
sx q[2];
rz(-0.59610072) q[2];
sx q[2];
rz(-1.7945012) q[2];
rz(-0.10284452) q[3];
sx q[3];
rz(-1.2343854) q[3];
sx q[3];
rz(-1.6245406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0031925072) q[0];
sx q[0];
rz(-1.6380558) q[0];
sx q[0];
rz(-0.060977161) q[0];
rz(1.2334476) q[1];
sx q[1];
rz(-1.7465218) q[1];
sx q[1];
rz(1.6159509) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71919854) q[0];
sx q[0];
rz(-2.6073747) q[0];
sx q[0];
rz(-0.13636503) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6055008) q[2];
sx q[2];
rz(-2.3922709) q[2];
sx q[2];
rz(-2.9108436) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.96887302) q[1];
sx q[1];
rz(-2.8668826) q[1];
sx q[1];
rz(-1.151888) q[1];
x q[2];
rz(0.91904007) q[3];
sx q[3];
rz(-1.7036333) q[3];
sx q[3];
rz(-0.024452535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.47152758) q[2];
sx q[2];
rz(-1.3482886) q[2];
sx q[2];
rz(-1.8111551) q[2];
rz(2.0165675) q[3];
sx q[3];
rz(-0.87441134) q[3];
sx q[3];
rz(3.103638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9576981) q[0];
sx q[0];
rz(-1.6989468) q[0];
sx q[0];
rz(-0.25699064) q[0];
rz(-0.43486241) q[1];
sx q[1];
rz(-2.3915274) q[1];
sx q[1];
rz(-0.90352568) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5241476) q[0];
sx q[0];
rz(-0.092369583) q[0];
sx q[0];
rz(1.6156107) q[0];
x q[1];
rz(0.14638325) q[2];
sx q[2];
rz(-2.3902262) q[2];
sx q[2];
rz(-0.20359542) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7419646) q[1];
sx q[1];
rz(-0.077721462) q[1];
sx q[1];
rz(2.4844954) q[1];
rz(-2.1210455) q[3];
sx q[3];
rz(-1.9686832) q[3];
sx q[3];
rz(0.83845058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32288185) q[2];
sx q[2];
rz(-1.5078397) q[2];
sx q[2];
rz(-1.7203077) q[2];
rz(-0.70982248) q[3];
sx q[3];
rz(-1.1387419) q[3];
sx q[3];
rz(0.53409725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98599559) q[0];
sx q[0];
rz(-2.7105712) q[0];
sx q[0];
rz(2.0565597) q[0];
rz(2.494508) q[1];
sx q[1];
rz(-1.1754464) q[1];
sx q[1];
rz(0.064402493) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2875117) q[0];
sx q[0];
rz(-0.72923822) q[0];
sx q[0];
rz(-0.50042787) q[0];
rz(-2.2930458) q[2];
sx q[2];
rz(-2.1050958) q[2];
sx q[2];
rz(-1.9398745) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.459584) q[1];
sx q[1];
rz(-1.1136076) q[1];
sx q[1];
rz(-2.3383635) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91592083) q[3];
sx q[3];
rz(-0.43817156) q[3];
sx q[3];
rz(2.8000591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25962466) q[2];
sx q[2];
rz(-0.3825863) q[2];
sx q[2];
rz(-0.043370334) q[2];
rz(-0.67443332) q[3];
sx q[3];
rz(-1.784227) q[3];
sx q[3];
rz(1.9969214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6944273) q[0];
sx q[0];
rz(-0.58037102) q[0];
sx q[0];
rz(0.32661435) q[0];
rz(1.0940374) q[1];
sx q[1];
rz(-2.2973165) q[1];
sx q[1];
rz(2.6577267) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9897468) q[0];
sx q[0];
rz(-2.2156179) q[0];
sx q[0];
rz(2.487464) q[0];
x q[1];
rz(-2.7174453) q[2];
sx q[2];
rz(-1.3636936) q[2];
sx q[2];
rz(1.0754367) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.30712515) q[1];
sx q[1];
rz(-0.79453429) q[1];
sx q[1];
rz(1.267872) q[1];
rz(-pi) q[2];
rz(-2.8362379) q[3];
sx q[3];
rz(-2.476446) q[3];
sx q[3];
rz(-2.0572452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6924374) q[2];
sx q[2];
rz(-2.1026976) q[2];
sx q[2];
rz(0.078744002) q[2];
rz(0.76830831) q[3];
sx q[3];
rz(-1.2945622) q[3];
sx q[3];
rz(2.9505742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67097265) q[0];
sx q[0];
rz(-0.60523954) q[0];
sx q[0];
rz(-0.11669267) q[0];
rz(2.281588) q[1];
sx q[1];
rz(-1.8000894) q[1];
sx q[1];
rz(-0.87342993) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3610342) q[0];
sx q[0];
rz(-1.3300899) q[0];
sx q[0];
rz(2.3248501) q[0];
x q[1];
rz(1.8928746) q[2];
sx q[2];
rz(-1.0828185) q[2];
sx q[2];
rz(2.8805594) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2817677) q[1];
sx q[1];
rz(-2.0689397) q[1];
sx q[1];
rz(-1.7120275) q[1];
rz(-pi) q[2];
rz(3.0449405) q[3];
sx q[3];
rz(-2.5279867) q[3];
sx q[3];
rz(2.8190316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7200835) q[2];
sx q[2];
rz(-2.8953711) q[2];
sx q[2];
rz(2.1923547) q[2];
rz(1.8947489) q[3];
sx q[3];
rz(-1.4693762) q[3];
sx q[3];
rz(-0.59103549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.128189) q[0];
sx q[0];
rz(-2.1407776) q[0];
sx q[0];
rz(1.5110973) q[0];
rz(0.36698256) q[1];
sx q[1];
rz(-1.3584247) q[1];
sx q[1];
rz(-0.57986837) q[1];
rz(-0.30324528) q[2];
sx q[2];
rz(-1.1228391) q[2];
sx q[2];
rz(-0.11930677) q[2];
rz(-0.69042023) q[3];
sx q[3];
rz(-1.2391117) q[3];
sx q[3];
rz(1.5285872) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
