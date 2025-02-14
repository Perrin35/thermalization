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
rz(-1.363938) q[0];
sx q[0];
rz(-2.7231556) q[0];
sx q[0];
rz(1.3016181) q[0];
rz(0.16297451) q[1];
sx q[1];
rz(-1.4459223) q[1];
sx q[1];
rz(0.016782848) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7201286) q[0];
sx q[0];
rz(-1.4961494) q[0];
sx q[0];
rz(1.9326841) q[0];
rz(-0.016525908) q[2];
sx q[2];
rz(-1.3912829) q[2];
sx q[2];
rz(-1.5915888) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.77474817) q[1];
sx q[1];
rz(-1.5758744) q[1];
sx q[1];
rz(-3.1366411) q[1];
x q[2];
rz(-0.34256012) q[3];
sx q[3];
rz(-2.9800219) q[3];
sx q[3];
rz(1.8752961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5825384) q[2];
sx q[2];
rz(-2.4808919) q[2];
sx q[2];
rz(-1.5622697) q[2];
rz(-2.9301379) q[3];
sx q[3];
rz(-3.1410757) q[3];
sx q[3];
rz(2.9805984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6120537) q[0];
sx q[0];
rz(-2.8654629) q[0];
sx q[0];
rz(1.8147234) q[0];
rz(-2.5621085) q[1];
sx q[1];
rz(-0.0038298413) q[1];
sx q[1];
rz(0.63900596) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9352183) q[0];
sx q[0];
rz(-0.91418302) q[0];
sx q[0];
rz(-2.4103863) q[0];
x q[1];
rz(-0.1370378) q[2];
sx q[2];
rz(-0.1220905) q[2];
sx q[2];
rz(3.0220495) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3957876) q[1];
sx q[1];
rz(-1.5823872) q[1];
sx q[1];
rz(-0.017000217) q[1];
rz(-pi) q[2];
rz(-0.08777471) q[3];
sx q[3];
rz(-2.3592161) q[3];
sx q[3];
rz(2.6915472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4352033) q[2];
sx q[2];
rz(-3.0053164) q[2];
sx q[2];
rz(1.603568) q[2];
rz(-1.5609353) q[3];
sx q[3];
rz(-3.1272562) q[3];
sx q[3];
rz(-3.1106136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8921709) q[0];
sx q[0];
rz(-0.51512655) q[0];
sx q[0];
rz(-2.7601335) q[0];
rz(2.4341266) q[1];
sx q[1];
rz(-0.019376945) q[1];
sx q[1];
rz(1.1245419) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3113925) q[0];
sx q[0];
rz(-1.8233577) q[0];
sx q[0];
rz(-0.24858944) q[0];
rz(-pi) q[1];
x q[1];
rz(0.025990268) q[2];
sx q[2];
rz(-1.6857237) q[2];
sx q[2];
rz(2.9851802) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.687011) q[1];
sx q[1];
rz(-3.0774766) q[1];
sx q[1];
rz(2.9518576) q[1];
x q[2];
rz(0.74060969) q[3];
sx q[3];
rz(-2.4428074) q[3];
sx q[3];
rz(1.6388338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4418929) q[2];
sx q[2];
rz(-3.1293588) q[2];
sx q[2];
rz(0.049467889) q[2];
rz(2.5345645) q[3];
sx q[3];
rz(-3.1403465) q[3];
sx q[3];
rz(-1.9342669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98168755) q[0];
sx q[0];
rz(-2.9732381) q[0];
sx q[0];
rz(-3.1244151) q[0];
rz(-0.29302868) q[1];
sx q[1];
rz(-0.79048645) q[1];
sx q[1];
rz(-1.5471829) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.806694) q[0];
sx q[0];
rz(-0.8282173) q[0];
sx q[0];
rz(2.3252299) q[0];
rz(-pi) q[1];
rz(-2.5998579) q[2];
sx q[2];
rz(-1.3770665) q[2];
sx q[2];
rz(1.3044747) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3694689) q[1];
sx q[1];
rz(-1.511206) q[1];
sx q[1];
rz(3.0173484) q[1];
rz(1.5805575) q[3];
sx q[3];
rz(-2.2167444) q[3];
sx q[3];
rz(-2.498359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8881417) q[2];
sx q[2];
rz(-2.6710822) q[2];
sx q[2];
rz(-0.68224254) q[2];
rz(0.055179723) q[3];
sx q[3];
rz(-3.1339055) q[3];
sx q[3];
rz(1.3050219) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3799915) q[0];
sx q[0];
rz(-0.43552265) q[0];
sx q[0];
rz(2.7165661) q[0];
rz(1.60166) q[1];
sx q[1];
rz(-0.48307499) q[1];
sx q[1];
rz(-0.81533283) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5402712) q[0];
sx q[0];
rz(-3.0794332) q[0];
sx q[0];
rz(0.17610733) q[0];
rz(-pi) q[1];
rz(-1.5940395) q[2];
sx q[2];
rz(-1.574062) q[2];
sx q[2];
rz(2.4583985) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.40774959) q[1];
sx q[1];
rz(-1.6899278) q[1];
sx q[1];
rz(1.6581737) q[1];
rz(-pi) q[2];
rz(-2.8132503) q[3];
sx q[3];
rz(-0.36188618) q[3];
sx q[3];
rz(2.7992424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1347947) q[2];
sx q[2];
rz(-0.012601348) q[2];
sx q[2];
rz(-1.4726144) q[2];
rz(-0.93004477) q[3];
sx q[3];
rz(-0.01447066) q[3];
sx q[3];
rz(2.3013733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3227661) q[0];
sx q[0];
rz(-3.1184986) q[0];
sx q[0];
rz(1.4267138) q[0];
rz(2.4115883) q[1];
sx q[1];
rz(-0.58861029) q[1];
sx q[1];
rz(2.0402562) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5752714) q[0];
sx q[0];
rz(-2.3744517) q[0];
sx q[0];
rz(1.9643892) q[0];
rz(1.3404714) q[2];
sx q[2];
rz(-1.4061951) q[2];
sx q[2];
rz(-2.4245461) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37820617) q[1];
sx q[1];
rz(-2.9855033) q[1];
sx q[1];
rz(-0.915145) q[1];
rz(-1.6953129) q[3];
sx q[3];
rz(-1.6027556) q[3];
sx q[3];
rz(0.65297844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.71658984) q[2];
sx q[2];
rz(-0.060881946) q[2];
sx q[2];
rz(1.8343743) q[2];
rz(-0.37846765) q[3];
sx q[3];
rz(-0.022947939) q[3];
sx q[3];
rz(-2.4881261) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3017479) q[0];
sx q[0];
rz(-1.2675588) q[0];
sx q[0];
rz(0.92754716) q[0];
rz(1.7842267) q[1];
sx q[1];
rz(-0.83186847) q[1];
sx q[1];
rz(1.5313139) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46857014) q[0];
sx q[0];
rz(-3.0971905) q[0];
sx q[0];
rz(-2.0973678) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8157829) q[2];
sx q[2];
rz(-1.9568482) q[2];
sx q[2];
rz(-0.99967848) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.011259638) q[1];
sx q[1];
rz(-1.5734871) q[1];
sx q[1];
rz(-1.6998768) q[1];
rz(-pi) q[2];
rz(2.82656) q[3];
sx q[3];
rz(-1.0803534) q[3];
sx q[3];
rz(-0.88446188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7808468) q[2];
sx q[2];
rz(-0.0043892269) q[2];
sx q[2];
rz(1.8878262) q[2];
rz(0.69416657) q[3];
sx q[3];
rz(-0.73918754) q[3];
sx q[3];
rz(0.23301253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53741443) q[0];
sx q[0];
rz(-1.0056714) q[0];
sx q[0];
rz(2.1042714) q[0];
rz(-1.6088156) q[1];
sx q[1];
rz(-0.2205801) q[1];
sx q[1];
rz(1.6745837) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5418422) q[0];
sx q[0];
rz(-1.5034202) q[0];
sx q[0];
rz(1.7662494) q[0];
rz(0.080341415) q[2];
sx q[2];
rz(-2.8626277) q[2];
sx q[2];
rz(2.7921048) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.26997494) q[1];
sx q[1];
rz(-1.5712357) q[1];
sx q[1];
rz(1.5696708) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6862304) q[3];
sx q[3];
rz(-1.4858559) q[3];
sx q[3];
rz(0.37173879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.011270114) q[2];
sx q[2];
rz(-0.20877561) q[2];
sx q[2];
rz(3.0976963) q[2];
rz(2.6210426) q[3];
sx q[3];
rz(-3.136941) q[3];
sx q[3];
rz(1.4588149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7875824) q[0];
sx q[0];
rz(-3.1391322) q[0];
sx q[0];
rz(-1.822923) q[0];
rz(1.4175381) q[1];
sx q[1];
rz(-0.28957614) q[1];
sx q[1];
rz(-1.5971378) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.330724) q[0];
sx q[0];
rz(-1.4557342) q[0];
sx q[0];
rz(2.6601391) q[0];
x q[1];
rz(0.6575281) q[2];
sx q[2];
rz(-1.7818461) q[2];
sx q[2];
rz(1.5653277) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.68543562) q[1];
sx q[1];
rz(-2.7815869) q[1];
sx q[1];
rz(-2.092157) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4539976) q[3];
sx q[3];
rz(-1.9904226) q[3];
sx q[3];
rz(0.18292038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.80457193) q[2];
sx q[2];
rz(-1.8420409) q[2];
sx q[2];
rz(-0.19680944) q[2];
rz(1.9491516) q[3];
sx q[3];
rz(-0.20659031) q[3];
sx q[3];
rz(-0.20096745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8404959) q[0];
sx q[0];
rz(-1.3571285) q[0];
sx q[0];
rz(1.949973) q[0];
rz(-1.5246897) q[1];
sx q[1];
rz(-0.646851) q[1];
sx q[1];
rz(1.5651388) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8218481) q[0];
sx q[0];
rz(-0.34043202) q[0];
sx q[0];
rz(-0.52297395) q[0];
x q[1];
rz(-1.1249816) q[2];
sx q[2];
rz(-0.49575323) q[2];
sx q[2];
rz(3.1204566) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8818156) q[1];
sx q[1];
rz(-1.5702973) q[1];
sx q[1];
rz(-0.00094836469) q[1];
rz(0.11348806) q[3];
sx q[3];
rz(-1.7120512) q[3];
sx q[3];
rz(0.065441386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9578751) q[2];
sx q[2];
rz(-0.58791939) q[2];
sx q[2];
rz(-1.4584165) q[2];
rz(3.1111187) q[3];
sx q[3];
rz(-3.1320429) q[3];
sx q[3];
rz(-0.20235801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6644345) q[0];
sx q[0];
rz(-1.8071334) q[0];
sx q[0];
rz(-1.4596756) q[0];
rz(1.5741813) q[1];
sx q[1];
rz(-1.8125143) q[1];
sx q[1];
rz(0.090851091) q[1];
rz(-1.6463584) q[2];
sx q[2];
rz(-1.5757558) q[2];
sx q[2];
rz(-1.4103665) q[2];
rz(-0.50441691) q[3];
sx q[3];
rz(-2.258959) q[3];
sx q[3];
rz(-2.7213982) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
