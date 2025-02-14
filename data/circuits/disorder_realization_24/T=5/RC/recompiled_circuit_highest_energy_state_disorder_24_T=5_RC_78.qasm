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
rz(0.46718207) q[0];
sx q[0];
rz(-1.2527569) q[0];
sx q[0];
rz(7.0897515) q[0];
rz(-0.84713495) q[1];
sx q[1];
rz(-0.48589125) q[1];
sx q[1];
rz(-2.2810305) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33968494) q[0];
sx q[0];
rz(-0.97754495) q[0];
sx q[0];
rz(0.16925933) q[0];
rz(-0.86463682) q[2];
sx q[2];
rz(-1.3792017) q[2];
sx q[2];
rz(-0.77048555) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7736176) q[1];
sx q[1];
rz(-2.8624075) q[1];
sx q[1];
rz(-0.18191819) q[1];
x q[2];
rz(-0.99765649) q[3];
sx q[3];
rz(-1.9966148) q[3];
sx q[3];
rz(-2.4906858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.90031558) q[2];
sx q[2];
rz(-0.52854717) q[2];
sx q[2];
rz(0.47145525) q[2];
rz(-2.6491162) q[3];
sx q[3];
rz(-1.9086647) q[3];
sx q[3];
rz(-1.2537778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8660698) q[0];
sx q[0];
rz(-0.38551426) q[0];
sx q[0];
rz(2.8634014) q[0];
rz(1.1909852) q[1];
sx q[1];
rz(-1.8089801) q[1];
sx q[1];
rz(-1.8461548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0012521) q[0];
sx q[0];
rz(-0.80632001) q[0];
sx q[0];
rz(2.4693523) q[0];
rz(1.5191742) q[2];
sx q[2];
rz(-1.4913017) q[2];
sx q[2];
rz(2.5148066) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.55838096) q[1];
sx q[1];
rz(-1.19046) q[1];
sx q[1];
rz(3.1391869) q[1];
rz(-pi) q[2];
rz(-0.94280394) q[3];
sx q[3];
rz(-2.0217388) q[3];
sx q[3];
rz(2.1357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.94464716) q[2];
sx q[2];
rz(-1.6702007) q[2];
sx q[2];
rz(-0.50986457) q[2];
rz(1.4465796) q[3];
sx q[3];
rz(-0.46049419) q[3];
sx q[3];
rz(-0.50486008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6079717) q[0];
sx q[0];
rz(-1.5190268) q[0];
sx q[0];
rz(2.1571889) q[0];
rz(-0.016117485) q[1];
sx q[1];
rz(-1.1382444) q[1];
sx q[1];
rz(1.3498397) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4348819) q[0];
sx q[0];
rz(-0.87315403) q[0];
sx q[0];
rz(1.7876704) q[0];
rz(-pi) q[1];
rz(2.1355633) q[2];
sx q[2];
rz(-1.2899219) q[2];
sx q[2];
rz(-1.8233521) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.74047725) q[1];
sx q[1];
rz(-2.4414032) q[1];
sx q[1];
rz(-3.0790331) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3758395) q[3];
sx q[3];
rz(-1.651305) q[3];
sx q[3];
rz(1.3020368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.71780378) q[2];
sx q[2];
rz(-0.70222792) q[2];
sx q[2];
rz(-0.88199893) q[2];
rz(-1.9574022) q[3];
sx q[3];
rz(-0.94112527) q[3];
sx q[3];
rz(0.73630303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0475567) q[0];
sx q[0];
rz(-1.9147669) q[0];
sx q[0];
rz(-0.077022821) q[0];
rz(0.19829622) q[1];
sx q[1];
rz(-1.501333) q[1];
sx q[1];
rz(0.1870627) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7884946) q[0];
sx q[0];
rz(-1.234963) q[0];
sx q[0];
rz(-1.737835) q[0];
rz(-pi) q[1];
rz(2.0346626) q[2];
sx q[2];
rz(-2.7527546) q[2];
sx q[2];
rz(2.1269682) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9976207) q[1];
sx q[1];
rz(-1.629843) q[1];
sx q[1];
rz(-2.4322492) q[1];
rz(1.2310394) q[3];
sx q[3];
rz(-1.2715142) q[3];
sx q[3];
rz(1.2359985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5932811) q[2];
sx q[2];
rz(-0.45479861) q[2];
sx q[2];
rz(-1.456267) q[2];
rz(0.71825394) q[3];
sx q[3];
rz(-1.7536438) q[3];
sx q[3];
rz(-2.3072402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2497571) q[0];
sx q[0];
rz(-1.3351853) q[0];
sx q[0];
rz(-0.43858132) q[0];
rz(-0.35722411) q[1];
sx q[1];
rz(-1.2780739) q[1];
sx q[1];
rz(1.2507778) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.253535) q[0];
sx q[0];
rz(-1.4708637) q[0];
sx q[0];
rz(1.1424095) q[0];
rz(1.8821472) q[2];
sx q[2];
rz(-1.6451453) q[2];
sx q[2];
rz(1.5386563) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.631613) q[1];
sx q[1];
rz(-2.5812979) q[1];
sx q[1];
rz(-2.3255682) q[1];
rz(2.881348) q[3];
sx q[3];
rz(-1.4422999) q[3];
sx q[3];
rz(-0.63333007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9095416) q[2];
sx q[2];
rz(-0.48018685) q[2];
sx q[2];
rz(-1.5117744) q[2];
rz(0.016228598) q[3];
sx q[3];
rz(-1.0954233) q[3];
sx q[3];
rz(0.79286638) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4474354) q[0];
sx q[0];
rz(-1.2801535) q[0];
sx q[0];
rz(-1.1496899) q[0];
rz(-0.42090526) q[1];
sx q[1];
rz(-1.4525388) q[1];
sx q[1];
rz(-0.044205753) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6290508) q[0];
sx q[0];
rz(-2.3568332) q[0];
sx q[0];
rz(0.46585887) q[0];
x q[1];
rz(-1.030613) q[2];
sx q[2];
rz(-1.9760796) q[2];
sx q[2];
rz(2.0342397) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.468535) q[1];
sx q[1];
rz(-1.6516764) q[1];
sx q[1];
rz(1.1799501) q[1];
rz(-1.8474691) q[3];
sx q[3];
rz(-1.0545316) q[3];
sx q[3];
rz(-0.54126213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1786903) q[2];
sx q[2];
rz(-2.2447605) q[2];
sx q[2];
rz(2.1964591) q[2];
rz(-1.6804228) q[3];
sx q[3];
rz(-1.1472568) q[3];
sx q[3];
rz(-0.51819658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1951676) q[0];
sx q[0];
rz(-2.8312046) q[0];
sx q[0];
rz(2.2921966) q[0];
rz(1.7646344) q[1];
sx q[1];
rz(-0.57146776) q[1];
sx q[1];
rz(-1.8570541) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3234069) q[0];
sx q[0];
rz(-2.0571097) q[0];
sx q[0];
rz(-2.0165927) q[0];
rz(-0.30751245) q[2];
sx q[2];
rz(-1.4685681) q[2];
sx q[2];
rz(-2.6999465) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9640935) q[1];
sx q[1];
rz(-2.0960137) q[1];
sx q[1];
rz(1.6119269) q[1];
rz(0.48619481) q[3];
sx q[3];
rz(-1.0665585) q[3];
sx q[3];
rz(-2.274226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4756502) q[2];
sx q[2];
rz(-1.7721756) q[2];
sx q[2];
rz(1.5744038) q[2];
rz(2.808908) q[3];
sx q[3];
rz(-1.0654457) q[3];
sx q[3];
rz(-2.1596597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5653) q[0];
sx q[0];
rz(-1.1197634) q[0];
sx q[0];
rz(2.0530307) q[0];
rz(-0.92900485) q[1];
sx q[1];
rz(-1.6477511) q[1];
sx q[1];
rz(-0.30379024) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32919183) q[0];
sx q[0];
rz(-2.441933) q[0];
sx q[0];
rz(-0.60075642) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4073885) q[2];
sx q[2];
rz(-1.4532928) q[2];
sx q[2];
rz(0.30899707) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6188083) q[1];
sx q[1];
rz(-1.6635824) q[1];
sx q[1];
rz(2.001808) q[1];
rz(2.0511384) q[3];
sx q[3];
rz(-1.0075724) q[3];
sx q[3];
rz(-2.8728812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2527689) q[2];
sx q[2];
rz(-1.4759651) q[2];
sx q[2];
rz(2.8775173) q[2];
rz(-2.0090328) q[3];
sx q[3];
rz(-2.8045636) q[3];
sx q[3];
rz(1.0549217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7842512) q[0];
sx q[0];
rz(-2.7989474) q[0];
sx q[0];
rz(-1.0145048) q[0];
rz(2.2134589) q[1];
sx q[1];
rz(-1.4200297) q[1];
sx q[1];
rz(2.6511505) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9336995) q[0];
sx q[0];
rz(-1.636278) q[0];
sx q[0];
rz(1.4235953) q[0];
rz(-pi) q[1];
rz(-2.7968452) q[2];
sx q[2];
rz(-1.9354068) q[2];
sx q[2];
rz(0.62275902) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0874493) q[1];
sx q[1];
rz(-1.6119401) q[1];
sx q[1];
rz(1.270134) q[1];
x q[2];
rz(1.0317627) q[3];
sx q[3];
rz(-1.3628873) q[3];
sx q[3];
rz(-0.62139213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8162615) q[2];
sx q[2];
rz(-1.0518495) q[2];
sx q[2];
rz(3.0391147) q[2];
rz(1.8898194) q[3];
sx q[3];
rz(-1.7779558) q[3];
sx q[3];
rz(1.4832835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3935299) q[0];
sx q[0];
rz(-3.0952125) q[0];
sx q[0];
rz(1.8512132) q[0];
rz(-2.6875467) q[1];
sx q[1];
rz(-1.8171277) q[1];
sx q[1];
rz(0.18347278) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35724872) q[0];
sx q[0];
rz(-1.5897449) q[0];
sx q[0];
rz(1.3878893) q[0];
rz(-pi) q[1];
x q[1];
rz(0.021963693) q[2];
sx q[2];
rz(-2.5434867) q[2];
sx q[2];
rz(2.919521) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80663348) q[1];
sx q[1];
rz(-2.0628961) q[1];
sx q[1];
rz(1.5978807) q[1];
x q[2];
rz(1.7137035) q[3];
sx q[3];
rz(-2.5710754) q[3];
sx q[3];
rz(-0.6541259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7196677) q[2];
sx q[2];
rz(-1.0591966) q[2];
sx q[2];
rz(-0.025040778) q[2];
rz(-2.5981564) q[3];
sx q[3];
rz(-2.4380324) q[3];
sx q[3];
rz(0.27633015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.339879) q[0];
sx q[0];
rz(-1.0152974) q[0];
sx q[0];
rz(0.80831084) q[0];
rz(-2.8080151) q[1];
sx q[1];
rz(-1.7671276) q[1];
sx q[1];
rz(-0.50552013) q[1];
rz(1.2407718) q[2];
sx q[2];
rz(-1.4494962) q[2];
sx q[2];
rz(-1.3818216) q[2];
rz(0.65350914) q[3];
sx q[3];
rz(-2.1353888) q[3];
sx q[3];
rz(0.63513811) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
