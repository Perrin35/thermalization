OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5912136) q[0];
sx q[0];
rz(-0.0033291078) q[0];
sx q[0];
rz(2.79628) q[0];
rz(-1.3357063) q[1];
sx q[1];
rz(-0.3392646) q[1];
sx q[1];
rz(0.27944922) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011443519) q[0];
sx q[0];
rz(-2.7675417) q[0];
sx q[0];
rz(-2.1610297) q[0];
rz(-pi) q[1];
rz(-1.2279534) q[2];
sx q[2];
rz(-0.81878412) q[2];
sx q[2];
rz(2.0057099) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.93912032) q[1];
sx q[1];
rz(-1.3807553) q[1];
sx q[1];
rz(0.01357667) q[1];
x q[2];
rz(-2.4996098) q[3];
sx q[3];
rz(-1.5390453) q[3];
sx q[3];
rz(2.8905405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.44148579) q[2];
sx q[2];
rz(-1.4725279) q[2];
sx q[2];
rz(-2.1842365) q[2];
rz(2.8422614) q[3];
sx q[3];
rz(-0.39756164) q[3];
sx q[3];
rz(-0.41199747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1608202) q[0];
sx q[0];
rz(-0.94962025) q[0];
sx q[0];
rz(2.7278996) q[0];
rz(-1.7970239) q[1];
sx q[1];
rz(-0.78318703) q[1];
sx q[1];
rz(2.5059674) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3755075) q[0];
sx q[0];
rz(-1.5264395) q[0];
sx q[0];
rz(-1.4384598) q[0];
rz(-pi) q[1];
rz(1.7027249) q[2];
sx q[2];
rz(-1.5243829) q[2];
sx q[2];
rz(-0.012243587) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7285068) q[1];
sx q[1];
rz(-1.8131078) q[1];
sx q[1];
rz(0.26279022) q[1];
rz(0.91931822) q[3];
sx q[3];
rz(-2.161536) q[3];
sx q[3];
rz(-0.35349333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.7522493) q[2];
sx q[2];
rz(-1.2717335) q[2];
sx q[2];
rz(-0.52865571) q[2];
rz(-1.2403437) q[3];
sx q[3];
rz(-2.7820008) q[3];
sx q[3];
rz(2.8958029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5813331) q[0];
sx q[0];
rz(-1.9868877) q[0];
sx q[0];
rz(2.8833959) q[0];
rz(-1.6351581) q[1];
sx q[1];
rz(-2.5841027) q[1];
sx q[1];
rz(-0.43446508) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2697849) q[0];
sx q[0];
rz(-1.9315533) q[0];
sx q[0];
rz(2.0478134) q[0];
rz(-3.0760879) q[2];
sx q[2];
rz(-1.0422049) q[2];
sx q[2];
rz(2.5158109) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9745969) q[1];
sx q[1];
rz(-1.3612862) q[1];
sx q[1];
rz(-2.6796883) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2779949) q[3];
sx q[3];
rz(-0.86175418) q[3];
sx q[3];
rz(1.1532702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.42052856) q[2];
sx q[2];
rz(-1.3238182) q[2];
sx q[2];
rz(3.0857962) q[2];
rz(0.93786401) q[3];
sx q[3];
rz(-0.28356975) q[3];
sx q[3];
rz(2.8985033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043512251) q[0];
sx q[0];
rz(-0.94399095) q[0];
sx q[0];
rz(3.0010624) q[0];
rz(0.17164104) q[1];
sx q[1];
rz(-1.3069897) q[1];
sx q[1];
rz(-2.8780639) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9517277) q[0];
sx q[0];
rz(-0.57706149) q[0];
sx q[0];
rz(-0.18738562) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61435917) q[2];
sx q[2];
rz(-2.376997) q[2];
sx q[2];
rz(-1.6524397) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4315223) q[1];
sx q[1];
rz(-2.353851) q[1];
sx q[1];
rz(2.6452933) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2372381) q[3];
sx q[3];
rz(-1.8309438) q[3];
sx q[3];
rz(0.0098269193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0984829) q[2];
sx q[2];
rz(-0.3442328) q[2];
sx q[2];
rz(1.8919224) q[2];
rz(1.194687) q[3];
sx q[3];
rz(-2.1134816) q[3];
sx q[3];
rz(2.4627114) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93976218) q[0];
sx q[0];
rz(-1.285346) q[0];
sx q[0];
rz(-0.029065954) q[0];
rz(-1.7395696) q[1];
sx q[1];
rz(-0.35750917) q[1];
sx q[1];
rz(-2.8129541) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5111361) q[0];
sx q[0];
rz(-1.6726613) q[0];
sx q[0];
rz(2.0589552) q[0];
x q[1];
rz(-0.26942307) q[2];
sx q[2];
rz(-1.7413483) q[2];
sx q[2];
rz(-1.5829057) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6824274) q[1];
sx q[1];
rz(-2.703856) q[1];
sx q[1];
rz(-2.9912205) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0786812) q[3];
sx q[3];
rz(-2.2552367) q[3];
sx q[3];
rz(-1.9198315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0236726) q[2];
sx q[2];
rz(-2.6866044) q[2];
sx q[2];
rz(2.9809791) q[2];
rz(1.1462071) q[3];
sx q[3];
rz(-1.3137772) q[3];
sx q[3];
rz(-2.5207991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49204957) q[0];
sx q[0];
rz(-0.8594802) q[0];
sx q[0];
rz(0.78654003) q[0];
rz(-2.7644764) q[1];
sx q[1];
rz(-0.84723324) q[1];
sx q[1];
rz(2.699111) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15765685) q[0];
sx q[0];
rz(-1.2120976) q[0];
sx q[0];
rz(1.962612) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2636678) q[2];
sx q[2];
rz(-2.1233635) q[2];
sx q[2];
rz(-0.27586684) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.46736273) q[1];
sx q[1];
rz(-1.2503887) q[1];
sx q[1];
rz(0.25735374) q[1];
rz(-pi) q[2];
rz(-1.8864185) q[3];
sx q[3];
rz(-1.481803) q[3];
sx q[3];
rz(-2.2790668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68142146) q[2];
sx q[2];
rz(-0.57839102) q[2];
sx q[2];
rz(-2.1703413) q[2];
rz(2.4462637) q[3];
sx q[3];
rz(-2.9054502) q[3];
sx q[3];
rz(-2.7142081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63715315) q[0];
sx q[0];
rz(-1.210286) q[0];
sx q[0];
rz(-1.0094281) q[0];
rz(-0.55074739) q[1];
sx q[1];
rz(-1.8661205) q[1];
sx q[1];
rz(-1.9783463) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35818415) q[0];
sx q[0];
rz(-1.9005214) q[0];
sx q[0];
rz(-0.14915906) q[0];
rz(-0.25068702) q[2];
sx q[2];
rz(-0.87826585) q[2];
sx q[2];
rz(0.78679774) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0102331) q[1];
sx q[1];
rz(-1.3512003) q[1];
sx q[1];
rz(1.4036199) q[1];
rz(-pi) q[2];
rz(0.87485119) q[3];
sx q[3];
rz(-1.8112744) q[3];
sx q[3];
rz(1.2547753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.358868) q[2];
sx q[2];
rz(-1.1512558) q[2];
sx q[2];
rz(-0.7981832) q[2];
rz(2.362137) q[3];
sx q[3];
rz(-0.53648406) q[3];
sx q[3];
rz(-1.9688169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9629795) q[0];
sx q[0];
rz(-2.6586752) q[0];
sx q[0];
rz(-3.0187507) q[0];
rz(-3.0154862) q[1];
sx q[1];
rz(-1.6364731) q[1];
sx q[1];
rz(1.2164446) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.018369) q[0];
sx q[0];
rz(-0.95802486) q[0];
sx q[0];
rz(-1.0606674) q[0];
rz(-pi) q[1];
rz(-0.33151303) q[2];
sx q[2];
rz(-2.6636332) q[2];
sx q[2];
rz(-2.8216459) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0876906) q[1];
sx q[1];
rz(-1.1864099) q[1];
sx q[1];
rz(-1.3575421) q[1];
rz(-pi) q[2];
rz(-0.27463953) q[3];
sx q[3];
rz(-2.5317149) q[3];
sx q[3];
rz(-0.21851893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.2609666) q[2];
sx q[2];
rz(-0.61769056) q[2];
sx q[2];
rz(-0.043126062) q[2];
rz(-0.17523266) q[3];
sx q[3];
rz(-0.8774811) q[3];
sx q[3];
rz(1.5568679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.50424987) q[0];
sx q[0];
rz(-3.0872587) q[0];
sx q[0];
rz(-0.20877008) q[0];
rz(1.4978706) q[1];
sx q[1];
rz(-1.6521963) q[1];
sx q[1];
rz(2.0170905) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61125206) q[0];
sx q[0];
rz(-1.7065305) q[0];
sx q[0];
rz(-1.432857) q[0];
x q[1];
rz(0.30377109) q[2];
sx q[2];
rz(-2.4035932) q[2];
sx q[2];
rz(2.8673025) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81793663) q[1];
sx q[1];
rz(-1.4374104) q[1];
sx q[1];
rz(3.1404085) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1379925) q[3];
sx q[3];
rz(-0.78046679) q[3];
sx q[3];
rz(-1.4844075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0008529) q[2];
sx q[2];
rz(-2.4177987) q[2];
sx q[2];
rz(0.17803426) q[2];
rz(1.5464276) q[3];
sx q[3];
rz(-1.833257) q[3];
sx q[3];
rz(-0.49062887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35266018) q[0];
sx q[0];
rz(-1.1505609) q[0];
sx q[0];
rz(1.0349405) q[0];
rz(2.3433698) q[1];
sx q[1];
rz(-0.98033506) q[1];
sx q[1];
rz(2.9916874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4231739) q[0];
sx q[0];
rz(-0.75773865) q[0];
sx q[0];
rz(0.54990479) q[0];
x q[1];
rz(3.0860076) q[2];
sx q[2];
rz(-0.28700799) q[2];
sx q[2];
rz(-2.3983948) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1356218) q[1];
sx q[1];
rz(-1.0711728) q[1];
sx q[1];
rz(0.47132229) q[1];
rz(-pi) q[2];
rz(0.67930119) q[3];
sx q[3];
rz(-1.3759817) q[3];
sx q[3];
rz(-0.65115813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7297111) q[2];
sx q[2];
rz(-1.231266) q[2];
sx q[2];
rz(2.7330772) q[2];
rz(-2.216693) q[3];
sx q[3];
rz(-2.3291406) q[3];
sx q[3];
rz(0.48172054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9020486) q[0];
sx q[0];
rz(-1.5938546) q[0];
sx q[0];
rz(-1.6123733) q[0];
rz(-1.7655903) q[1];
sx q[1];
rz(-1.9648432) q[1];
sx q[1];
rz(1.2479938) q[1];
rz(-0.79509905) q[2];
sx q[2];
rz(-1.8021402) q[2];
sx q[2];
rz(-2.6593593) q[2];
rz(0.33384791) q[3];
sx q[3];
rz(-2.4945989) q[3];
sx q[3];
rz(0.37024959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];