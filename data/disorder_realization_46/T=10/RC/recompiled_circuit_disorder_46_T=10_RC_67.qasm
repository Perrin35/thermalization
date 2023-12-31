OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5503791) q[0];
sx q[0];
rz(3.1449218) q[0];
sx q[0];
rz(9.7700906) q[0];
rz(1.8058864) q[1];
sx q[1];
rz(3.4808573) q[1];
sx q[1];
rz(9.1453287) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0246668) q[0];
sx q[0];
rz(-1.3660087) q[0];
sx q[0];
rz(-1.2555528) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7819671) q[2];
sx q[2];
rz(-1.8188393) q[2];
sx q[2];
rz(2.9458407) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2742259) q[1];
sx q[1];
rz(-0.19051954) q[1];
sx q[1];
rz(-1.5003367) q[1];
rz(1.5311601) q[3];
sx q[3];
rz(-2.2124024) q[3];
sx q[3];
rz(-1.2960145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.44148579) q[2];
sx q[2];
rz(-1.6690648) q[2];
sx q[2];
rz(-2.1842365) q[2];
rz(-2.8422614) q[3];
sx q[3];
rz(-2.744031) q[3];
sx q[3];
rz(-0.41199747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1608202) q[0];
sx q[0];
rz(-2.1919724) q[0];
sx q[0];
rz(0.41369307) q[0];
rz(1.7970239) q[1];
sx q[1];
rz(-2.3584056) q[1];
sx q[1];
rz(-0.63562524) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5168415) q[0];
sx q[0];
rz(-0.13953129) q[0];
sx q[0];
rz(-1.8952888) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.046819709) q[2];
sx q[2];
rz(-1.702582) q[2];
sx q[2];
rz(-1.552396) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2555274) q[1];
sx q[1];
rz(-0.35554245) q[1];
sx q[1];
rz(-0.76053263) q[1];
rz(-pi) q[2];
rz(2.440968) q[3];
sx q[3];
rz(-1.0430338) q[3];
sx q[3];
rz(0.81567314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3893434) q[2];
sx q[2];
rz(-1.8698591) q[2];
sx q[2];
rz(-0.52865571) q[2];
rz(1.901249) q[3];
sx q[3];
rz(-2.7820008) q[3];
sx q[3];
rz(2.8958029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56025958) q[0];
sx q[0];
rz(-1.154705) q[0];
sx q[0];
rz(-0.2581968) q[0];
rz(-1.5064346) q[1];
sx q[1];
rz(-0.55748993) q[1];
sx q[1];
rz(2.7071276) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2697849) q[0];
sx q[0];
rz(-1.9315533) q[0];
sx q[0];
rz(2.0478134) q[0];
rz(-1.6824109) q[2];
sx q[2];
rz(-2.6093405) q[2];
sx q[2];
rz(0.75512952) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3420521) q[1];
sx q[1];
rz(-2.6375348) q[1];
sx q[1];
rz(0.44517681) q[1];
rz(-pi) q[2];
rz(2.4934019) q[3];
sx q[3];
rz(-2.18581) q[3];
sx q[3];
rz(2.0730413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7210641) q[2];
sx q[2];
rz(-1.8177744) q[2];
sx q[2];
rz(-3.0857962) q[2];
rz(0.93786401) q[3];
sx q[3];
rz(-0.28356975) q[3];
sx q[3];
rz(2.8985033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0980804) q[0];
sx q[0];
rz(-2.1976017) q[0];
sx q[0];
rz(0.14053024) q[0];
rz(2.9699516) q[1];
sx q[1];
rz(-1.834603) q[1];
sx q[1];
rz(0.26352873) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9517277) q[0];
sx q[0];
rz(-0.57706149) q[0];
sx q[0];
rz(-2.954207) q[0];
rz(-pi) q[1];
rz(-2.5272335) q[2];
sx q[2];
rz(-0.76459568) q[2];
sx q[2];
rz(1.6524397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7100704) q[1];
sx q[1];
rz(-2.353851) q[1];
sx q[1];
rz(2.6452933) q[1];
x q[2];
rz(2.2534091) q[3];
sx q[3];
rz(-2.7215951) q[3];
sx q[3];
rz(-2.1995467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0984829) q[2];
sx q[2];
rz(-2.7973599) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2018305) q[0];
sx q[0];
rz(-1.285346) q[0];
sx q[0];
rz(-3.1125267) q[0];
rz(1.7395696) q[1];
sx q[1];
rz(-2.7840835) q[1];
sx q[1];
rz(-2.8129541) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12954457) q[0];
sx q[0];
rz(-2.64376) q[0];
sx q[0];
rz(1.356202) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57428898) q[2];
sx q[2];
rz(-2.8238378) q[2];
sx q[2];
rz(2.6025835) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29341104) q[1];
sx q[1];
rz(-1.1383346) q[1];
sx q[1];
rz(1.5007988) q[1];
x q[2];
rz(-0.7469437) q[3];
sx q[3];
rz(-1.1960104) q[3];
sx q[3];
rz(-0.67583109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0236726) q[2];
sx q[2];
rz(-2.6866044) q[2];
sx q[2];
rz(0.16061352) q[2];
rz(1.1462071) q[3];
sx q[3];
rz(-1.8278154) q[3];
sx q[3];
rz(-0.62079352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49204957) q[0];
sx q[0];
rz(-2.2821125) q[0];
sx q[0];
rz(-2.3550526) q[0];
rz(-2.7644764) q[1];
sx q[1];
rz(-2.2943594) q[1];
sx q[1];
rz(0.4424817) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15765685) q[0];
sx q[0];
rz(-1.929495) q[0];
sx q[0];
rz(1.962612) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8779249) q[2];
sx q[2];
rz(-1.0182292) q[2];
sx q[2];
rz(-0.27586684) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.22873951) q[1];
sx q[1];
rz(-2.733426) q[1];
sx q[1];
rz(-2.2250882) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2551742) q[3];
sx q[3];
rz(-1.6597897) q[3];
sx q[3];
rz(-2.2790668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4601712) q[2];
sx q[2];
rz(-0.57839102) q[2];
sx q[2];
rz(2.1703413) q[2];
rz(-0.69532895) q[3];
sx q[3];
rz(-2.9054502) q[3];
sx q[3];
rz(-2.7142081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5044395) q[0];
sx q[0];
rz(-1.9313066) q[0];
sx q[0];
rz(-2.1321645) q[0];
rz(-2.5908453) q[1];
sx q[1];
rz(-1.2754722) q[1];
sx q[1];
rz(-1.9783463) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.163994) q[0];
sx q[0];
rz(-1.4297276) q[0];
sx q[0];
rz(1.903957) q[0];
rz(-pi) q[1];
rz(-2.2789754) q[2];
sx q[2];
rz(-1.3786945) q[2];
sx q[2];
rz(-2.5196599) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5442859) q[1];
sx q[1];
rz(-1.7339216) q[1];
sx q[1];
rz(0.2225999) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30927741) q[3];
sx q[3];
rz(-0.89865548) q[3];
sx q[3];
rz(2.6291763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.78272468) q[2];
sx q[2];
rz(-1.9903368) q[2];
sx q[2];
rz(-2.3434095) q[2];
rz(-0.77945566) q[3];
sx q[3];
rz(-2.6051086) q[3];
sx q[3];
rz(1.9688169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9629795) q[0];
sx q[0];
rz(-0.48291746) q[0];
sx q[0];
rz(3.0187507) q[0];
rz(3.0154862) q[1];
sx q[1];
rz(-1.5051196) q[1];
sx q[1];
rz(1.2164446) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2467263) q[0];
sx q[0];
rz(-0.77573949) q[0];
sx q[0];
rz(0.6070444) q[0];
x q[1];
rz(2.8100796) q[2];
sx q[2];
rz(-0.47795948) q[2];
sx q[2];
rz(2.8216459) q[2];
rz(-pi/2) q[3];
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
rz(2.8669531) q[3];
sx q[3];
rz(-2.5317149) q[3];
sx q[3];
rz(2.9230737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8806261) q[2];
sx q[2];
rz(-0.61769056) q[2];
sx q[2];
rz(0.043126062) q[2];
rz(2.96636) q[3];
sx q[3];
rz(-0.8774811) q[3];
sx q[3];
rz(-1.5847248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50424987) q[0];
sx q[0];
rz(-3.0872587) q[0];
sx q[0];
rz(-0.20877008) q[0];
rz(-1.6437221) q[1];
sx q[1];
rz(-1.4893963) q[1];
sx q[1];
rz(-2.0170905) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9546684) q[0];
sx q[0];
rz(-0.19321975) q[0];
sx q[0];
rz(-0.78878553) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71473177) q[2];
sx q[2];
rz(-1.7734314) q[2];
sx q[2];
rz(-1.5243901) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.323656) q[1];
sx q[1];
rz(-1.4374104) q[1];
sx q[1];
rz(3.1404085) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1379925) q[3];
sx q[3];
rz(-0.78046679) q[3];
sx q[3];
rz(-1.6571852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1407397) q[2];
sx q[2];
rz(-2.4177987) q[2];
sx q[2];
rz(-2.9635584) q[2];
rz(-1.5464276) q[3];
sx q[3];
rz(-1.833257) q[3];
sx q[3];
rz(-2.6509638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35266018) q[0];
sx q[0];
rz(-1.1505609) q[0];
sx q[0];
rz(-2.1066522) q[0];
rz(0.79822284) q[1];
sx q[1];
rz(-2.1612576) q[1];
sx q[1];
rz(-0.14990526) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56652727) q[0];
sx q[0];
rz(-1.9381822) q[0];
sx q[0];
rz(-0.67879403) q[0];
rz(-pi) q[1];
rz(2.855004) q[2];
sx q[2];
rz(-1.5865241) q[2];
sx q[2];
rz(0.77428267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18969892) q[1];
sx q[1];
rz(-0.67283291) q[1];
sx q[1];
rz(-2.2646907) q[1];
x q[2];
rz(0.30431872) q[3];
sx q[3];
rz(-2.439194) q[3];
sx q[3];
rz(1.9866634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4118816) q[2];
sx q[2];
rz(-1.9103266) q[2];
sx q[2];
rz(-0.40851545) q[2];
rz(2.216693) q[3];
sx q[3];
rz(-0.81245208) q[3];
sx q[3];
rz(-2.6598721) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2395441) q[0];
sx q[0];
rz(-1.5938546) q[0];
sx q[0];
rz(-1.6123733) q[0];
rz(-1.3760024) q[1];
sx q[1];
rz(-1.1767495) q[1];
sx q[1];
rz(-1.8935988) q[1];
rz(2.8228912) q[2];
sx q[2];
rz(-2.3206884) q[2];
sx q[2];
rz(-0.86736292) q[2];
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
