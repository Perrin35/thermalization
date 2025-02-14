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
rz(5.7972941) q[1];
sx q[1];
rz(13.426933) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33968494) q[0];
sx q[0];
rz(-0.97754495) q[0];
sx q[0];
rz(-2.9723333) q[0];
x q[1];
rz(-1.8612618) q[2];
sx q[2];
rz(-2.4142401) q[2];
sx q[2];
rz(1.0199821) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.027801188) q[1];
sx q[1];
rz(-1.6206726) q[1];
sx q[1];
rz(-2.8667843) q[1];
rz(-0.49500449) q[3];
sx q[3];
rz(-2.0873063) q[3];
sx q[3];
rz(-0.65935307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.90031558) q[2];
sx q[2];
rz(-2.6130455) q[2];
sx q[2];
rz(0.47145525) q[2];
rz(0.49247646) q[3];
sx q[3];
rz(-1.2329279) q[3];
sx q[3];
rz(1.2537778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8660698) q[0];
sx q[0];
rz(-0.38551426) q[0];
sx q[0];
rz(-0.27819124) q[0];
rz(1.9506075) q[1];
sx q[1];
rz(-1.8089801) q[1];
sx q[1];
rz(1.8461548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71466509) q[0];
sx q[0];
rz(-0.97070995) q[0];
sx q[0];
rz(0.99487181) q[0];
x q[1];
rz(-2.5668199) q[2];
sx q[2];
rz(-0.094755562) q[2];
sx q[2];
rz(-0.049959838) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5832117) q[1];
sx q[1];
rz(-1.9511327) q[1];
sx q[1];
rz(0.0024057121) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53923082) q[3];
sx q[3];
rz(-2.1279716) q[3];
sx q[3];
rz(2.2702366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1969455) q[2];
sx q[2];
rz(-1.4713919) q[2];
sx q[2];
rz(0.50986457) q[2];
rz(-1.4465796) q[3];
sx q[3];
rz(-0.46049419) q[3];
sx q[3];
rz(-2.6367326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.533621) q[0];
sx q[0];
rz(-1.5190268) q[0];
sx q[0];
rz(-0.98440379) q[0];
rz(0.016117485) q[1];
sx q[1];
rz(-2.0033483) q[1];
sx q[1];
rz(-1.791753) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72347092) q[0];
sx q[0];
rz(-1.7364565) q[0];
sx q[0];
rz(-2.4322574) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8124766) q[2];
sx q[2];
rz(-1.0306685) q[2];
sx q[2];
rz(-0.42641668) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3193839) q[1];
sx q[1];
rz(-2.2693386) q[1];
sx q[1];
rz(-1.5181659) q[1];
rz(-pi) q[2];
rz(-0.11589072) q[3];
sx q[3];
rz(-0.76911622) q[3];
sx q[3];
rz(-2.9562841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.71780378) q[2];
sx q[2];
rz(-0.70222792) q[2];
sx q[2];
rz(-0.88199893) q[2];
rz(1.9574022) q[3];
sx q[3];
rz(-2.2004674) q[3];
sx q[3];
rz(-2.4052896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0475567) q[0];
sx q[0];
rz(-1.2268257) q[0];
sx q[0];
rz(-3.0645698) q[0];
rz(0.19829622) q[1];
sx q[1];
rz(-1.501333) q[1];
sx q[1];
rz(-2.95453) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8683851) q[0];
sx q[0];
rz(-1.4131695) q[0];
sx q[0];
rz(-0.34021838) q[0];
rz(-pi) q[1];
rz(2.9603029) q[2];
sx q[2];
rz(-1.224887) q[2];
sx q[2];
rz(0.51900253) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.143972) q[1];
sx q[1];
rz(-1.629843) q[1];
sx q[1];
rz(-0.70934341) q[1];
x q[2];
rz(0.31627215) q[3];
sx q[3];
rz(-1.894884) q[3];
sx q[3];
rz(-2.9106331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5483115) q[2];
sx q[2];
rz(-0.45479861) q[2];
sx q[2];
rz(-1.6853257) q[2];
rz(2.4233387) q[3];
sx q[3];
rz(-1.7536438) q[3];
sx q[3];
rz(-0.83435241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2497571) q[0];
sx q[0];
rz(-1.8064073) q[0];
sx q[0];
rz(-2.7030113) q[0];
rz(0.35722411) q[1];
sx q[1];
rz(-1.2780739) q[1];
sx q[1];
rz(1.8908148) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8880576) q[0];
sx q[0];
rz(-1.4708637) q[0];
sx q[0];
rz(-1.9991831) q[0];
x q[1];
rz(-1.2594455) q[2];
sx q[2];
rz(-1.4964474) q[2];
sx q[2];
rz(1.6029364) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5099796) q[1];
sx q[1];
rz(-0.56029472) q[1];
sx q[1];
rz(-2.3255682) q[1];
x q[2];
rz(2.6762371) q[3];
sx q[3];
rz(-0.28959238) q[3];
sx q[3];
rz(0.48894879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9095416) q[2];
sx q[2];
rz(-0.48018685) q[2];
sx q[2];
rz(-1.5117744) q[2];
rz(0.016228598) q[3];
sx q[3];
rz(-2.0461693) q[3];
sx q[3];
rz(-0.79286638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(0.42090526) q[1];
sx q[1];
rz(-1.6890539) q[1];
sx q[1];
rz(3.0973869) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6290508) q[0];
sx q[0];
rz(-2.3568332) q[0];
sx q[0];
rz(0.46585887) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8755336) q[2];
sx q[2];
rz(-0.6630156) q[2];
sx q[2];
rz(-0.11817486) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.068983745) q[1];
sx q[1];
rz(-1.9602959) q[1];
sx q[1];
rz(-3.054148) q[1];
x q[2];
rz(-2.6085195) q[3];
sx q[3];
rz(-1.3309475) q[3];
sx q[3];
rz(-1.9727954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9629024) q[2];
sx q[2];
rz(-2.2447605) q[2];
sx q[2];
rz(0.94513354) q[2];
rz(1.4611698) q[3];
sx q[3];
rz(-1.1472568) q[3];
sx q[3];
rz(-0.51819658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1951676) q[0];
sx q[0];
rz(-2.8312046) q[0];
sx q[0];
rz(2.2921966) q[0];
rz(-1.7646344) q[1];
sx q[1];
rz(-0.57146776) q[1];
sx q[1];
rz(-1.2845385) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3234069) q[0];
sx q[0];
rz(-1.0844829) q[0];
sx q[0];
rz(-2.0165927) q[0];
rz(-2.8148267) q[2];
sx q[2];
rz(-0.32354718) q[2];
sx q[2];
rz(-1.7015333) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9640935) q[1];
sx q[1];
rz(-2.0960137) q[1];
sx q[1];
rz(-1.5296658) q[1];
rz(-pi) q[2];
rz(2.6553978) q[3];
sx q[3];
rz(-1.0665585) q[3];
sx q[3];
rz(2.274226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.66594243) q[2];
sx q[2];
rz(-1.3694171) q[2];
sx q[2];
rz(1.5671889) q[2];
rz(0.33268467) q[3];
sx q[3];
rz(-2.0761469) q[3];
sx q[3];
rz(0.98193297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57629267) q[0];
sx q[0];
rz(-1.1197634) q[0];
sx q[0];
rz(2.0530307) q[0];
rz(0.92900485) q[1];
sx q[1];
rz(-1.4938415) q[1];
sx q[1];
rz(-0.30379024) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7244686) q[0];
sx q[0];
rz(-1.9433634) q[0];
sx q[0];
rz(2.5346816) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94307282) q[2];
sx q[2];
rz(-2.94063) q[2];
sx q[2];
rz(-0.64370868) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6188083) q[1];
sx q[1];
rz(-1.4780103) q[1];
sx q[1];
rz(-2.001808) q[1];
rz(-pi) q[2];
rz(2.522842) q[3];
sx q[3];
rz(-1.9721974) q[3];
sx q[3];
rz(-1.5681745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2527689) q[2];
sx q[2];
rz(-1.6656275) q[2];
sx q[2];
rz(2.8775173) q[2];
rz(-1.1325599) q[3];
sx q[3];
rz(-2.8045636) q[3];
sx q[3];
rz(-1.0549217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7842512) q[0];
sx q[0];
rz(-2.7989474) q[0];
sx q[0];
rz(-1.0145048) q[0];
rz(-0.92813379) q[1];
sx q[1];
rz(-1.4200297) q[1];
sx q[1];
rz(-0.49044213) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7883915) q[0];
sx q[0];
rz(-1.4239131) q[0];
sx q[0];
rz(-0.066195458) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3447475) q[2];
sx q[2];
rz(-1.9354068) q[2];
sx q[2];
rz(0.62275902) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.52940581) q[1];
sx q[1];
rz(-1.2703964) q[1];
sx q[1];
rz(-3.0985188) q[1];
rz(-1.1808628) q[3];
sx q[3];
rz(-2.56757) q[3];
sx q[3];
rz(-2.5244978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.32533112) q[2];
sx q[2];
rz(-2.0897431) q[2];
sx q[2];
rz(3.0391147) q[2];
rz(-1.2517733) q[3];
sx q[3];
rz(-1.3636369) q[3];
sx q[3];
rz(-1.4832835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3935299) q[0];
sx q[0];
rz(-0.04638014) q[0];
sx q[0];
rz(-1.8512132) q[0];
rz(-2.6875467) q[1];
sx q[1];
rz(-1.3244649) q[1];
sx q[1];
rz(-0.18347278) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7843439) q[0];
sx q[0];
rz(-1.5518477) q[0];
sx q[0];
rz(-1.3878893) q[0];
rz(-pi) q[1];
rz(1.5558335) q[2];
sx q[2];
rz(-0.97285473) q[2];
sx q[2];
rz(0.19549616) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2776838) q[1];
sx q[1];
rz(-0.49278345) q[1];
sx q[1];
rz(3.0911195) q[1];
rz(1.0049263) q[3];
sx q[3];
rz(-1.4938032) q[3];
sx q[3];
rz(-2.1044097) q[3];
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
rz(0.54343623) q[3];
sx q[3];
rz(-0.70356026) q[3];
sx q[3];
rz(-0.27633015) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.339879) q[0];
sx q[0];
rz(-2.1262953) q[0];
sx q[0];
rz(-2.3332818) q[0];
rz(-0.33357757) q[1];
sx q[1];
rz(-1.3744651) q[1];
sx q[1];
rz(2.6360725) q[1];
rz(1.2110151) q[2];
sx q[2];
rz(-0.3508437) q[2];
sx q[2];
rz(2.99101) q[2];
rz(0.80584851) q[3];
sx q[3];
rz(-2.3060006) q[3];
sx q[3];
rz(-1.5455442) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
