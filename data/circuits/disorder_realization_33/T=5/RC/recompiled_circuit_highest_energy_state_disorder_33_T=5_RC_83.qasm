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
rz(1.1627816) q[0];
sx q[0];
rz(-0.94036189) q[0];
sx q[0];
rz(-2.9475687) q[0];
rz(-3.6858978) q[1];
sx q[1];
rz(1.5679918) q[1];
sx q[1];
rz(12.369649) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08840522) q[0];
sx q[0];
rz(-2.1413099) q[0];
sx q[0];
rz(-0.14630099) q[0];
x q[1];
rz(-1.5926378) q[2];
sx q[2];
rz(-1.7882203) q[2];
sx q[2];
rz(-2.3453051) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0086872) q[1];
sx q[1];
rz(-2.7186435) q[1];
sx q[1];
rz(-2.4407331) q[1];
rz(-2.1473925) q[3];
sx q[3];
rz(-1.4399035) q[3];
sx q[3];
rz(-1.5792221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0611614) q[2];
sx q[2];
rz(-1.4208527) q[2];
sx q[2];
rz(3.1283992) q[2];
rz(2.0773928) q[3];
sx q[3];
rz(-1.0366169) q[3];
sx q[3];
rz(0.17087759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6212293) q[0];
sx q[0];
rz(-0.55084387) q[0];
sx q[0];
rz(-2.6589822) q[0];
rz(-2.6699325) q[1];
sx q[1];
rz(-0.99716798) q[1];
sx q[1];
rz(-2.4536123) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2575098) q[0];
sx q[0];
rz(-0.36195746) q[0];
sx q[0];
rz(-2.5779526) q[0];
x q[1];
rz(-1.2763763) q[2];
sx q[2];
rz(-0.79508077) q[2];
sx q[2];
rz(2.9816422) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6370442) q[1];
sx q[1];
rz(-2.6682862) q[1];
sx q[1];
rz(3.0454703) q[1];
rz(-pi) q[2];
x q[2];
rz(0.08556871) q[3];
sx q[3];
rz(-1.7734756) q[3];
sx q[3];
rz(2.6406555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1241577) q[2];
sx q[2];
rz(-1.0729125) q[2];
sx q[2];
rz(-2.9166481) q[2];
rz(-2.0624835) q[3];
sx q[3];
rz(-2.2033117) q[3];
sx q[3];
rz(-0.45043954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17063046) q[0];
sx q[0];
rz(-2.5522975) q[0];
sx q[0];
rz(-1.851409) q[0];
rz(-1.8633441) q[1];
sx q[1];
rz(-1.9854913) q[1];
sx q[1];
rz(-1.8589171) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0874003) q[0];
sx q[0];
rz(-0.71119672) q[0];
sx q[0];
rz(1.0048701) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.01802) q[2];
sx q[2];
rz(-1.6500435) q[2];
sx q[2];
rz(1.0956956) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0209109) q[1];
sx q[1];
rz(-1.5889165) q[1];
sx q[1];
rz(-1.5906035) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77529737) q[3];
sx q[3];
rz(-2.2140045) q[3];
sx q[3];
rz(0.51706367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.67905417) q[2];
sx q[2];
rz(-1.5306229) q[2];
sx q[2];
rz(-0.47492096) q[2];
rz(1.9654407) q[3];
sx q[3];
rz(-0.85634309) q[3];
sx q[3];
rz(-0.28158751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7684286) q[0];
sx q[0];
rz(-1.7730862) q[0];
sx q[0];
rz(2.9737293) q[0];
rz(0.31790512) q[1];
sx q[1];
rz(-1.2891506) q[1];
sx q[1];
rz(1.0344523) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8694591) q[0];
sx q[0];
rz(-1.3087166) q[0];
sx q[0];
rz(-3.1054405) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9475161) q[2];
sx q[2];
rz(-2.2206306) q[2];
sx q[2];
rz(-1.5776933) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.36970475) q[1];
sx q[1];
rz(-2.2834816) q[1];
sx q[1];
rz(-0.94526498) q[1];
rz(-0.61180964) q[3];
sx q[3];
rz(-1.4880848) q[3];
sx q[3];
rz(-1.8066607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0340524) q[2];
sx q[2];
rz(-0.62959051) q[2];
sx q[2];
rz(-2.8752987) q[2];
rz(-0.1693503) q[3];
sx q[3];
rz(-2.0784056) q[3];
sx q[3];
rz(-2.9625986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87194815) q[0];
sx q[0];
rz(-0.18505159) q[0];
sx q[0];
rz(-1.2931152) q[0];
rz(-3.0710908) q[1];
sx q[1];
rz(-0.87424707) q[1];
sx q[1];
rz(2.3354796) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7775775) q[0];
sx q[0];
rz(-2.8991472) q[0];
sx q[0];
rz(1.0274946) q[0];
rz(0.057439645) q[2];
sx q[2];
rz(-1.4220793) q[2];
sx q[2];
rz(-1.1401389) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8580839) q[1];
sx q[1];
rz(-1.3779989) q[1];
sx q[1];
rz(-0.60528614) q[1];
rz(-pi) q[2];
x q[2];
rz(1.351089) q[3];
sx q[3];
rz(-1.0423805) q[3];
sx q[3];
rz(2.3912969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.21141323) q[2];
sx q[2];
rz(-1.4750007) q[2];
sx q[2];
rz(-0.30964568) q[2];
rz(-0.51112255) q[3];
sx q[3];
rz(-2.5680254) q[3];
sx q[3];
rz(-0.82093325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9822134) q[0];
sx q[0];
rz(-2.5143304) q[0];
sx q[0];
rz(-0.43701592) q[0];
rz(-0.45477319) q[1];
sx q[1];
rz(-0.79650703) q[1];
sx q[1];
rz(-2.2589267) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012271492) q[0];
sx q[0];
rz(-2.8573749) q[0];
sx q[0];
rz(1.2491262) q[0];
x q[1];
rz(-1.6505372) q[2];
sx q[2];
rz(-0.3796595) q[2];
sx q[2];
rz(-0.14282957) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3040124) q[1];
sx q[1];
rz(-1.092957) q[1];
sx q[1];
rz(0.3752075) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7035376) q[3];
sx q[3];
rz(-1.0133303) q[3];
sx q[3];
rz(-1.5749951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5060045) q[2];
sx q[2];
rz(-1.7134066) q[2];
sx q[2];
rz(-0.44090718) q[2];
rz(2.0697557) q[3];
sx q[3];
rz(-0.25522885) q[3];
sx q[3];
rz(0.60060874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1275682) q[0];
sx q[0];
rz(-1.9381645) q[0];
sx q[0];
rz(-1.2822275) q[0];
rz(-2.0558689) q[1];
sx q[1];
rz(-1.6174569) q[1];
sx q[1];
rz(-1.5132743) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8305539) q[0];
sx q[0];
rz(-0.74493248) q[0];
sx q[0];
rz(5*pi/9) q[0];
rz(-pi) q[1];
rz(2.4413928) q[2];
sx q[2];
rz(-1.8196897) q[2];
sx q[2];
rz(2.2867416) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8095743) q[1];
sx q[1];
rz(-0.93832131) q[1];
sx q[1];
rz(0.47398403) q[1];
rz(-pi) q[2];
rz(1.9508771) q[3];
sx q[3];
rz(-0.8506368) q[3];
sx q[3];
rz(2.5270568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6699803) q[2];
sx q[2];
rz(-2.4358304) q[2];
sx q[2];
rz(2.3213049) q[2];
rz(-0.39508501) q[3];
sx q[3];
rz(-2.3494651) q[3];
sx q[3];
rz(0.61409605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.87913269) q[0];
sx q[0];
rz(-1.4985871) q[0];
sx q[0];
rz(0.60687989) q[0];
rz(-2.795769) q[1];
sx q[1];
rz(-1.9551829) q[1];
sx q[1];
rz(0.80723673) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2169289) q[0];
sx q[0];
rz(-1.47424) q[0];
sx q[0];
rz(1.2256304) q[0];
rz(-1.973602) q[2];
sx q[2];
rz(-1.323408) q[2];
sx q[2];
rz(1.1198662) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8816384) q[1];
sx q[1];
rz(-1.7865766) q[1];
sx q[1];
rz(-0.27387932) q[1];
x q[2];
rz(-2.5446078) q[3];
sx q[3];
rz(-1.0506949) q[3];
sx q[3];
rz(0.4515243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3240933) q[2];
sx q[2];
rz(-1.1439415) q[2];
sx q[2];
rz(2.6452046) q[2];
rz(0.084970623) q[3];
sx q[3];
rz(-2.7111263) q[3];
sx q[3];
rz(0.5808723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9896511) q[0];
sx q[0];
rz(-1.713151) q[0];
sx q[0];
rz(0.38452837) q[0];
rz(-2.8893068) q[1];
sx q[1];
rz(-1.3832046) q[1];
sx q[1];
rz(1.5581473) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66551946) q[0];
sx q[0];
rz(-1.5613149) q[0];
sx q[0];
rz(-0.46561636) q[0];
rz(1.6884638) q[2];
sx q[2];
rz(-1.298438) q[2];
sx q[2];
rz(-0.54335574) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0817683) q[1];
sx q[1];
rz(-1.5579954) q[1];
sx q[1];
rz(-0.98155419) q[1];
x q[2];
rz(2.4676314) q[3];
sx q[3];
rz(-2.1836851) q[3];
sx q[3];
rz(-0.4919258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.79534507) q[2];
sx q[2];
rz(-1.2582422) q[2];
sx q[2];
rz(-1.1400878) q[2];
rz(1.7831066) q[3];
sx q[3];
rz(-1.2525109) q[3];
sx q[3];
rz(-2.8235249) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7110905) q[0];
sx q[0];
rz(-1.6999812) q[0];
sx q[0];
rz(2.7464113) q[0];
rz(1.2218062) q[1];
sx q[1];
rz(-0.8539353) q[1];
sx q[1];
rz(-2.2538259) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2632653) q[0];
sx q[0];
rz(-2.3774231) q[0];
sx q[0];
rz(0.40212888) q[0];
rz(-pi) q[1];
rz(1.0958065) q[2];
sx q[2];
rz(-1.5676427) q[2];
sx q[2];
rz(-2.9649865) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5612167) q[1];
sx q[1];
rz(-1.9264915) q[1];
sx q[1];
rz(0.027173398) q[1];
rz(-1.8923081) q[3];
sx q[3];
rz(-2.8148533) q[3];
sx q[3];
rz(-2.660291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6937574) q[2];
sx q[2];
rz(-1.3112023) q[2];
sx q[2];
rz(-2.5698404) q[2];
rz(-1.5104431) q[3];
sx q[3];
rz(-1.4331199) q[3];
sx q[3];
rz(1.3408287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1258662) q[0];
sx q[0];
rz(-1.3073574) q[0];
sx q[0];
rz(0.14412185) q[0];
rz(1.1852599) q[1];
sx q[1];
rz(-2.2846501) q[1];
sx q[1];
rz(1.2774998) q[1];
rz(1.324426) q[2];
sx q[2];
rz(-1.9472716) q[2];
sx q[2];
rz(-1.7903259) q[2];
rz(-0.64438728) q[3];
sx q[3];
rz(-1.3105262) q[3];
sx q[3];
rz(-2.5986828) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
