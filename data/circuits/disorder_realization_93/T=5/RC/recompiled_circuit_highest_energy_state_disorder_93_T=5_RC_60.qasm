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
rz(-0.45577249) q[0];
sx q[0];
rz(-1.1205641) q[0];
sx q[0];
rz(1.4135345) q[0];
rz(0.36355525) q[1];
sx q[1];
rz(4.9699291) q[1];
sx q[1];
rz(15.416754) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6585892) q[0];
sx q[0];
rz(-1.8908672) q[0];
sx q[0];
rz(-3.1204289) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2313263) q[2];
sx q[2];
rz(-1.6874773) q[2];
sx q[2];
rz(-1.9872023) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.061569842) q[1];
sx q[1];
rz(-1.0008308) q[1];
sx q[1];
rz(1.391173) q[1];
rz(0.81756567) q[3];
sx q[3];
rz(-1.0888087) q[3];
sx q[3];
rz(-0.40460247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.0061965813) q[2];
sx q[2];
rz(-1.5634894) q[2];
sx q[2];
rz(0.023999365) q[2];
rz(-2.7855347) q[3];
sx q[3];
rz(-2.4672716) q[3];
sx q[3];
rz(3.1168028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085670797) q[0];
sx q[0];
rz(-2.4503158) q[0];
sx q[0];
rz(2.5066277) q[0];
rz(-2.3786646) q[1];
sx q[1];
rz(-1.4632016) q[1];
sx q[1];
rz(-0.91711226) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35058366) q[0];
sx q[0];
rz(-1.5361333) q[0];
sx q[0];
rz(-3.123218) q[0];
rz(-pi) q[1];
rz(-0.47562648) q[2];
sx q[2];
rz(-0.94671072) q[2];
sx q[2];
rz(1.2211354) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8134384) q[1];
sx q[1];
rz(-2.0106689) q[1];
sx q[1];
rz(-0.085170345) q[1];
x q[2];
rz(-1.2767601) q[3];
sx q[3];
rz(-2.5007479) q[3];
sx q[3];
rz(-1.8185735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.94747535) q[2];
sx q[2];
rz(-0.52299356) q[2];
sx q[2];
rz(2.4083162) q[2];
rz(-1.0670079) q[3];
sx q[3];
rz(-1.7206444) q[3];
sx q[3];
rz(0.1799306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5697524) q[0];
sx q[0];
rz(-1.1886007) q[0];
sx q[0];
rz(2.6133614) q[0];
rz(-2.275548) q[1];
sx q[1];
rz(-1.3287611) q[1];
sx q[1];
rz(-0.51308647) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7932803) q[0];
sx q[0];
rz(-1.6133119) q[0];
sx q[0];
rz(2.0094064) q[0];
rz(-2.2202291) q[2];
sx q[2];
rz(-1.3201534) q[2];
sx q[2];
rz(1.373137) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9218623) q[1];
sx q[1];
rz(-2.7693536) q[1];
sx q[1];
rz(-1.2138741) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5422899) q[3];
sx q[3];
rz(-2.0703201) q[3];
sx q[3];
rz(2.0716425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.036006007) q[2];
sx q[2];
rz(-1.0893818) q[2];
sx q[2];
rz(-0.038724381) q[2];
rz(2.683908) q[3];
sx q[3];
rz(-2.3364412) q[3];
sx q[3];
rz(-1.800644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4337346) q[0];
sx q[0];
rz(-0.41446328) q[0];
sx q[0];
rz(1.1210972) q[0];
rz(2.3838249) q[1];
sx q[1];
rz(-1.6242177) q[1];
sx q[1];
rz(0.68414348) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30267495) q[0];
sx q[0];
rz(-0.7644628) q[0];
sx q[0];
rz(1.8643127) q[0];
x q[1];
rz(2.5252417) q[2];
sx q[2];
rz(-3.1292097) q[2];
sx q[2];
rz(-1.3146666) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63057478) q[1];
sx q[1];
rz(-1.2666432) q[1];
sx q[1];
rz(1.263878) q[1];
rz(-pi) q[2];
rz(-2.6287967) q[3];
sx q[3];
rz(-1.0362451) q[3];
sx q[3];
rz(0.057138047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4844369) q[2];
sx q[2];
rz(-1.5269205) q[2];
sx q[2];
rz(0.26051513) q[2];
rz(-2.303463) q[3];
sx q[3];
rz(-1.9939491) q[3];
sx q[3];
rz(3.1269791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.218006) q[0];
sx q[0];
rz(-1.3777233) q[0];
sx q[0];
rz(-0.51061428) q[0];
rz(-1.8922197) q[1];
sx q[1];
rz(-1.9793972) q[1];
sx q[1];
rz(-1.6872663) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.175486) q[0];
sx q[0];
rz(-1.5425872) q[0];
sx q[0];
rz(-2.9278446) q[0];
rz(-0.13975039) q[2];
sx q[2];
rz(-1.6183766) q[2];
sx q[2];
rz(2.6222677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2203387) q[1];
sx q[1];
rz(-2.5272213) q[1];
sx q[1];
rz(-2.6831616) q[1];
rz(-pi) q[2];
rz(2.9972059) q[3];
sx q[3];
rz(-1.4098415) q[3];
sx q[3];
rz(-0.99100366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5877567) q[2];
sx q[2];
rz(-2.7601056) q[2];
sx q[2];
rz(0.39056632) q[2];
rz(-2.8824814) q[3];
sx q[3];
rz(-1.6389537) q[3];
sx q[3];
rz(-2.794877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3542341) q[0];
sx q[0];
rz(-2.7095095) q[0];
sx q[0];
rz(0.64467347) q[0];
rz(-1.3020172) q[1];
sx q[1];
rz(-1.6141067) q[1];
sx q[1];
rz(0.34861809) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2520599) q[0];
sx q[0];
rz(-0.67399287) q[0];
sx q[0];
rz(2.3537082) q[0];
x q[1];
rz(-2.5041144) q[2];
sx q[2];
rz(-0.63796746) q[2];
sx q[2];
rz(-1.6055941) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7182838) q[1];
sx q[1];
rz(-1.2759142) q[1];
sx q[1];
rz(-1.1319833) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5362605) q[3];
sx q[3];
rz(-1.4538962) q[3];
sx q[3];
rz(1.7021321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1031441) q[2];
sx q[2];
rz(-0.95547533) q[2];
sx q[2];
rz(3.1032584) q[2];
rz(2.1740055) q[3];
sx q[3];
rz(-1.7318232) q[3];
sx q[3];
rz(0.35965317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2751145) q[0];
sx q[0];
rz(-2.622128) q[0];
sx q[0];
rz(-3.0679829) q[0];
rz(1.7346802) q[1];
sx q[1];
rz(-0.55714566) q[1];
sx q[1];
rz(0.73307347) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.276537) q[0];
sx q[0];
rz(-0.77651327) q[0];
sx q[0];
rz(-2.3797101) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.302519) q[2];
sx q[2];
rz(-1.0986549) q[2];
sx q[2];
rz(-0.034914645) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9736874) q[1];
sx q[1];
rz(-1.2406772) q[1];
sx q[1];
rz(-2.7658505) q[1];
rz(-0.81871512) q[3];
sx q[3];
rz(-0.33740852) q[3];
sx q[3];
rz(2.329987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0343895) q[2];
sx q[2];
rz(-1.723039) q[2];
sx q[2];
rz(-3.1305195) q[2];
rz(-2.8360046) q[3];
sx q[3];
rz(-2.0162851) q[3];
sx q[3];
rz(1.8886458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3333862) q[0];
sx q[0];
rz(-2.5337063) q[0];
sx q[0];
rz(2.407684) q[0];
rz(-0.58087307) q[1];
sx q[1];
rz(-1.0092694) q[1];
sx q[1];
rz(-0.098800585) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0694612) q[0];
sx q[0];
rz(-2.0939576) q[0];
sx q[0];
rz(-2.4510372) q[0];
rz(-pi) q[1];
rz(-2.4315028) q[2];
sx q[2];
rz(-2.3259893) q[2];
sx q[2];
rz(2.2888773) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.172625) q[1];
sx q[1];
rz(-1.3208814) q[1];
sx q[1];
rz(2.6843798) q[1];
x q[2];
rz(0.55136724) q[3];
sx q[3];
rz(-2.9893162) q[3];
sx q[3];
rz(2.4874529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.29844478) q[2];
sx q[2];
rz(-1.2848022) q[2];
sx q[2];
rz(-2.0419545) q[2];
rz(-0.034865033) q[3];
sx q[3];
rz(-0.74368447) q[3];
sx q[3];
rz(-2.5449424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-3.1140887) q[0];
sx q[0];
rz(-1.9845668) q[0];
sx q[0];
rz(1.1771359) q[0];
rz(1.4741395) q[1];
sx q[1];
rz(-0.93282229) q[1];
sx q[1];
rz(1.4449545) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6048353) q[0];
sx q[0];
rz(-0.49237456) q[0];
sx q[0];
rz(2.5144666) q[0];
rz(-pi) q[1];
rz(2.5553702) q[2];
sx q[2];
rz(-2.2811523) q[2];
sx q[2];
rz(1.1733871) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.260238) q[1];
sx q[1];
rz(-1.5020276) q[1];
sx q[1];
rz(-1.0555022) q[1];
rz(-pi) q[2];
rz(2.2463465) q[3];
sx q[3];
rz(-1.352868) q[3];
sx q[3];
rz(-0.15897863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9745447) q[2];
sx q[2];
rz(-1.5844774) q[2];
sx q[2];
rz(-0.17189279) q[2];
rz(-0.78957549) q[3];
sx q[3];
rz(-2.004345) q[3];
sx q[3];
rz(-1.3353698) q[3];
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
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0852785) q[0];
sx q[0];
rz(-2.9264937) q[0];
sx q[0];
rz(-0.54157448) q[0];
rz(-1.1594634) q[1];
sx q[1];
rz(-1.3075202) q[1];
sx q[1];
rz(2.3084739) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7689) q[0];
sx q[0];
rz(-1.8463377) q[0];
sx q[0];
rz(-0.616347) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9160742) q[2];
sx q[2];
rz(-1.0547027) q[2];
sx q[2];
rz(-2.6654625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1320859) q[1];
sx q[1];
rz(-0.888538) q[1];
sx q[1];
rz(-2.8481774) q[1];
rz(-3.1396486) q[3];
sx q[3];
rz(-0.84315791) q[3];
sx q[3];
rz(1.1465285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1200166) q[2];
sx q[2];
rz(-0.77777255) q[2];
sx q[2];
rz(-2.1982101) q[2];
rz(2.0639482) q[3];
sx q[3];
rz(-1.5543944) q[3];
sx q[3];
rz(-2.8130892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42644603) q[0];
sx q[0];
rz(-1.0363415) q[0];
sx q[0];
rz(2.8847726) q[0];
rz(0.23964755) q[1];
sx q[1];
rz(-2.2365204) q[1];
sx q[1];
rz(2.0047275) q[1];
rz(0.93076651) q[2];
sx q[2];
rz(-2.1590744) q[2];
sx q[2];
rz(-2.7958913) q[2];
rz(-2.4619674) q[3];
sx q[3];
rz(-2.7503895) q[3];
sx q[3];
rz(2.214307) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
