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
rz(5.1626212) q[0];
sx q[0];
rz(10.838312) q[0];
rz(0.36355525) q[1];
sx q[1];
rz(-1.3132562) q[1];
sx q[1];
rz(-0.29120905) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094452642) q[0];
sx q[0];
rz(-1.5908851) q[0];
sx q[0];
rz(-1.8909341) q[0];
rz(-0.12367308) q[2];
sx q[2];
rz(-1.9078662) q[2];
sx q[2];
rz(-2.6840984) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0800228) q[1];
sx q[1];
rz(-1.0008308) q[1];
sx q[1];
rz(1.391173) q[1];
rz(-pi) q[2];
rz(-2.324027) q[3];
sx q[3];
rz(-2.052784) q[3];
sx q[3];
rz(0.40460247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1353961) q[2];
sx q[2];
rz(-1.5634894) q[2];
sx q[2];
rz(-3.1175933) q[2];
rz(2.7855347) q[3];
sx q[3];
rz(-2.4672716) q[3];
sx q[3];
rz(-3.1168028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0559219) q[0];
sx q[0];
rz(-0.69127685) q[0];
sx q[0];
rz(2.5066277) q[0];
rz(2.3786646) q[1];
sx q[1];
rz(-1.4632016) q[1];
sx q[1];
rz(0.91711226) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2208495) q[0];
sx q[0];
rz(-1.5524327) q[0];
sx q[0];
rz(1.6054652) q[0];
rz(-2.2516052) q[2];
sx q[2];
rz(-1.9515079) q[2];
sx q[2];
rz(-0.057304545) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32815427) q[1];
sx q[1];
rz(-1.1309237) q[1];
sx q[1];
rz(-0.085170345) q[1];
rz(-pi) q[2];
rz(0.95083745) q[3];
sx q[3];
rz(-1.7449494) q[3];
sx q[3];
rz(-3.1319428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1941173) q[2];
sx q[2];
rz(-2.6185991) q[2];
sx q[2];
rz(0.73327649) q[2];
rz(-2.0745847) q[3];
sx q[3];
rz(-1.7206444) q[3];
sx q[3];
rz(-0.1799306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5718403) q[0];
sx q[0];
rz(-1.952992) q[0];
sx q[0];
rz(-2.6133614) q[0];
rz(0.86604467) q[1];
sx q[1];
rz(-1.3287611) q[1];
sx q[1];
rz(-0.51308647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13212261) q[0];
sx q[0];
rz(-0.44053253) q[0];
sx q[0];
rz(1.6706353) q[0];
rz(0.31103525) q[2];
sx q[2];
rz(-2.1967109) q[2];
sx q[2];
rz(-3.1300822) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5410903) q[1];
sx q[1];
rz(-1.9185431) q[1];
sx q[1];
rz(0.13558023) q[1];
x q[2];
rz(3.0894075) q[3];
sx q[3];
rz(-0.5002678) q[3];
sx q[3];
rz(-2.1310998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.036006007) q[2];
sx q[2];
rz(-2.0522108) q[2];
sx q[2];
rz(-3.1028683) q[2];
rz(-2.683908) q[3];
sx q[3];
rz(-0.80515146) q[3];
sx q[3];
rz(-1.800644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7078581) q[0];
sx q[0];
rz(-0.41446328) q[0];
sx q[0];
rz(1.1210972) q[0];
rz(-2.3838249) q[1];
sx q[1];
rz(-1.6242177) q[1];
sx q[1];
rz(2.4574492) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4829041) q[0];
sx q[0];
rz(-1.7724121) q[0];
sx q[0];
rz(-0.828142) q[0];
rz(-pi) q[1];
rz(3.1314881) q[2];
sx q[2];
rz(-1.5636383) q[2];
sx q[2];
rz(-2.7814076) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5110179) q[1];
sx q[1];
rz(-1.8749494) q[1];
sx q[1];
rz(1.263878) q[1];
rz(-pi) q[2];
rz(-0.512796) q[3];
sx q[3];
rz(-2.1053475) q[3];
sx q[3];
rz(0.057138047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4844369) q[2];
sx q[2];
rz(-1.5269205) q[2];
sx q[2];
rz(-0.26051513) q[2];
rz(-2.303463) q[3];
sx q[3];
rz(-1.9939491) q[3];
sx q[3];
rz(3.1269791) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.218006) q[0];
sx q[0];
rz(-1.3777233) q[0];
sx q[0];
rz(2.6309784) q[0];
rz(1.8922197) q[1];
sx q[1];
rz(-1.1621954) q[1];
sx q[1];
rz(1.4543264) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59856738) q[0];
sx q[0];
rz(-1.784458) q[0];
sx q[0];
rz(-1.5419307) q[0];
rz(-pi) q[1];
rz(1.5227484) q[2];
sx q[2];
rz(-1.7103875) q[2];
sx q[2];
rz(-1.0447811) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.92125398) q[1];
sx q[1];
rz(-2.5272213) q[1];
sx q[1];
rz(-0.45843108) q[1];
rz(1.4081786) q[3];
sx q[3];
rz(-1.713304) q[3];
sx q[3];
rz(0.60309013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5877567) q[2];
sx q[2];
rz(-0.38148701) q[2];
sx q[2];
rz(0.39056632) q[2];
rz(0.25911123) q[3];
sx q[3];
rz(-1.6389537) q[3];
sx q[3];
rz(-2.794877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7873586) q[0];
sx q[0];
rz(-2.7095095) q[0];
sx q[0];
rz(2.4969192) q[0];
rz(1.8395754) q[1];
sx q[1];
rz(-1.6141067) q[1];
sx q[1];
rz(0.34861809) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88953274) q[0];
sx q[0];
rz(-0.67399287) q[0];
sx q[0];
rz(-2.3537082) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5041144) q[2];
sx q[2];
rz(-0.63796746) q[2];
sx q[2];
rz(-1.6055941) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.28305123) q[1];
sx q[1];
rz(-1.1521473) q[1];
sx q[1];
rz(-0.323723) q[1];
rz(-pi) q[2];
rz(2.8556267) q[3];
sx q[3];
rz(-0.12187258) q[3];
sx q[3];
rz(-1.7274477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1031441) q[2];
sx q[2];
rz(-0.95547533) q[2];
sx q[2];
rz(0.038334282) q[2];
rz(-0.96758715) q[3];
sx q[3];
rz(-1.7318232) q[3];
sx q[3];
rz(0.35965317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2751145) q[0];
sx q[0];
rz(-0.51946467) q[0];
sx q[0];
rz(-0.073609784) q[0];
rz(-1.4069125) q[1];
sx q[1];
rz(-2.584447) q[1];
sx q[1];
rz(2.4085192) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063696472) q[0];
sx q[0];
rz(-2.1025582) q[0];
sx q[0];
rz(2.1666906) q[0];
rz(0.91821155) q[2];
sx q[2];
rz(-2.2951153) q[2];
sx q[2];
rz(2.0744963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.16790529) q[1];
sx q[1];
rz(-1.9009155) q[1];
sx q[1];
rz(-2.7658505) q[1];
rz(-pi) q[2];
rz(-2.3228775) q[3];
sx q[3];
rz(-0.33740852) q[3];
sx q[3];
rz(0.81160566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0343895) q[2];
sx q[2];
rz(-1.723039) q[2];
sx q[2];
rz(-3.1305195) q[2];
rz(-0.30558807) q[3];
sx q[3];
rz(-2.0162851) q[3];
sx q[3];
rz(-1.8886458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3333862) q[0];
sx q[0];
rz(-0.60788637) q[0];
sx q[0];
rz(-2.407684) q[0];
rz(-0.58087307) q[1];
sx q[1];
rz(-1.0092694) q[1];
sx q[1];
rz(3.0427921) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89284586) q[0];
sx q[0];
rz(-0.98631421) q[0];
sx q[0];
rz(-2.2131323) q[0];
x q[1];
rz(-2.4634741) q[2];
sx q[2];
rz(-2.0653915) q[2];
sx q[2];
rz(-1.8910318) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.1360201) q[1];
sx q[1];
rz(-2.6248049) q[1];
sx q[1];
rz(2.6173475) q[1];
x q[2];
rz(1.6510165) q[3];
sx q[3];
rz(-1.4412243) q[3];
sx q[3];
rz(3.0440273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8431479) q[2];
sx q[2];
rz(-1.8567905) q[2];
sx q[2];
rz(-1.0996381) q[2];
rz(3.1067276) q[3];
sx q[3];
rz(-2.3979082) q[3];
sx q[3];
rz(-0.59665027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027503969) q[0];
sx q[0];
rz(-1.1570258) q[0];
sx q[0];
rz(1.1771359) q[0];
rz(1.4741395) q[1];
sx q[1];
rz(-0.93282229) q[1];
sx q[1];
rz(-1.6966381) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53430492) q[0];
sx q[0];
rz(-1.8518847) q[0];
sx q[0];
rz(-2.7318153) q[0];
rz(2.3722052) q[2];
sx q[2];
rz(-1.1379998) q[2];
sx q[2];
rz(0.80608778) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18981537) q[1];
sx q[1];
rz(-2.6221402) q[1];
sx q[1];
rz(1.7096667) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8650896) q[3];
sx q[3];
rz(-0.91405896) q[3];
sx q[3];
rz(1.2402676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16704796) q[2];
sx q[2];
rz(-1.5844774) q[2];
sx q[2];
rz(-0.17189279) q[2];
rz(-0.78957549) q[3];
sx q[3];
rz(-2.004345) q[3];
sx q[3];
rz(1.8062228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0852785) q[0];
sx q[0];
rz(-2.9264937) q[0];
sx q[0];
rz(-2.6000182) q[0];
rz(1.9821292) q[1];
sx q[1];
rz(-1.8340725) q[1];
sx q[1];
rz(-2.3084739) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3726927) q[0];
sx q[0];
rz(-1.295255) q[0];
sx q[0];
rz(0.616347) q[0];
rz(-2.59899) q[2];
sx q[2];
rz(-1.869598) q[2];
sx q[2];
rz(-1.2703368) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6852452) q[1];
sx q[1];
rz(-0.73328555) q[1];
sx q[1];
rz(1.228778) q[1];
rz(-pi) q[2];
rz(-3.1396486) q[3];
sx q[3];
rz(-2.2984347) q[3];
sx q[3];
rz(-1.1465285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1200166) q[2];
sx q[2];
rz(-2.3638201) q[2];
sx q[2];
rz(-2.1982101) q[2];
rz(-2.0639482) q[3];
sx q[3];
rz(-1.5871983) q[3];
sx q[3];
rz(0.32850346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.7151466) q[0];
sx q[0];
rz(-1.0363415) q[0];
sx q[0];
rz(2.8847726) q[0];
rz(-2.9019451) q[1];
sx q[1];
rz(-2.2365204) q[1];
sx q[1];
rz(2.0047275) q[1];
rz(-0.73020936) q[2];
sx q[2];
rz(-0.84029031) q[2];
sx q[2];
rz(2.5569722) q[2];
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
