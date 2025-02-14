OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3136723) q[0];
sx q[0];
rz(-1.496324) q[0];
sx q[0];
rz(-2.9287455) q[0];
rz(-1.4053474) q[1];
sx q[1];
rz(-1.5264629) q[1];
sx q[1];
rz(2.5636173) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4302464) q[0];
sx q[0];
rz(-1.0763729) q[0];
sx q[0];
rz(1.113282) q[0];
rz(-pi) q[1];
rz(-0.87066381) q[2];
sx q[2];
rz(-1.7144702) q[2];
sx q[2];
rz(-1.1149017) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4446775) q[1];
sx q[1];
rz(-1.5648876) q[1];
sx q[1];
rz(1.9828933) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7226572) q[3];
sx q[3];
rz(-1.9351601) q[3];
sx q[3];
rz(-1.0797105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2026244) q[2];
sx q[2];
rz(-1.5105931) q[2];
sx q[2];
rz(-0.6187588) q[2];
rz(2.7135811) q[3];
sx q[3];
rz(-1.908554) q[3];
sx q[3];
rz(1.7091883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0851058) q[0];
sx q[0];
rz(-0.036186941) q[0];
sx q[0];
rz(-2.4888743) q[0];
rz(-2.808049) q[1];
sx q[1];
rz(-0.8077375) q[1];
sx q[1];
rz(-2.6127167) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7031731) q[0];
sx q[0];
rz(-2.4932917) q[0];
sx q[0];
rz(2.2080883) q[0];
rz(-pi) q[1];
rz(-1.400176) q[2];
sx q[2];
rz(-0.58310027) q[2];
sx q[2];
rz(1.0379593) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6436504) q[1];
sx q[1];
rz(-1.4406246) q[1];
sx q[1];
rz(-1.4243717) q[1];
x q[2];
rz(2.3102898) q[3];
sx q[3];
rz(-1.6171494) q[3];
sx q[3];
rz(1.6024737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7876106) q[2];
sx q[2];
rz(-1.748961) q[2];
sx q[2];
rz(0.37484136) q[2];
rz(1.2104872) q[3];
sx q[3];
rz(-0.46606246) q[3];
sx q[3];
rz(-1.9207641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1061363) q[0];
sx q[0];
rz(-1.9566256) q[0];
sx q[0];
rz(-0.55738604) q[0];
rz(-2.4222971) q[1];
sx q[1];
rz(-0.44984111) q[1];
sx q[1];
rz(2.5296899) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.681088) q[0];
sx q[0];
rz(-1.9040445) q[0];
sx q[0];
rz(1.2171576) q[0];
rz(-pi) q[1];
rz(2.923392) q[2];
sx q[2];
rz(-2.2434776) q[2];
sx q[2];
rz(-0.30029485) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0590089) q[1];
sx q[1];
rz(-1.6477668) q[1];
sx q[1];
rz(-0.99748912) q[1];
x q[2];
rz(1.7541998) q[3];
sx q[3];
rz(-1.1737744) q[3];
sx q[3];
rz(-0.77172905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3376075) q[2];
sx q[2];
rz(-1.3472202) q[2];
sx q[2];
rz(2.1862629) q[2];
rz(-2.5721926) q[3];
sx q[3];
rz(-1.1207213) q[3];
sx q[3];
rz(1.3620522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3595381) q[0];
sx q[0];
rz(-0.29490092) q[0];
sx q[0];
rz(0.094757946) q[0];
rz(-1.6057728) q[1];
sx q[1];
rz(-1.1631807) q[1];
sx q[1];
rz(-0.43704978) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9511918) q[0];
sx q[0];
rz(-3.0110208) q[0];
sx q[0];
rz(2.269362) q[0];
rz(-pi) q[1];
rz(-1.7520467) q[2];
sx q[2];
rz(-2.0006575) q[2];
sx q[2];
rz(-1.4409122) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1126276) q[1];
sx q[1];
rz(-1.3092759) q[1];
sx q[1];
rz(-0.00021832988) q[1];
rz(-pi) q[2];
rz(-0.38744827) q[3];
sx q[3];
rz(-2.193938) q[3];
sx q[3];
rz(0.98549313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7729418) q[2];
sx q[2];
rz(-2.3029885) q[2];
sx q[2];
rz(1.3402026) q[2];
rz(-2.7784427) q[3];
sx q[3];
rz(-2.4996417) q[3];
sx q[3];
rz(2.6974881) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7813613) q[0];
sx q[0];
rz(-0.64557993) q[0];
sx q[0];
rz(3.0907104) q[0];
rz(0.52967349) q[1];
sx q[1];
rz(-0.61057225) q[1];
sx q[1];
rz(-3.1231336) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3878903) q[0];
sx q[0];
rz(-1.2648996) q[0];
sx q[0];
rz(1.0359156) q[0];
rz(-pi) q[1];
rz(0.71645488) q[2];
sx q[2];
rz(-1.7229417) q[2];
sx q[2];
rz(-0.97454754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5010656) q[1];
sx q[1];
rz(-0.94087761) q[1];
sx q[1];
rz(2.7940511) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6191447) q[3];
sx q[3];
rz(-1.7237437) q[3];
sx q[3];
rz(-0.82081036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9336443) q[2];
sx q[2];
rz(-1.7174145) q[2];
sx q[2];
rz(-1.2733744) q[2];
rz(0.4452855) q[3];
sx q[3];
rz(-0.71355009) q[3];
sx q[3];
rz(-2.2860897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91525602) q[0];
sx q[0];
rz(-2.7610918) q[0];
sx q[0];
rz(2.8476727) q[0];
rz(1.8982915) q[1];
sx q[1];
rz(-1.2970122) q[1];
sx q[1];
rz(-1.4909202) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4478233) q[0];
sx q[0];
rz(-2.4943879) q[0];
sx q[0];
rz(-1.291841) q[0];
x q[1];
rz(-2.8155509) q[2];
sx q[2];
rz(-1.6614104) q[2];
sx q[2];
rz(-2.2590274) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6676644) q[1];
sx q[1];
rz(-1.3753624) q[1];
sx q[1];
rz(-0.21904314) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22533234) q[3];
sx q[3];
rz(-1.1374416) q[3];
sx q[3];
rz(3.0399377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9320429) q[2];
sx q[2];
rz(-1.2806634) q[2];
sx q[2];
rz(0.10188421) q[2];
rz(-2.3515676) q[3];
sx q[3];
rz(-2.4378889) q[3];
sx q[3];
rz(1.8858645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1387966) q[0];
sx q[0];
rz(-3.0468472) q[0];
sx q[0];
rz(-0.62762678) q[0];
rz(-1.7812642) q[1];
sx q[1];
rz(-1.5831213) q[1];
sx q[1];
rz(0.21242177) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0414114) q[0];
sx q[0];
rz(-0.79754058) q[0];
sx q[0];
rz(1.0486288) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27842267) q[2];
sx q[2];
rz(-2.1489193) q[2];
sx q[2];
rz(-0.18243624) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1362939) q[1];
sx q[1];
rz(-1.7583656) q[1];
sx q[1];
rz(-2.0779528) q[1];
x q[2];
rz(2.3644629) q[3];
sx q[3];
rz(-2.5087959) q[3];
sx q[3];
rz(1.5416886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3820496) q[2];
sx q[2];
rz(-2.2727649) q[2];
sx q[2];
rz(0.36054912) q[2];
rz(-2.3709295) q[3];
sx q[3];
rz(-1.2040851) q[3];
sx q[3];
rz(-2.8785021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3537972) q[0];
sx q[0];
rz(-1.0582835) q[0];
sx q[0];
rz(-1.8428006) q[0];
rz(2.6020488) q[1];
sx q[1];
rz(-1.455247) q[1];
sx q[1];
rz(-1.4849327) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77991262) q[0];
sx q[0];
rz(-2.49278) q[0];
sx q[0];
rz(-0.54260962) q[0];
x q[1];
rz(0.46631821) q[2];
sx q[2];
rz(-0.46874415) q[2];
sx q[2];
rz(-0.14097225) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7228884) q[1];
sx q[1];
rz(-2.8378928) q[1];
sx q[1];
rz(1.610135) q[1];
rz(-2.4267107) q[3];
sx q[3];
rz(-1.2586888) q[3];
sx q[3];
rz(1.1955137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3761042) q[2];
sx q[2];
rz(-1.8643943) q[2];
sx q[2];
rz(-0.45846024) q[2];
rz(1.1485398) q[3];
sx q[3];
rz(-3.0117922) q[3];
sx q[3];
rz(-0.54500088) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.604436) q[0];
sx q[0];
rz(-0.24188365) q[0];
sx q[0];
rz(0.81445527) q[0];
rz(0.7704598) q[1];
sx q[1];
rz(-1.6043681) q[1];
sx q[1];
rz(-1.3260215) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48416177) q[0];
sx q[0];
rz(-1.7586305) q[0];
sx q[0];
rz(-1.8321048) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70416038) q[2];
sx q[2];
rz(-1.4920782) q[2];
sx q[2];
rz(-0.12264473) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.55521132) q[1];
sx q[1];
rz(-1.5813236) q[1];
sx q[1];
rz(-1.4137244) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9853346) q[3];
sx q[3];
rz(-1.7516633) q[3];
sx q[3];
rz(-0.52494088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0731547) q[2];
sx q[2];
rz(-2.1775553) q[2];
sx q[2];
rz(0.55802074) q[2];
rz(-0.26992118) q[3];
sx q[3];
rz(-2.8840265) q[3];
sx q[3];
rz(-0.71815467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.54616) q[0];
sx q[0];
rz(-1.3256925) q[0];
sx q[0];
rz(0.85920715) q[0];
rz(2.4190306) q[1];
sx q[1];
rz(-0.95046202) q[1];
sx q[1];
rz(-1.7019466) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7116657) q[0];
sx q[0];
rz(-1.5769122) q[0];
sx q[0];
rz(1.5790918) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7669575) q[2];
sx q[2];
rz(-1.3578547) q[2];
sx q[2];
rz(2.8952026) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5240092) q[1];
sx q[1];
rz(-1.7477711) q[1];
sx q[1];
rz(2.3664775) q[1];
rz(0.52442162) q[3];
sx q[3];
rz(-0.45634142) q[3];
sx q[3];
rz(3.0681036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.667111) q[2];
sx q[2];
rz(-1.2789395) q[2];
sx q[2];
rz(1.4943592) q[2];
rz(2.7254368) q[3];
sx q[3];
rz(-1.4978724) q[3];
sx q[3];
rz(0.070629899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-3.1304929) q[0];
sx q[0];
rz(-0.53642219) q[0];
sx q[0];
rz(-0.82905967) q[0];
rz(2.2566055) q[1];
sx q[1];
rz(-0.61274715) q[1];
sx q[1];
rz(-1.3428584) q[1];
rz(0.045447363) q[2];
sx q[2];
rz(-2.5927967) q[2];
sx q[2];
rz(1.1517699) q[2];
rz(2.8881843) q[3];
sx q[3];
rz(-0.80066917) q[3];
sx q[3];
rz(-1.0193326) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
