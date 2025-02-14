OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7665793) q[0];
sx q[0];
rz(-1.2393247) q[0];
sx q[0];
rz(1.0898606) q[0];
rz(1.0640979) q[1];
sx q[1];
rz(-1.8884594) q[1];
sx q[1];
rz(0.90322948) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8860461) q[0];
sx q[0];
rz(-2.1849792) q[0];
sx q[0];
rz(0.79844676) q[0];
rz(-pi) q[1];
rz(1.5888693) q[2];
sx q[2];
rz(-2.4583092) q[2];
sx q[2];
rz(0.23427134) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9574388) q[1];
sx q[1];
rz(-1.9136962) q[1];
sx q[1];
rz(2.5162906) q[1];
rz(-pi) q[2];
rz(-2.0704449) q[3];
sx q[3];
rz(-1.4649434) q[3];
sx q[3];
rz(-1.203804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4452867) q[2];
sx q[2];
rz(-3.004965) q[2];
sx q[2];
rz(-0.43329263) q[2];
rz(0.28996921) q[3];
sx q[3];
rz(-2.2765997) q[3];
sx q[3];
rz(0.32052952) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9669773) q[0];
sx q[0];
rz(-2.2791635) q[0];
sx q[0];
rz(-1.8509266) q[0];
rz(2.4621452) q[1];
sx q[1];
rz(-0.77630711) q[1];
sx q[1];
rz(-2.2490833) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7793625) q[0];
sx q[0];
rz(-2.9625872) q[0];
sx q[0];
rz(-1.711861) q[0];
rz(-pi) q[1];
rz(2.4878534) q[2];
sx q[2];
rz(-1.8934665) q[2];
sx q[2];
rz(1.8316837) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.034169056) q[1];
sx q[1];
rz(-2.4236107) q[1];
sx q[1];
rz(0.42426829) q[1];
rz(1.1970911) q[3];
sx q[3];
rz(-1.804816) q[3];
sx q[3];
rz(1.4732052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.14900011) q[2];
sx q[2];
rz(-2.0626103) q[2];
sx q[2];
rz(-1.1326257) q[2];
rz(0.66631404) q[3];
sx q[3];
rz(-1.7814813) q[3];
sx q[3];
rz(-1.2247491) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3742974) q[0];
sx q[0];
rz(-0.39300028) q[0];
sx q[0];
rz(2.1597916) q[0];
rz(0.94217316) q[1];
sx q[1];
rz(-1.3092382) q[1];
sx q[1];
rz(-0.91032496) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019429723) q[0];
sx q[0];
rz(-2.0519901) q[0];
sx q[0];
rz(1.0896171) q[0];
rz(-pi) q[1];
rz(1.6907755) q[2];
sx q[2];
rz(-1.8455122) q[2];
sx q[2];
rz(0.071823013) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1012464) q[1];
sx q[1];
rz(-1.9034428) q[1];
sx q[1];
rz(-2.8221647) q[1];
x q[2];
rz(1.0873454) q[3];
sx q[3];
rz(-0.48763613) q[3];
sx q[3];
rz(-0.82043649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.94074774) q[2];
sx q[2];
rz(-0.36483279) q[2];
sx q[2];
rz(2.9546837) q[2];
rz(-2.9777891) q[3];
sx q[3];
rz(-2.2713594) q[3];
sx q[3];
rz(-0.12266172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0201037) q[0];
sx q[0];
rz(-1.7262456) q[0];
sx q[0];
rz(-0.69766587) q[0];
rz(-2.5949219) q[1];
sx q[1];
rz(-2.8488939) q[1];
sx q[1];
rz(-1.1950511) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.228926) q[0];
sx q[0];
rz(-2.1536835) q[0];
sx q[0];
rz(0.41622644) q[0];
rz(-pi) q[1];
rz(-2.9505492) q[2];
sx q[2];
rz(-1.4378387) q[2];
sx q[2];
rz(-1.7482479) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.89151284) q[1];
sx q[1];
rz(-1.7367474) q[1];
sx q[1];
rz(3.0338698) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.324034) q[3];
sx q[3];
rz(-2.8487848) q[3];
sx q[3];
rz(2.4532401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.67808548) q[2];
sx q[2];
rz(-1.7184075) q[2];
sx q[2];
rz(0.19742337) q[2];
rz(-3.0366963) q[3];
sx q[3];
rz(-2.0320804) q[3];
sx q[3];
rz(0.088002861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.5474434) q[0];
sx q[0];
rz(-2.6056885) q[0];
sx q[0];
rz(2.1323668) q[0];
rz(0.028845305) q[1];
sx q[1];
rz(-1.4776769) q[1];
sx q[1];
rz(0.65863329) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.665425) q[0];
sx q[0];
rz(-1.1575677) q[0];
sx q[0];
rz(0.40979235) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20435996) q[2];
sx q[2];
rz(-1.0931226) q[2];
sx q[2];
rz(-0.43276873) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.62477905) q[1];
sx q[1];
rz(-2.9492231) q[1];
sx q[1];
rz(-2.6386847) q[1];
rz(-pi) q[2];
rz(-2.4464408) q[3];
sx q[3];
rz(-1.7489232) q[3];
sx q[3];
rz(0.069308829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.37904) q[2];
sx q[2];
rz(-2.5950044) q[2];
sx q[2];
rz(-0.23400447) q[2];
rz(-2.0512569) q[3];
sx q[3];
rz(-1.5666311) q[3];
sx q[3];
rz(2.0719297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91319084) q[0];
sx q[0];
rz(-0.15616067) q[0];
sx q[0];
rz(-1.8432023) q[0];
rz(0.78650728) q[1];
sx q[1];
rz(-2.2857917) q[1];
sx q[1];
rz(-1.3589121) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6494736) q[0];
sx q[0];
rz(-2.2303061) q[0];
sx q[0];
rz(1.6716206) q[0];
rz(2.4491389) q[2];
sx q[2];
rz(-2.7992749) q[2];
sx q[2];
rz(2.3794425) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.003215) q[1];
sx q[1];
rz(-1.5691461) q[1];
sx q[1];
rz(1.5695851) q[1];
rz(-pi) q[2];
rz(1.8904866) q[3];
sx q[3];
rz(-2.0123008) q[3];
sx q[3];
rz(2.6755345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4141554) q[2];
sx q[2];
rz(-1.4667908) q[2];
sx q[2];
rz(0.1725014) q[2];
rz(0.48834673) q[3];
sx q[3];
rz(-1.1427053) q[3];
sx q[3];
rz(-1.1138227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1682424) q[0];
sx q[0];
rz(-0.872648) q[0];
sx q[0];
rz(-2.8016222) q[0];
rz(0.32644692) q[1];
sx q[1];
rz(-2.4531334) q[1];
sx q[1];
rz(1.2350157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8064855) q[0];
sx q[0];
rz(-0.77101427) q[0];
sx q[0];
rz(2.5946027) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6831362) q[2];
sx q[2];
rz(-1.8962911) q[2];
sx q[2];
rz(0.95409648) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3843879) q[1];
sx q[1];
rz(-1.565462) q[1];
sx q[1];
rz(-2.6621006) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39374473) q[3];
sx q[3];
rz(-1.3662158) q[3];
sx q[3];
rz(1.0130628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1116703) q[2];
sx q[2];
rz(-0.80283529) q[2];
sx q[2];
rz(-0.55533448) q[2];
rz(-0.16767821) q[3];
sx q[3];
rz(-1.6043112) q[3];
sx q[3];
rz(2.663747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4324206) q[0];
sx q[0];
rz(-0.92996159) q[0];
sx q[0];
rz(-1.1093371) q[0];
rz(-1.1322016) q[1];
sx q[1];
rz(-2.7448476) q[1];
sx q[1];
rz(2.1999377) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96977329) q[0];
sx q[0];
rz(-2.1226774) q[0];
sx q[0];
rz(-0.18751796) q[0];
x q[1];
rz(-1.2852816) q[2];
sx q[2];
rz(-1.7196349) q[2];
sx q[2];
rz(-0.51644737) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.68384993) q[1];
sx q[1];
rz(-2.2096429) q[1];
sx q[1];
rz(-0.15775494) q[1];
x q[2];
rz(-0.48896472) q[3];
sx q[3];
rz(-0.43114907) q[3];
sx q[3];
rz(-2.2167689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38810101) q[2];
sx q[2];
rz(-1.0270303) q[2];
sx q[2];
rz(-1.6027742) q[2];
rz(-1.7704891) q[3];
sx q[3];
rz(-1.95581) q[3];
sx q[3];
rz(-3.0901618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082315363) q[0];
sx q[0];
rz(-0.61510724) q[0];
sx q[0];
rz(2.1737461) q[0];
rz(-1.8702501) q[1];
sx q[1];
rz(-1.2920734) q[1];
sx q[1];
rz(-0.06180067) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74154749) q[0];
sx q[0];
rz(-1.1843268) q[0];
sx q[0];
rz(-0.58048141) q[0];
rz(-pi) q[1];
rz(-1.7713304) q[2];
sx q[2];
rz(-1.8248744) q[2];
sx q[2];
rz(-0.24607436) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.73229746) q[1];
sx q[1];
rz(-1.4753072) q[1];
sx q[1];
rz(-0.68640253) q[1];
rz(-0.7155719) q[3];
sx q[3];
rz(-0.4288097) q[3];
sx q[3];
rz(2.9675067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9515848) q[2];
sx q[2];
rz(-0.20415674) q[2];
sx q[2];
rz(0.07587138) q[2];
rz(1.9742981) q[3];
sx q[3];
rz(-1.1018402) q[3];
sx q[3];
rz(-0.30437881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17700125) q[0];
sx q[0];
rz(-2.8622506) q[0];
sx q[0];
rz(-2.8569073) q[0];
rz(3.0668861) q[1];
sx q[1];
rz(-1.9120522) q[1];
sx q[1];
rz(-1.7237192) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3475932) q[0];
sx q[0];
rz(-0.82836223) q[0];
sx q[0];
rz(1.8338127) q[0];
x q[1];
rz(1.333513) q[2];
sx q[2];
rz(-2.5717989) q[2];
sx q[2];
rz(2.2462318) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.98445856) q[1];
sx q[1];
rz(-0.95370953) q[1];
sx q[1];
rz(-1.0544712) q[1];
rz(-1.5261493) q[3];
sx q[3];
rz(-0.64638019) q[3];
sx q[3];
rz(2.165739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8031926) q[2];
sx q[2];
rz(-1.1407547) q[2];
sx q[2];
rz(-1.8624064) q[2];
rz(-0.98215669) q[3];
sx q[3];
rz(-1.6833865) q[3];
sx q[3];
rz(-0.88472432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90047705) q[0];
sx q[0];
rz(-1.8178839) q[0];
sx q[0];
rz(-1.5771014) q[0];
rz(1.0152394) q[1];
sx q[1];
rz(-1.992234) q[1];
sx q[1];
rz(2.0043859) q[1];
rz(1.3080636) q[2];
sx q[2];
rz(-1.6752401) q[2];
sx q[2];
rz(-0.40950767) q[2];
rz(-2.1463263) q[3];
sx q[3];
rz(-2.223816) q[3];
sx q[3];
rz(1.73903) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
