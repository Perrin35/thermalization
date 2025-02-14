OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4910645) q[0];
sx q[0];
rz(-2.129038) q[0];
sx q[0];
rz(-0.92232409) q[0];
rz(2.2946279) q[1];
sx q[1];
rz(-1.4743409) q[1];
sx q[1];
rz(2.9234731) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51061741) q[0];
sx q[0];
rz(-1.6625064) q[0];
sx q[0];
rz(2.6789078) q[0];
rz(-0.10317219) q[2];
sx q[2];
rz(-0.37984797) q[2];
sx q[2];
rz(2.0774297) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.328124) q[1];
sx q[1];
rz(-1.7901929) q[1];
sx q[1];
rz(-3.0982927) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.078212528) q[3];
sx q[3];
rz(-1.999475) q[3];
sx q[3];
rz(-3.0826867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.073079022) q[2];
sx q[2];
rz(-1.2129236) q[2];
sx q[2];
rz(0.53654137) q[2];
rz(2.1327298) q[3];
sx q[3];
rz(-0.77102414) q[3];
sx q[3];
rz(2.3276276) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5040078) q[0];
sx q[0];
rz(-1.3115839) q[0];
sx q[0];
rz(-3.1245533) q[0];
rz(-2.0783966) q[1];
sx q[1];
rz(-0.76779643) q[1];
sx q[1];
rz(2.1868736) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9794036) q[0];
sx q[0];
rz(-1.6260176) q[0];
sx q[0];
rz(2.9121141) q[0];
rz(-2.9770697) q[2];
sx q[2];
rz(-1.8662819) q[2];
sx q[2];
rz(2.4128259) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9318051) q[1];
sx q[1];
rz(-1.3452736) q[1];
sx q[1];
rz(-2.8978845) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.223079) q[3];
sx q[3];
rz(-3.0849815) q[3];
sx q[3];
rz(-2.268102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.22975989) q[2];
sx q[2];
rz(-0.91270295) q[2];
sx q[2];
rz(-2.0274053) q[2];
rz(-2.0071425) q[3];
sx q[3];
rz(-0.19616923) q[3];
sx q[3];
rz(-0.1964868) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.590362) q[0];
sx q[0];
rz(-2.1901665) q[0];
sx q[0];
rz(2.9587342) q[0];
rz(-0.081347801) q[1];
sx q[1];
rz(-0.75129879) q[1];
sx q[1];
rz(-2.523211) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8326022) q[0];
sx q[0];
rz(-1.497735) q[0];
sx q[0];
rz(-2.354524) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9419627) q[2];
sx q[2];
rz(-2.3702894) q[2];
sx q[2];
rz(2.9735801) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8152541) q[1];
sx q[1];
rz(-2.7732121) q[1];
sx q[1];
rz(-1.1756363) q[1];
x q[2];
rz(-0.86966536) q[3];
sx q[3];
rz(-2.0732862) q[3];
sx q[3];
rz(2.0521856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0028093) q[2];
sx q[2];
rz(-1.323779) q[2];
sx q[2];
rz(0.36661026) q[2];
rz(1.9256516) q[3];
sx q[3];
rz(-1.9462908) q[3];
sx q[3];
rz(2.478157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65130305) q[0];
sx q[0];
rz(-1.2325352) q[0];
sx q[0];
rz(2.41462) q[0];
rz(2.2333249) q[1];
sx q[1];
rz(-2.5636702) q[1];
sx q[1];
rz(-1.04331) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20651992) q[0];
sx q[0];
rz(-1.5653725) q[0];
sx q[0];
rz(0.66117735) q[0];
rz(-pi) q[1];
rz(0.06177549) q[2];
sx q[2];
rz(-1.1959658) q[2];
sx q[2];
rz(2.0448562) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21884313) q[1];
sx q[1];
rz(-0.51783872) q[1];
sx q[1];
rz(2.0940368) q[1];
x q[2];
rz(1.0049694) q[3];
sx q[3];
rz(-0.64098606) q[3];
sx q[3];
rz(-1.0571277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72994453) q[2];
sx q[2];
rz(-2.8204155) q[2];
sx q[2];
rz(1.8357065) q[2];
rz(0.11792396) q[3];
sx q[3];
rz(-2.6730461) q[3];
sx q[3];
rz(1.0606631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0333198) q[0];
sx q[0];
rz(-2.1554027) q[0];
sx q[0];
rz(-1.6075851) q[0];
rz(2.8458505) q[1];
sx q[1];
rz(-1.60138) q[1];
sx q[1];
rz(-1.4588413) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61304257) q[0];
sx q[0];
rz(-2.3146221) q[0];
sx q[0];
rz(3.0463329) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2590911) q[2];
sx q[2];
rz(-0.46445981) q[2];
sx q[2];
rz(-1.2683753) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0843868) q[1];
sx q[1];
rz(-2.4090892) q[1];
sx q[1];
rz(-1.8671579) q[1];
rz(2.9107575) q[3];
sx q[3];
rz(-2.1394205) q[3];
sx q[3];
rz(-0.057536803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.38833388) q[2];
sx q[2];
rz(-2.7497141) q[2];
sx q[2];
rz(-2.4957116) q[2];
rz(1.215747) q[3];
sx q[3];
rz(-1.682621) q[3];
sx q[3];
rz(1.7997883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011768613) q[0];
sx q[0];
rz(-2.3079066) q[0];
sx q[0];
rz(0.60761333) q[0];
rz(-1.9629078) q[1];
sx q[1];
rz(-1.8557529) q[1];
sx q[1];
rz(-0.36373055) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0640489) q[0];
sx q[0];
rz(-1.0564959) q[0];
sx q[0];
rz(-0.27405996) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.046437101) q[2];
sx q[2];
rz(-1.2383467) q[2];
sx q[2];
rz(-2.9118371) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.44139578) q[1];
sx q[1];
rz(-1.4508985) q[1];
sx q[1];
rz(-2.3223206) q[1];
rz(2.167114) q[3];
sx q[3];
rz(-1.5913561) q[3];
sx q[3];
rz(-1.009915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5522449) q[2];
sx q[2];
rz(-0.10422464) q[2];
sx q[2];
rz(1.9488526) q[2];
rz(-2.4097811) q[3];
sx q[3];
rz(-1.7164427) q[3];
sx q[3];
rz(1.4881136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8530387) q[0];
sx q[0];
rz(-2.9712501) q[0];
sx q[0];
rz(-1.5227675) q[0];
rz(-1.6743926) q[1];
sx q[1];
rz(-1.8312981) q[1];
sx q[1];
rz(2.4640962) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38792363) q[0];
sx q[0];
rz(-0.69783995) q[0];
sx q[0];
rz(-1.4927255) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23394312) q[2];
sx q[2];
rz(-1.3877227) q[2];
sx q[2];
rz(1.2884566) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.65645331) q[1];
sx q[1];
rz(-0.59018007) q[1];
sx q[1];
rz(2.3132669) q[1];
x q[2];
rz(-0.73245184) q[3];
sx q[3];
rz(-1.4518514) q[3];
sx q[3];
rz(2.4710187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7447394) q[2];
sx q[2];
rz(-0.91529673) q[2];
sx q[2];
rz(1.8358561) q[2];
rz(-2.8030677) q[3];
sx q[3];
rz(-2.0207696) q[3];
sx q[3];
rz(0.6849851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6029538) q[0];
sx q[0];
rz(-0.21553497) q[0];
sx q[0];
rz(-0.18390528) q[0];
rz(-1.6869847) q[1];
sx q[1];
rz(-0.55211663) q[1];
sx q[1];
rz(-2.6457381) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60431474) q[0];
sx q[0];
rz(-2.0576982) q[0];
sx q[0];
rz(-0.30211289) q[0];
rz(-pi) q[1];
rz(0.83424904) q[2];
sx q[2];
rz(-0.96820346) q[2];
sx q[2];
rz(1.955223) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.4019805) q[1];
sx q[1];
rz(-0.53369265) q[1];
sx q[1];
rz(0.81783612) q[1];
rz(-pi) q[2];
rz(0.11885507) q[3];
sx q[3];
rz(-1.0215852) q[3];
sx q[3];
rz(-2.1191747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0869861) q[2];
sx q[2];
rz(-2.2973674) q[2];
sx q[2];
rz(1.9795214) q[2];
rz(1.453513) q[3];
sx q[3];
rz(-2.1789357) q[3];
sx q[3];
rz(-1.2801722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42089713) q[0];
sx q[0];
rz(-1.1524042) q[0];
sx q[0];
rz(-2.5757117) q[0];
rz(2.7512918) q[1];
sx q[1];
rz(-1.5719599) q[1];
sx q[1];
rz(1.3077259) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85912624) q[0];
sx q[0];
rz(-2.0604366) q[0];
sx q[0];
rz(-2.7895176) q[0];
rz(0.21251596) q[2];
sx q[2];
rz(-1.2005998) q[2];
sx q[2];
rz(-0.0070564673) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.93995982) q[1];
sx q[1];
rz(-1.8324513) q[1];
sx q[1];
rz(2.5152339) q[1];
rz(1.4212178) q[3];
sx q[3];
rz(-1.4523376) q[3];
sx q[3];
rz(-1.3988023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.86964837) q[2];
sx q[2];
rz(-2.7248236) q[2];
sx q[2];
rz(1.0375674) q[2];
rz(-1.8218254) q[3];
sx q[3];
rz(-0.82873738) q[3];
sx q[3];
rz(0.64798361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8808402) q[0];
sx q[0];
rz(-0.97761959) q[0];
sx q[0];
rz(2.4993437) q[0];
rz(-1.3308659) q[1];
sx q[1];
rz(-0.86527491) q[1];
sx q[1];
rz(-2.9790402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88961381) q[0];
sx q[0];
rz(-1.9601213) q[0];
sx q[0];
rz(0.84521238) q[0];
x q[1];
rz(-2.3847975) q[2];
sx q[2];
rz(-0.50163236) q[2];
sx q[2];
rz(-1.0197786) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7127258) q[1];
sx q[1];
rz(-1.7460632) q[1];
sx q[1];
rz(2.0574942) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4914042) q[3];
sx q[3];
rz(-2.1526436) q[3];
sx q[3];
rz(-2.1285076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6315397) q[2];
sx q[2];
rz(-1.8659464) q[2];
sx q[2];
rz(-2.1007382) q[2];
rz(-0.96796525) q[3];
sx q[3];
rz(-1.0699882) q[3];
sx q[3];
rz(0.66106838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80617245) q[0];
sx q[0];
rz(-1.5652884) q[0];
sx q[0];
rz(-0.096927222) q[0];
rz(0.23282911) q[1];
sx q[1];
rz(-1.0284582) q[1];
sx q[1];
rz(-3.1340541) q[1];
rz(-2.5977124) q[2];
sx q[2];
rz(-0.82533045) q[2];
sx q[2];
rz(-2.5173204) q[2];
rz(2.4416853) q[3];
sx q[3];
rz(-1.2716765) q[3];
sx q[3];
rz(-1.8636462) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
