OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.65052819) q[0];
sx q[0];
rz(-1.0125546) q[0];
sx q[0];
rz(0.92232409) q[0];
rz(2.2946279) q[1];
sx q[1];
rz(-1.4743409) q[1];
sx q[1];
rz(-0.21811952) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87859652) q[0];
sx q[0];
rz(-2.6705526) q[0];
sx q[0];
rz(-2.9383927) q[0];
rz(3.0384205) q[2];
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
rz(-2.328124) q[1];
sx q[1];
rz(-1.3513997) q[1];
sx q[1];
rz(3.0982927) q[1];
rz(-0.078212528) q[3];
sx q[3];
rz(-1.999475) q[3];
sx q[3];
rz(-3.0826867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.073079022) q[2];
sx q[2];
rz(-1.928669) q[2];
sx q[2];
rz(-2.6050513) q[2];
rz(-1.0088629) q[3];
sx q[3];
rz(-2.3705685) q[3];
sx q[3];
rz(-2.3276276) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5040078) q[0];
sx q[0];
rz(-1.8300087) q[0];
sx q[0];
rz(-3.1245533) q[0];
rz(1.0631961) q[1];
sx q[1];
rz(-0.76779643) q[1];
sx q[1];
rz(-0.95471901) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5009618) q[0];
sx q[0];
rz(-0.23591536) q[0];
sx q[0];
rz(-2.9032034) q[0];
x q[1];
rz(0.164523) q[2];
sx q[2];
rz(-1.2753107) q[2];
sx q[2];
rz(0.72876677) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9318051) q[1];
sx q[1];
rz(-1.3452736) q[1];
sx q[1];
rz(0.24370812) q[1];
x q[2];
rz(-1.6240261) q[3];
sx q[3];
rz(-1.5900776) q[3];
sx q[3];
rz(0.35010168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9118328) q[2];
sx q[2];
rz(-2.2288897) q[2];
sx q[2];
rz(2.0274053) q[2];
rz(-1.1344502) q[3];
sx q[3];
rz(-2.9454234) q[3];
sx q[3];
rz(-0.1964868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.5512307) q[0];
sx q[0];
rz(-2.1901665) q[0];
sx q[0];
rz(2.9587342) q[0];
rz(0.081347801) q[1];
sx q[1];
rz(-2.3902939) q[1];
sx q[1];
rz(0.61838165) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9524117) q[0];
sx q[0];
rz(-2.3518686) q[0];
sx q[0];
rz(-0.10297063) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76128086) q[2];
sx q[2];
rz(-1.4321186) q[2];
sx q[2];
rz(-1.5468462) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.8733682) q[1];
sx q[1];
rz(-1.4317239) q[1];
sx q[1];
rz(-1.2285608) q[1];
x q[2];
rz(0.86966536) q[3];
sx q[3];
rz(-1.0683064) q[3];
sx q[3];
rz(-1.089407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1387834) q[2];
sx q[2];
rz(-1.8178136) q[2];
sx q[2];
rz(0.36661026) q[2];
rz(1.9256516) q[3];
sx q[3];
rz(-1.1953019) q[3];
sx q[3];
rz(-2.478157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4902896) q[0];
sx q[0];
rz(-1.2325352) q[0];
sx q[0];
rz(2.41462) q[0];
rz(0.90826774) q[1];
sx q[1];
rz(-2.5636702) q[1];
sx q[1];
rz(1.04331) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9350727) q[0];
sx q[0];
rz(-1.5653725) q[0];
sx q[0];
rz(-0.66117735) q[0];
rz(-pi) q[1];
rz(-1.1953148) q[2];
sx q[2];
rz(-1.6282778) q[2];
sx q[2];
rz(0.45141803) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.21884313) q[1];
sx q[1];
rz(-0.51783872) q[1];
sx q[1];
rz(-2.0940368) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0049694) q[3];
sx q[3];
rz(-2.5006066) q[3];
sx q[3];
rz(2.084465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4116481) q[2];
sx q[2];
rz(-2.8204155) q[2];
sx q[2];
rz(-1.3058861) q[2];
rz(0.11792396) q[3];
sx q[3];
rz(-0.4685466) q[3];
sx q[3];
rz(2.0809295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1082728) q[0];
sx q[0];
rz(-2.1554027) q[0];
sx q[0];
rz(1.6075851) q[0];
rz(-0.29574212) q[1];
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
rz(-0.47285493) q[0];
sx q[0];
rz(-0.74873996) q[0];
sx q[0];
rz(1.6738026) q[0];
rz(-pi) q[1];
rz(2.2590911) q[2];
sx q[2];
rz(-2.6771328) q[2];
sx q[2];
rz(1.2683753) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4319358) q[1];
sx q[1];
rz(-1.767358) q[1];
sx q[1];
rz(0.86039575) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1516861) q[3];
sx q[3];
rz(-1.3767929) q[3];
sx q[3];
rz(-1.7542183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.38833388) q[2];
sx q[2];
rz(-2.7497141) q[2];
sx q[2];
rz(0.64588109) q[2];
rz(1.215747) q[3];
sx q[3];
rz(-1.4589717) q[3];
sx q[3];
rz(1.3418044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.011768613) q[0];
sx q[0];
rz(-2.3079066) q[0];
sx q[0];
rz(-0.60761333) q[0];
rz(1.9629078) q[1];
sx q[1];
rz(-1.8557529) q[1];
sx q[1];
rz(-2.7778621) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0640489) q[0];
sx q[0];
rz(-2.0850967) q[0];
sx q[0];
rz(-2.8675327) q[0];
x q[1];
rz(-3.0951556) q[2];
sx q[2];
rz(-1.2383467) q[2];
sx q[2];
rz(2.9118371) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8848829) q[1];
sx q[1];
rz(-2.3824168) q[1];
sx q[1];
rz(1.745454) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.167114) q[3];
sx q[3];
rz(-1.5913561) q[3];
sx q[3];
rz(1.009915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5522449) q[2];
sx q[2];
rz(-3.037368) q[2];
sx q[2];
rz(1.9488526) q[2];
rz(-0.73181152) q[3];
sx q[3];
rz(-1.4251499) q[3];
sx q[3];
rz(-1.653479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28855395) q[0];
sx q[0];
rz(-0.1703425) q[0];
sx q[0];
rz(1.6188251) q[0];
rz(-1.4672) q[1];
sx q[1];
rz(-1.3102945) q[1];
sx q[1];
rz(-0.67749643) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2427432) q[0];
sx q[0];
rz(-1.5206608) q[0];
sx q[0];
rz(-0.87445796) q[0];
rz(-pi) q[1];
rz(1.3827146) q[2];
sx q[2];
rz(-1.3408337) q[2];
sx q[2];
rz(-0.32569685) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4851393) q[1];
sx q[1];
rz(-0.59018007) q[1];
sx q[1];
rz(-2.3132669) q[1];
x q[2];
rz(1.7301637) q[3];
sx q[3];
rz(-2.2969118) q[3];
sx q[3];
rz(2.1350525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39685321) q[2];
sx q[2];
rz(-0.91529673) q[2];
sx q[2];
rz(-1.8358561) q[2];
rz(2.8030677) q[3];
sx q[3];
rz(-1.1208231) q[3];
sx q[3];
rz(0.6849851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6029538) q[0];
sx q[0];
rz(-0.21553497) q[0];
sx q[0];
rz(-2.9576874) q[0];
rz(-1.454608) q[1];
sx q[1];
rz(-0.55211663) q[1];
sx q[1];
rz(-0.49585453) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9496807) q[0];
sx q[0];
rz(-2.5750468) q[0];
sx q[0];
rz(-1.0587771) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77347254) q[2];
sx q[2];
rz(-0.91433734) q[2];
sx q[2];
rz(-2.9673849) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.293893) q[1];
sx q[1];
rz(-1.2155079) q[1];
sx q[1];
rz(-1.1637264) q[1];
rz(-pi) q[2];
rz(-1.3794231) q[3];
sx q[3];
rz(-2.5809605) q[3];
sx q[3];
rz(-2.3440685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0869861) q[2];
sx q[2];
rz(-0.84422529) q[2];
sx q[2];
rz(-1.1620713) q[2];
rz(-1.6880796) q[3];
sx q[3];
rz(-0.96265692) q[3];
sx q[3];
rz(-1.8614205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42089713) q[0];
sx q[0];
rz(-1.9891885) q[0];
sx q[0];
rz(2.5757117) q[0];
rz(-0.3903009) q[1];
sx q[1];
rz(-1.5719599) q[1];
sx q[1];
rz(-1.8338667) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2824664) q[0];
sx q[0];
rz(-2.0604366) q[0];
sx q[0];
rz(2.7895176) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21251596) q[2];
sx q[2];
rz(-1.9409928) q[2];
sx q[2];
rz(0.0070564673) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2016328) q[1];
sx q[1];
rz(-1.3091413) q[1];
sx q[1];
rz(0.6263588) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0218089) q[3];
sx q[3];
rz(-1.7193188) q[3];
sx q[3];
rz(2.95179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2719443) q[2];
sx q[2];
rz(-0.41676909) q[2];
sx q[2];
rz(1.0375674) q[2];
rz(1.8218254) q[3];
sx q[3];
rz(-2.3128553) q[3];
sx q[3];
rz(-2.493609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8808402) q[0];
sx q[0];
rz(-2.1639731) q[0];
sx q[0];
rz(2.4993437) q[0];
rz(1.3308659) q[1];
sx q[1];
rz(-0.86527491) q[1];
sx q[1];
rz(-0.16255249) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0560682) q[0];
sx q[0];
rz(-0.80636388) q[0];
sx q[0];
rz(1.0170522) q[0];
x q[1];
rz(2.3847975) q[2];
sx q[2];
rz(-0.50163236) q[2];
sx q[2];
rz(1.0197786) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.46025113) q[1];
sx q[1];
rz(-2.6266909) q[1];
sx q[1];
rz(1.9327362) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58329669) q[3];
sx q[3];
rz(-1.5044893) q[3];
sx q[3];
rz(0.60140677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6315397) q[2];
sx q[2];
rz(-1.8659464) q[2];
sx q[2];
rz(-2.1007382) q[2];
rz(0.96796525) q[3];
sx q[3];
rz(-1.0699882) q[3];
sx q[3];
rz(-0.66106838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
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
rz(-2.0816879) q[2];
sx q[2];
rz(-0.89085205) q[2];
sx q[2];
rz(-0.10382626) q[2];
rz(2.6951058) q[3];
sx q[3];
rz(-0.75110186) q[3];
sx q[3];
rz(2.5121381) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
