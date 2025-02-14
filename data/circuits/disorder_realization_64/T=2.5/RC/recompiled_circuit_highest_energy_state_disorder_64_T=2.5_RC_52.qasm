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
rz(-0.023802726) q[0];
sx q[0];
rz(-2.0355712) q[0];
sx q[0];
rz(0.7769146) q[0];
rz(1.39224) q[1];
sx q[1];
rz(-1.3148146) q[1];
sx q[1];
rz(-0.97631747) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38589729) q[0];
sx q[0];
rz(-1.160826) q[0];
sx q[0];
rz(0.14571054) q[0];
x q[1];
rz(-0.10586057) q[2];
sx q[2];
rz(-2.5753394) q[2];
sx q[2];
rz(0.98041269) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1433409) q[1];
sx q[1];
rz(-0.51330459) q[1];
sx q[1];
rz(-2.1480888) q[1];
x q[2];
rz(-0.032194897) q[3];
sx q[3];
rz(-2.1436084) q[3];
sx q[3];
rz(2.9762852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6355847) q[2];
sx q[2];
rz(-3.0878461) q[2];
sx q[2];
rz(-2.7347943) q[2];
rz(0.16945101) q[3];
sx q[3];
rz(-2.6120766) q[3];
sx q[3];
rz(-1.0725526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88103831) q[0];
sx q[0];
rz(-0.23242234) q[0];
sx q[0];
rz(3.1378003) q[0];
rz(-3.0637528) q[1];
sx q[1];
rz(-2.4816315) q[1];
sx q[1];
rz(-2.8357764) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77340375) q[0];
sx q[0];
rz(-2.0819252) q[0];
sx q[0];
rz(-2.9298733) q[0];
rz(2.3897116) q[2];
sx q[2];
rz(-1.4724178) q[2];
sx q[2];
rz(1.4954612) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7485436) q[1];
sx q[1];
rz(-2.3094588) q[1];
sx q[1];
rz(-2.119675) q[1];
rz(-0.81021328) q[3];
sx q[3];
rz(-0.088220291) q[3];
sx q[3];
rz(-2.719413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1514312) q[2];
sx q[2];
rz(-1.3913245) q[2];
sx q[2];
rz(-2.4988417) q[2];
rz(-3.0691872) q[3];
sx q[3];
rz(-1.0591155) q[3];
sx q[3];
rz(1.36093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-3.0533503) q[0];
sx q[0];
rz(-0.14277661) q[0];
sx q[0];
rz(2.9852168) q[0];
rz(-3.1001672) q[1];
sx q[1];
rz(-0.62774575) q[1];
sx q[1];
rz(1.5511537) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0254733) q[0];
sx q[0];
rz(-1.0295233) q[0];
sx q[0];
rz(1.812029) q[0];
rz(-pi) q[1];
rz(2.4855108) q[2];
sx q[2];
rz(-1.6245884) q[2];
sx q[2];
rz(-2.057586) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9030021) q[1];
sx q[1];
rz(-2.6005473) q[1];
sx q[1];
rz(-0.21958406) q[1];
x q[2];
rz(0.66371347) q[3];
sx q[3];
rz(-2.1557689) q[3];
sx q[3];
rz(3.0970517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40858832) q[2];
sx q[2];
rz(-2.6027347) q[2];
sx q[2];
rz(-0.020922529) q[2];
rz(-2.9544592) q[3];
sx q[3];
rz(-2.9360866) q[3];
sx q[3];
rz(0.11370295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9631831) q[0];
sx q[0];
rz(-2.8091176) q[0];
sx q[0];
rz(0.43854976) q[0];
rz(1.6167538) q[1];
sx q[1];
rz(-2.8068145) q[1];
sx q[1];
rz(-0.24510342) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.995718) q[0];
sx q[0];
rz(-2.2882713) q[0];
sx q[0];
rz(-0.11086734) q[0];
rz(-pi) q[1];
rz(1.4133873) q[2];
sx q[2];
rz(-1.206996) q[2];
sx q[2];
rz(0.25870332) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4192701) q[1];
sx q[1];
rz(-2.5825204) q[1];
sx q[1];
rz(1.8354227) q[1];
x q[2];
rz(1.2538337) q[3];
sx q[3];
rz(-2.2989103) q[3];
sx q[3];
rz(1.1945981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.601292) q[2];
sx q[2];
rz(-0.43593323) q[2];
sx q[2];
rz(-2.8098246) q[2];
rz(0.48745421) q[3];
sx q[3];
rz(-1.0468227) q[3];
sx q[3];
rz(0.99307466) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29193923) q[0];
sx q[0];
rz(-1.4537469) q[0];
sx q[0];
rz(-0.77350235) q[0];
rz(1.9603112) q[1];
sx q[1];
rz(-0.14207323) q[1];
sx q[1];
rz(-1.7519417) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31598241) q[0];
sx q[0];
rz(-1.1706691) q[0];
sx q[0];
rz(-1.9670301) q[0];
rz(-0.90221407) q[2];
sx q[2];
rz(-0.83219516) q[2];
sx q[2];
rz(0.26022831) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.70374278) q[1];
sx q[1];
rz(-1.4885167) q[1];
sx q[1];
rz(2.6668496) q[1];
rz(-1.3895274) q[3];
sx q[3];
rz(-1.3366404) q[3];
sx q[3];
rz(-0.93269809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2191849) q[2];
sx q[2];
rz(-1.9318523) q[2];
sx q[2];
rz(2.6694471) q[2];
rz(-1.8426497) q[3];
sx q[3];
rz(-1.2932212) q[3];
sx q[3];
rz(2.3310272) q[3];
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
rz(pi/2) q[0];
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
rz(2.2209114) q[0];
sx q[0];
rz(-2.8784316) q[0];
sx q[0];
rz(0.26350185) q[0];
rz(2.0384516) q[1];
sx q[1];
rz(-1.8270854) q[1];
sx q[1];
rz(-2.7679494) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5140681) q[0];
sx q[0];
rz(-1.4980982) q[0];
sx q[0];
rz(0.060427314) q[0];
rz(-pi) q[1];
rz(-1.7079855) q[2];
sx q[2];
rz(-0.73167668) q[2];
sx q[2];
rz(-0.2889932) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9615104) q[1];
sx q[1];
rz(-2.1282548) q[1];
sx q[1];
rz(-2.7649759) q[1];
x q[2];
rz(-0.5881891) q[3];
sx q[3];
rz(-1.08687) q[3];
sx q[3];
rz(1.5671135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.11334795) q[2];
sx q[2];
rz(-2.9714606) q[2];
sx q[2];
rz(-2.587758) q[2];
rz(-1.3977741) q[3];
sx q[3];
rz(-2.5388986) q[3];
sx q[3];
rz(-0.3127313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58436191) q[0];
sx q[0];
rz(-1.0065684) q[0];
sx q[0];
rz(1.2297909) q[0];
rz(2.9025485) q[1];
sx q[1];
rz(-1.6300647) q[1];
sx q[1];
rz(0.30034932) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2984788) q[0];
sx q[0];
rz(-1.686578) q[0];
sx q[0];
rz(2.9779469) q[0];
rz(-2.4392468) q[2];
sx q[2];
rz(-1.9611729) q[2];
sx q[2];
rz(0.74884383) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9508267) q[1];
sx q[1];
rz(-0.015282282) q[1];
sx q[1];
rz(-1.5453668) q[1];
x q[2];
rz(-0.024552931) q[3];
sx q[3];
rz(-2.2835287) q[3];
sx q[3];
rz(2.7803382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.131669) q[2];
sx q[2];
rz(-1.6841623) q[2];
sx q[2];
rz(0.24492502) q[2];
rz(2.6217672) q[3];
sx q[3];
rz(-0.86331415) q[3];
sx q[3];
rz(2.4533217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35172611) q[0];
sx q[0];
rz(-0.34857294) q[0];
sx q[0];
rz(1.1630195) q[0];
rz(3.0746958) q[1];
sx q[1];
rz(-1.6480548) q[1];
sx q[1];
rz(1.012872) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6352309) q[0];
sx q[0];
rz(-1.073045) q[0];
sx q[0];
rz(0.91301163) q[0];
x q[1];
rz(2.9628721) q[2];
sx q[2];
rz(-1.0105437) q[2];
sx q[2];
rz(2.5541039) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0433181) q[1];
sx q[1];
rz(-1.2407082) q[1];
sx q[1];
rz(1.2374452) q[1];
rz(2.4550405) q[3];
sx q[3];
rz(-1.2408537) q[3];
sx q[3];
rz(0.8984962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.11671994) q[2];
sx q[2];
rz(-0.97815424) q[2];
sx q[2];
rz(-0.32279521) q[2];
rz(2.5358477) q[3];
sx q[3];
rz(-2.3492458) q[3];
sx q[3];
rz(-2.7927223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054319687) q[0];
sx q[0];
rz(-0.1467341) q[0];
sx q[0];
rz(-3.1291381) q[0];
rz(2.3948578) q[1];
sx q[1];
rz(-2.2181999) q[1];
sx q[1];
rz(-2.8616203) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0664205) q[0];
sx q[0];
rz(-3.018258) q[0];
sx q[0];
rz(1.9955817) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93918903) q[2];
sx q[2];
rz(-0.52131182) q[2];
sx q[2];
rz(2.0445532) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2586883) q[1];
sx q[1];
rz(-1.8778442) q[1];
sx q[1];
rz(2.3767002) q[1];
rz(-1.9334698) q[3];
sx q[3];
rz(-1.8020523) q[3];
sx q[3];
rz(1.4303007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3296457) q[2];
sx q[2];
rz(-0.75355607) q[2];
sx q[2];
rz(-2.9296056) q[2];
rz(-2.3181465) q[3];
sx q[3];
rz(-1.6971089) q[3];
sx q[3];
rz(0.25920355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14281808) q[0];
sx q[0];
rz(-0.057567216) q[0];
sx q[0];
rz(0.69277358) q[0];
rz(-0.57299262) q[1];
sx q[1];
rz(-1.7968105) q[1];
sx q[1];
rz(2.7105892) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7438223) q[0];
sx q[0];
rz(-1.6293238) q[0];
sx q[0];
rz(0.55214793) q[0];
rz(-2.0153322) q[2];
sx q[2];
rz(-2.1241786) q[2];
sx q[2];
rz(1.044342) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9622456) q[1];
sx q[1];
rz(-1.3192466) q[1];
sx q[1];
rz(-2.4169282) q[1];
rz(2.1206843) q[3];
sx q[3];
rz(-1.8326933) q[3];
sx q[3];
rz(-0.49347116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7196322) q[2];
sx q[2];
rz(-2.8675291) q[2];
sx q[2];
rz(-2.5813622) q[2];
rz(0.46960056) q[3];
sx q[3];
rz(-2.738939) q[3];
sx q[3];
rz(0.71389055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2847168) q[0];
sx q[0];
rz(-1.7354043) q[0];
sx q[0];
rz(2.0866557) q[0];
rz(-2.3003385) q[1];
sx q[1];
rz(-1.1019191) q[1];
sx q[1];
rz(-0.046774653) q[1];
rz(2.8259059) q[2];
sx q[2];
rz(-1.9236947) q[2];
sx q[2];
rz(-2.535939) q[2];
rz(-0.026997707) q[3];
sx q[3];
rz(-0.91201966) q[3];
sx q[3];
rz(-2.1569679) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
