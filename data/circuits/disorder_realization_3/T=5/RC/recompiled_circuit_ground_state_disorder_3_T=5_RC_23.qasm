OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71896267) q[0];
sx q[0];
rz(-0.29932061) q[0];
sx q[0];
rz(0.49459767) q[0];
rz(-1.9994796) q[1];
sx q[1];
rz(-2.1358868) q[1];
sx q[1];
rz(-1.1297273) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62837961) q[0];
sx q[0];
rz(-0.90115479) q[0];
sx q[0];
rz(1.3684526) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0725934) q[2];
sx q[2];
rz(-1.9747726) q[2];
sx q[2];
rz(-0.4046658) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.59334457) q[1];
sx q[1];
rz(-1.8234089) q[1];
sx q[1];
rz(-2.3139364) q[1];
rz(-1.0315597) q[3];
sx q[3];
rz(-0.72537106) q[3];
sx q[3];
rz(1.6417208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.31892458) q[2];
sx q[2];
rz(-1.8225887) q[2];
sx q[2];
rz(0.85282105) q[2];
rz(-1.8114113) q[3];
sx q[3];
rz(-2.4460402) q[3];
sx q[3];
rz(-3.0947963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.3332719) q[0];
sx q[0];
rz(-3.0905753) q[0];
sx q[0];
rz(1.6774696) q[0];
rz(-1.6630215) q[1];
sx q[1];
rz(-1.9824948) q[1];
sx q[1];
rz(1.2082072) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0883597) q[0];
sx q[0];
rz(-2.3908983) q[0];
sx q[0];
rz(0.73221598) q[0];
rz(-pi) q[1];
rz(1.0958395) q[2];
sx q[2];
rz(-1.7647226) q[2];
sx q[2];
rz(-0.79481193) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7866655) q[1];
sx q[1];
rz(-1.3043205) q[1];
sx q[1];
rz(-3.0158426) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7977116) q[3];
sx q[3];
rz(-2.6626427) q[3];
sx q[3];
rz(-2.6987181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.53604424) q[2];
sx q[2];
rz(-0.41993419) q[2];
sx q[2];
rz(-1.4844683) q[2];
rz(0.015497192) q[3];
sx q[3];
rz(-1.9280547) q[3];
sx q[3];
rz(1.1789471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.991796) q[0];
sx q[0];
rz(-1.3301671) q[0];
sx q[0];
rz(2.8699744) q[0];
rz(-0.89667165) q[1];
sx q[1];
rz(-0.50183693) q[1];
sx q[1];
rz(-1.2976049) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0947726) q[0];
sx q[0];
rz(-1.5808788) q[0];
sx q[0];
rz(-2.0338414) q[0];
rz(1.1900558) q[2];
sx q[2];
rz(-2.0443161) q[2];
sx q[2];
rz(-1.8434332) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1101164) q[1];
sx q[1];
rz(-1.8175442) q[1];
sx q[1];
rz(3.0241248) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5758697) q[3];
sx q[3];
rz(-0.80883615) q[3];
sx q[3];
rz(0.047732959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4012332) q[2];
sx q[2];
rz(-0.77784246) q[2];
sx q[2];
rz(-3.0878301) q[2];
rz(-1.2939804) q[3];
sx q[3];
rz(-2.3432178) q[3];
sx q[3];
rz(2.9062041) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9409222) q[0];
sx q[0];
rz(-2.5735452) q[0];
sx q[0];
rz(-1.4642375) q[0];
rz(-2.0109406) q[1];
sx q[1];
rz(-0.73431763) q[1];
sx q[1];
rz(-0.11437036) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0972041) q[0];
sx q[0];
rz(-2.1425793) q[0];
sx q[0];
rz(-1.0685789) q[0];
rz(0.61538265) q[2];
sx q[2];
rz(-2.3856731) q[2];
sx q[2];
rz(-2.4494954) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.83437551) q[1];
sx q[1];
rz(-2.465472) q[1];
sx q[1];
rz(-1.5285049) q[1];
rz(-pi) q[2];
rz(-0.19802494) q[3];
sx q[3];
rz(-1.5643483) q[3];
sx q[3];
rz(-2.7355268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6681246) q[2];
sx q[2];
rz(-1.5359842) q[2];
sx q[2];
rz(-3.065897) q[2];
rz(2.7068052) q[3];
sx q[3];
rz(-1.1554759) q[3];
sx q[3];
rz(-3.0619612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39744034) q[0];
sx q[0];
rz(-1.3120774) q[0];
sx q[0];
rz(0.42006668) q[0];
rz(-2.9282667) q[1];
sx q[1];
rz(-1.7330287) q[1];
sx q[1];
rz(1.8992281) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38753375) q[0];
sx q[0];
rz(-0.84897572) q[0];
sx q[0];
rz(2.5639064) q[0];
rz(-0.43102805) q[2];
sx q[2];
rz(-1.3044453) q[2];
sx q[2];
rz(-2.9961627) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2957569) q[1];
sx q[1];
rz(-1.6713665) q[1];
sx q[1];
rz(-0.015458903) q[1];
rz(-2.9202634) q[3];
sx q[3];
rz(-1.6175272) q[3];
sx q[3];
rz(-2.5525301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8215948) q[2];
sx q[2];
rz(-0.91511202) q[2];
sx q[2];
rz(-0.37459174) q[2];
rz(2.1290667) q[3];
sx q[3];
rz(-2.7495224) q[3];
sx q[3];
rz(-2.4969126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0300765) q[0];
sx q[0];
rz(-2.6426297) q[0];
sx q[0];
rz(-0.43055713) q[0];
rz(-2.997609) q[1];
sx q[1];
rz(-1.0824243) q[1];
sx q[1];
rz(0.63124257) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78203653) q[0];
sx q[0];
rz(-0.61114531) q[0];
sx q[0];
rz(3.0474328) q[0];
rz(1.1503136) q[2];
sx q[2];
rz(-1.4373206) q[2];
sx q[2];
rz(-1.6949289) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4858125) q[1];
sx q[1];
rz(-2.1596599) q[1];
sx q[1];
rz(-0.9649802) q[1];
x q[2];
rz(-0.33511929) q[3];
sx q[3];
rz(-1.1405754) q[3];
sx q[3];
rz(0.29725257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.97079078) q[2];
sx q[2];
rz(-2.0826714) q[2];
sx q[2];
rz(1.7010752) q[2];
rz(1.3614281) q[3];
sx q[3];
rz(-2.0378588) q[3];
sx q[3];
rz(-0.48719278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68706566) q[0];
sx q[0];
rz(-0.28959689) q[0];
sx q[0];
rz(-0.92426306) q[0];
rz(-0.29280064) q[1];
sx q[1];
rz(-1.9127138) q[1];
sx q[1];
rz(2.3407095) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7445114) q[0];
sx q[0];
rz(-1.9430706) q[0];
sx q[0];
rz(0.27639322) q[0];
rz(-0.34556324) q[2];
sx q[2];
rz(-0.93831944) q[2];
sx q[2];
rz(-1.8008302) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2877601) q[1];
sx q[1];
rz(-1.4656855) q[1];
sx q[1];
rz(-1.5413324) q[1];
rz(-0.022014736) q[3];
sx q[3];
rz(-1.5979294) q[3];
sx q[3];
rz(-2.5282945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.31876365) q[2];
sx q[2];
rz(-1.9209361) q[2];
sx q[2];
rz(-1.2804383) q[2];
rz(0.66796962) q[3];
sx q[3];
rz(-0.48798713) q[3];
sx q[3];
rz(2.7282696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8845344) q[0];
sx q[0];
rz(-1.9030544) q[0];
sx q[0];
rz(-1.1827693) q[0];
rz(-0.23712748) q[1];
sx q[1];
rz(-0.10219899) q[1];
sx q[1];
rz(-3.0822486) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0421365) q[0];
sx q[0];
rz(-1.2980532) q[0];
sx q[0];
rz(-3.0173561) q[0];
rz(-pi) q[1];
rz(1.0419215) q[2];
sx q[2];
rz(-2.2239074) q[2];
sx q[2];
rz(-0.83335857) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8000477) q[1];
sx q[1];
rz(-2.8570685) q[1];
sx q[1];
rz(-1.9519898) q[1];
rz(-pi) q[2];
rz(-1.8944727) q[3];
sx q[3];
rz(-1.6442862) q[3];
sx q[3];
rz(-0.86548129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5817029) q[2];
sx q[2];
rz(-0.54055944) q[2];
sx q[2];
rz(1.3524559) q[2];
rz(-1.7743568) q[3];
sx q[3];
rz(-1.4227899) q[3];
sx q[3];
rz(2.2860315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8597813) q[0];
sx q[0];
rz(-2.3563522) q[0];
sx q[0];
rz(-0.45968858) q[0];
rz(0.028060878) q[1];
sx q[1];
rz(-1.9919845) q[1];
sx q[1];
rz(-1.185816) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59212771) q[0];
sx q[0];
rz(-0.26927265) q[0];
sx q[0];
rz(0.16945355) q[0];
x q[1];
rz(1.6721647) q[2];
sx q[2];
rz(-1.4918054) q[2];
sx q[2];
rz(1.7192507) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3188022) q[1];
sx q[1];
rz(-1.0981961) q[1];
sx q[1];
rz(-0.54293718) q[1];
x q[2];
rz(2.5236058) q[3];
sx q[3];
rz(-1.8027961) q[3];
sx q[3];
rz(0.31832507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.49843732) q[2];
sx q[2];
rz(-0.87673059) q[2];
sx q[2];
rz(0.15677491) q[2];
rz(-0.59569851) q[3];
sx q[3];
rz(-0.34606338) q[3];
sx q[3];
rz(3.116385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(-2.6237685) q[0];
sx q[0];
rz(-1.0305923) q[0];
sx q[0];
rz(2.9577756) q[0];
rz(3.0635762) q[1];
sx q[1];
rz(-2.0842431) q[1];
sx q[1];
rz(-0.26430166) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3027356) q[0];
sx q[0];
rz(-0.6386916) q[0];
sx q[0];
rz(0.46052082) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0980706) q[2];
sx q[2];
rz(-0.95607483) q[2];
sx q[2];
rz(1.5103024) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7122927) q[1];
sx q[1];
rz(-1.4136864) q[1];
sx q[1];
rz(0.18204851) q[1];
rz(-pi) q[2];
rz(2.6464858) q[3];
sx q[3];
rz(-1.6551842) q[3];
sx q[3];
rz(1.6364678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2901624) q[2];
sx q[2];
rz(-2.1992079) q[2];
sx q[2];
rz(0.22872049) q[2];
rz(1.4043572) q[3];
sx q[3];
rz(-0.93969932) q[3];
sx q[3];
rz(-3.0586045) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9901154) q[0];
sx q[0];
rz(-0.84828068) q[0];
sx q[0];
rz(-1.620851) q[0];
rz(-2.4304541) q[1];
sx q[1];
rz(-1.8322721) q[1];
sx q[1];
rz(-2.4087404) q[1];
rz(3.0971211) q[2];
sx q[2];
rz(-1.0723249) q[2];
sx q[2];
rz(1.9981801) q[2];
rz(-2.5855999) q[3];
sx q[3];
rz(-2.0423732) q[3];
sx q[3];
rz(-1.1788551) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
