OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39855555) q[0];
sx q[0];
rz(-0.86657137) q[0];
sx q[0];
rz(0.33696365) q[0];
rz(0.26951867) q[1];
sx q[1];
rz(-2.1989006) q[1];
sx q[1];
rz(1.2595133) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058858697) q[0];
sx q[0];
rz(-2.132651) q[0];
sx q[0];
rz(-0.64942067) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3607698) q[2];
sx q[2];
rz(-2.2968596) q[2];
sx q[2];
rz(-0.10440102) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4866529) q[1];
sx q[1];
rz(-1.6702685) q[1];
sx q[1];
rz(2.59511) q[1];
x q[2];
rz(-0.41347031) q[3];
sx q[3];
rz(-1.0116825) q[3];
sx q[3];
rz(-2.632189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.17244615) q[2];
sx q[2];
rz(-1.5643876) q[2];
sx q[2];
rz(1.2338314) q[2];
rz(-1.6110169) q[3];
sx q[3];
rz(-1.4065892) q[3];
sx q[3];
rz(-0.56510258) q[3];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81561404) q[0];
sx q[0];
rz(-0.9650721) q[0];
sx q[0];
rz(2.8976029) q[0];
rz(-1.9832393) q[1];
sx q[1];
rz(-2.127425) q[1];
sx q[1];
rz(0.44833952) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8566003) q[0];
sx q[0];
rz(-2.0584848) q[0];
sx q[0];
rz(1.3707121) q[0];
x q[1];
rz(-2.3116287) q[2];
sx q[2];
rz(-1.1051851) q[2];
sx q[2];
rz(2.5453886) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2177451) q[1];
sx q[1];
rz(-2.4051146) q[1];
sx q[1];
rz(2.272241) q[1];
rz(2.6579992) q[3];
sx q[3];
rz(-1.075287) q[3];
sx q[3];
rz(-2.7413975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1925194) q[2];
sx q[2];
rz(-3.101888) q[2];
sx q[2];
rz(0.81083361) q[2];
rz(-3.132498) q[3];
sx q[3];
rz(-1.5261212) q[3];
sx q[3];
rz(2.8240589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72967616) q[0];
sx q[0];
rz(-1.7971973) q[0];
sx q[0];
rz(0.94773951) q[0];
rz(-0.89836994) q[1];
sx q[1];
rz(-0.71304524) q[1];
sx q[1];
rz(-0.20634849) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2284828) q[0];
sx q[0];
rz(-0.75487274) q[0];
sx q[0];
rz(1.5819527) q[0];
x q[1];
rz(2.3487665) q[2];
sx q[2];
rz(-1.8603051) q[2];
sx q[2];
rz(-3.135856) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5134597) q[1];
sx q[1];
rz(-2.0230789) q[1];
sx q[1];
rz(2.9123405) q[1];
x q[2];
rz(-3.0037155) q[3];
sx q[3];
rz(-2.8123224) q[3];
sx q[3];
rz(-0.38789685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.72948939) q[2];
sx q[2];
rz(-0.99200839) q[2];
sx q[2];
rz(-0.98423973) q[2];
rz(0.49736831) q[3];
sx q[3];
rz(-1.8822742) q[3];
sx q[3];
rz(-0.74500144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3618149) q[0];
sx q[0];
rz(-3.0497157) q[0];
sx q[0];
rz(2.9352557) q[0];
rz(1.6774718) q[1];
sx q[1];
rz(-0.76281491) q[1];
sx q[1];
rz(0.62613553) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3040309) q[0];
sx q[0];
rz(-1.921424) q[0];
sx q[0];
rz(2.3239345) q[0];
x q[1];
rz(-0.10278662) q[2];
sx q[2];
rz(-1.572552) q[2];
sx q[2];
rz(0.87550113) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5474769) q[1];
sx q[1];
rz(-1.0924589) q[1];
sx q[1];
rz(-2.5376471) q[1];
rz(-1.9575084) q[3];
sx q[3];
rz(-1.6624647) q[3];
sx q[3];
rz(2.693585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3760486) q[2];
sx q[2];
rz(-0.77771336) q[2];
sx q[2];
rz(-2.1518339) q[2];
rz(-2.1073714) q[3];
sx q[3];
rz(-1.1690305) q[3];
sx q[3];
rz(-0.52556747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.702221) q[0];
sx q[0];
rz(-1.7218497) q[0];
sx q[0];
rz(1.9787582) q[0];
rz(-2.974466) q[1];
sx q[1];
rz(-1.6163328) q[1];
sx q[1];
rz(2.4662245) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1646368) q[0];
sx q[0];
rz(-2.2222493) q[0];
sx q[0];
rz(-1.2950334) q[0];
rz(-pi) q[1];
rz(1.950267) q[2];
sx q[2];
rz(-1.8877875) q[2];
sx q[2];
rz(1.6940885) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31084222) q[1];
sx q[1];
rz(-1.9988339) q[1];
sx q[1];
rz(0.38106783) q[1];
rz(2.8699304) q[3];
sx q[3];
rz(-2.144648) q[3];
sx q[3];
rz(2.7558307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2165788) q[2];
sx q[2];
rz(-1.719097) q[2];
sx q[2];
rz(0.50755802) q[2];
rz(-3.1223068) q[3];
sx q[3];
rz(-2.3465893) q[3];
sx q[3];
rz(-1.9222586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5386706) q[0];
sx q[0];
rz(-0.62748533) q[0];
sx q[0];
rz(-0.70166171) q[0];
rz(-0.64274669) q[1];
sx q[1];
rz(-1.1593436) q[1];
sx q[1];
rz(1.3210375) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3747201) q[0];
sx q[0];
rz(-1.5581616) q[0];
sx q[0];
rz(0.94324525) q[0];
x q[1];
rz(2.7847544) q[2];
sx q[2];
rz(-2.2388487) q[2];
sx q[2];
rz(-0.84260637) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5112662) q[1];
sx q[1];
rz(-2.0980586) q[1];
sx q[1];
rz(-2.7106337) q[1];
rz(1.2679638) q[3];
sx q[3];
rz(-2.3583671) q[3];
sx q[3];
rz(0.86963213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.295149) q[2];
sx q[2];
rz(-0.48131338) q[2];
sx q[2];
rz(-2.5018137) q[2];
rz(-2.4844737) q[3];
sx q[3];
rz(-1.6491456) q[3];
sx q[3];
rz(0.28321987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039728634) q[0];
sx q[0];
rz(-1.855408) q[0];
sx q[0];
rz(-0.78045994) q[0];
rz(1.9038433) q[1];
sx q[1];
rz(-1.8433808) q[1];
sx q[1];
rz(-2.9775528) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4481159) q[0];
sx q[0];
rz(-2.0882029) q[0];
sx q[0];
rz(2.8679001) q[0];
rz(0.41167792) q[2];
sx q[2];
rz(-1.1601761) q[2];
sx q[2];
rz(-0.1146929) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2795231) q[1];
sx q[1];
rz(-1.0633103) q[1];
sx q[1];
rz(-2.1699778) q[1];
rz(-pi) q[2];
rz(0.45312552) q[3];
sx q[3];
rz(-2.2966566) q[3];
sx q[3];
rz(2.912584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26874545) q[2];
sx q[2];
rz(-1.1970604) q[2];
sx q[2];
rz(0.8806814) q[2];
rz(1.667048) q[3];
sx q[3];
rz(-2.0737952) q[3];
sx q[3];
rz(1.4145981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47377652) q[0];
sx q[0];
rz(-0.18874636) q[0];
sx q[0];
rz(2.012398) q[0];
rz(-2.1881564) q[1];
sx q[1];
rz(-2.2013181) q[1];
sx q[1];
rz(-2.5530691) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.584483) q[0];
sx q[0];
rz(-0.63562387) q[0];
sx q[0];
rz(2.0899713) q[0];
rz(-2.1641974) q[2];
sx q[2];
rz(-0.39476141) q[2];
sx q[2];
rz(0.66937689) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2413832) q[1];
sx q[1];
rz(-1.8892678) q[1];
sx q[1];
rz(-2.4584944) q[1];
x q[2];
rz(-1.3465914) q[3];
sx q[3];
rz(-1.5191169) q[3];
sx q[3];
rz(1.3286852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9863161) q[2];
sx q[2];
rz(-1.0078112) q[2];
sx q[2];
rz(1.0922095) q[2];
rz(0.80882597) q[3];
sx q[3];
rz(-2.0598965) q[3];
sx q[3];
rz(1.4441747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5876193) q[0];
sx q[0];
rz(-2.0646136) q[0];
sx q[0];
rz(-2.4699566) q[0];
rz(2.6486168) q[1];
sx q[1];
rz(-0.88075811) q[1];
sx q[1];
rz(1.902045) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.400378) q[0];
sx q[0];
rz(-1.5733871) q[0];
sx q[0];
rz(1.3335614) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0918248) q[2];
sx q[2];
rz(-2.9441212) q[2];
sx q[2];
rz(1.2021499) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.904027) q[1];
sx q[1];
rz(-2.3859302) q[1];
sx q[1];
rz(0.66046884) q[1];
rz(-1.4080234) q[3];
sx q[3];
rz(-2.3331169) q[3];
sx q[3];
rz(2.2089434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.12941831) q[2];
sx q[2];
rz(-1.6678026) q[2];
sx q[2];
rz(-2.5950281) q[2];
rz(-2.2187345) q[3];
sx q[3];
rz(-2.512629) q[3];
sx q[3];
rz(1.9950689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4839812) q[0];
sx q[0];
rz(-2.681499) q[0];
sx q[0];
rz(-0.42903236) q[0];
rz(-1.2826762) q[1];
sx q[1];
rz(-1.7762643) q[1];
sx q[1];
rz(3.0603337) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0711533) q[0];
sx q[0];
rz(-2.9550046) q[0];
sx q[0];
rz(2.7036315) q[0];
rz(-pi) q[1];
rz(-0.0325412) q[2];
sx q[2];
rz(-2.0306808) q[2];
sx q[2];
rz(-0.77600098) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.024549896) q[1];
sx q[1];
rz(-2.3094842) q[1];
sx q[1];
rz(1.0067654) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0662543) q[3];
sx q[3];
rz(-1.9849376) q[3];
sx q[3];
rz(0.6739236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73605865) q[2];
sx q[2];
rz(-2.2106876) q[2];
sx q[2];
rz(2.0299358) q[2];
rz(-1.1692272) q[3];
sx q[3];
rz(-2.4284095) q[3];
sx q[3];
rz(2.9986103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1375167) q[0];
sx q[0];
rz(-0.62340323) q[0];
sx q[0];
rz(2.2829983) q[0];
rz(0.89523347) q[1];
sx q[1];
rz(-1.3139071) q[1];
sx q[1];
rz(0.039176686) q[1];
rz(-0.011717144) q[2];
sx q[2];
rz(-1.313971) q[2];
sx q[2];
rz(1.32709) q[2];
rz(0.17651996) q[3];
sx q[3];
rz(-0.64809496) q[3];
sx q[3];
rz(2.5213836) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
