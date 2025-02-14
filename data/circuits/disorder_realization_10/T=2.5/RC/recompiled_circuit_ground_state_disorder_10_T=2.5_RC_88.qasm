OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.77914971) q[0];
sx q[0];
rz(2.7326512) q[0];
sx q[0];
rz(10.41806) q[0];
rz(2.9945057) q[1];
sx q[1];
rz(-2.1451201) q[1];
sx q[1];
rz(1.7239404) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.111747) q[0];
sx q[0];
rz(-0.8518712) q[0];
sx q[0];
rz(1.0875888) q[0];
rz(1.0101914) q[2];
sx q[2];
rz(-2.1935138) q[2];
sx q[2];
rz(1.6499008) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5681709) q[1];
sx q[1];
rz(-2.4470098) q[1];
sx q[1];
rz(0.43412368) q[1];
rz(-pi) q[2];
x q[2];
rz(0.032790498) q[3];
sx q[3];
rz(-1.7049978) q[3];
sx q[3];
rz(0.57574948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.950497) q[2];
sx q[2];
rz(-1.3206626) q[2];
sx q[2];
rz(-1.0547868) q[2];
rz(-2.6705006) q[3];
sx q[3];
rz(-1.6867009) q[3];
sx q[3];
rz(2.2733222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7489814) q[0];
sx q[0];
rz(-1.457021) q[0];
sx q[0];
rz(-2.9066322) q[0];
rz(-1.8670392) q[1];
sx q[1];
rz(-1.8320558) q[1];
sx q[1];
rz(0.68769208) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10768453) q[0];
sx q[0];
rz(-1.8982197) q[0];
sx q[0];
rz(-1.641724) q[0];
rz(-2.2830487) q[2];
sx q[2];
rz(-1.1919824) q[2];
sx q[2];
rz(2.130098) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.54476538) q[1];
sx q[1];
rz(-3.1278962) q[1];
sx q[1];
rz(1.9910452) q[1];
rz(3.0491203) q[3];
sx q[3];
rz(-0.95917976) q[3];
sx q[3];
rz(-0.05449748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71622744) q[2];
sx q[2];
rz(-1.8391515) q[2];
sx q[2];
rz(-0.20935527) q[2];
rz(2.1685205) q[3];
sx q[3];
rz(-0.89997411) q[3];
sx q[3];
rz(0.59276855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3000325) q[0];
sx q[0];
rz(-0.40599269) q[0];
sx q[0];
rz(-5/(11*pi)) q[0];
rz(-1.2380098) q[1];
sx q[1];
rz(-0.7716476) q[1];
sx q[1];
rz(3.0498116) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9301355) q[0];
sx q[0];
rz(-1.9152904) q[0];
sx q[0];
rz(0.99665595) q[0];
x q[1];
rz(-0.94987671) q[2];
sx q[2];
rz(-0.11813049) q[2];
sx q[2];
rz(2.0565513) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5621693) q[1];
sx q[1];
rz(-1.7618638) q[1];
sx q[1];
rz(-1.3538989) q[1];
rz(1.7157406) q[3];
sx q[3];
rz(-0.43625375) q[3];
sx q[3];
rz(3.1097092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7564275) q[2];
sx q[2];
rz(-1.6331208) q[2];
sx q[2];
rz(-2.7623994) q[2];
rz(-2.3181629) q[3];
sx q[3];
rz(-2.7354666) q[3];
sx q[3];
rz(2.8520975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32651153) q[0];
sx q[0];
rz(-0.87711763) q[0];
sx q[0];
rz(2.1429578) q[0];
rz(0.48209349) q[1];
sx q[1];
rz(-2.1908052) q[1];
sx q[1];
rz(-1.8236209) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4559193) q[0];
sx q[0];
rz(-0.64252526) q[0];
sx q[0];
rz(2.5701017) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7779269) q[2];
sx q[2];
rz(-2.5245856) q[2];
sx q[2];
rz(-1.3777767) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.91693976) q[1];
sx q[1];
rz(-0.5782776) q[1];
sx q[1];
rz(-2.7507332) q[1];
rz(-pi) q[2];
rz(-0.51709922) q[3];
sx q[3];
rz(-0.86315599) q[3];
sx q[3];
rz(1.4003889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2202997) q[2];
sx q[2];
rz(-0.90196323) q[2];
sx q[2];
rz(-2.656142) q[2];
rz(2.7576533) q[3];
sx q[3];
rz(-1.5555614) q[3];
sx q[3];
rz(0.78401047) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0601198) q[0];
sx q[0];
rz(-0.10426846) q[0];
sx q[0];
rz(3.1299348) q[0];
rz(-3.0043789) q[1];
sx q[1];
rz(-1.5882746) q[1];
sx q[1];
rz(-1.8395909) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0948167) q[0];
sx q[0];
rz(-1.6794378) q[0];
sx q[0];
rz(-0.18414761) q[0];
rz(-pi) q[1];
rz(-0.9658034) q[2];
sx q[2];
rz(-2.0070672) q[2];
sx q[2];
rz(0.80096132) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.80655386) q[1];
sx q[1];
rz(-1.2590027) q[1];
sx q[1];
rz(3.0278518) q[1];
rz(-0.75826606) q[3];
sx q[3];
rz(-0.50095424) q[3];
sx q[3];
rz(-1.6538594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4181218) q[2];
sx q[2];
rz(-1.301441) q[2];
sx q[2];
rz(2.8483025) q[2];
rz(-0.063118525) q[3];
sx q[3];
rz(-1.2226353) q[3];
sx q[3];
rz(-0.89490923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1191331) q[0];
sx q[0];
rz(-2.8526511) q[0];
sx q[0];
rz(-0.038507842) q[0];
rz(-2.0571845) q[1];
sx q[1];
rz(-0.42250982) q[1];
sx q[1];
rz(2.1748621) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9597067) q[0];
sx q[0];
rz(-0.27508914) q[0];
sx q[0];
rz(-0.77232124) q[0];
rz(-pi) q[1];
rz(1.837908) q[2];
sx q[2];
rz(-1.583364) q[2];
sx q[2];
rz(-1.723701) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.048592) q[1];
sx q[1];
rz(-1.9680259) q[1];
sx q[1];
rz(-2.7723958) q[1];
rz(1.8943664) q[3];
sx q[3];
rz(-0.67646356) q[3];
sx q[3];
rz(0.25581365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66095573) q[2];
sx q[2];
rz(-0.2350685) q[2];
sx q[2];
rz(-0.98141518) q[2];
rz(-1.6437982) q[3];
sx q[3];
rz(-1.8794329) q[3];
sx q[3];
rz(1.5952544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82866955) q[0];
sx q[0];
rz(-0.68083119) q[0];
sx q[0];
rz(-0.31461) q[0];
rz(-2.9529052) q[1];
sx q[1];
rz(-0.43470946) q[1];
sx q[1];
rz(-2.6729118) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73891) q[0];
sx q[0];
rz(-0.27640585) q[0];
sx q[0];
rz(-0.49582793) q[0];
rz(-pi) q[1];
rz(1.5596703) q[2];
sx q[2];
rz(-2.9738148) q[2];
sx q[2];
rz(1.1388701) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3568253) q[1];
sx q[1];
rz(-1.873849) q[1];
sx q[1];
rz(-2.3222938) q[1];
rz(1.8953034) q[3];
sx q[3];
rz(-1.1238928) q[3];
sx q[3];
rz(0.75845065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.20588747) q[2];
sx q[2];
rz(-2.2998655) q[2];
sx q[2];
rz(0.97071281) q[2];
rz(-0.5361706) q[3];
sx q[3];
rz(-1.2090809) q[3];
sx q[3];
rz(-2.4255588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4416874) q[0];
sx q[0];
rz(-2.5869885) q[0];
sx q[0];
rz(1.6957138) q[0];
rz(-2.1031759) q[1];
sx q[1];
rz(-1.1035792) q[1];
sx q[1];
rz(0.081092484) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86085194) q[0];
sx q[0];
rz(-0.74374712) q[0];
sx q[0];
rz(2.70983) q[0];
x q[1];
rz(-1.8683226) q[2];
sx q[2];
rz(-2.2045481) q[2];
sx q[2];
rz(0.81281205) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9633741) q[1];
sx q[1];
rz(-1.4556754) q[1];
sx q[1];
rz(-1.2769775) q[1];
x q[2];
rz(-0.57402667) q[3];
sx q[3];
rz(-1.9732631) q[3];
sx q[3];
rz(-2.5695679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9424092) q[2];
sx q[2];
rz(-2.3641391) q[2];
sx q[2];
rz(-2.8988163) q[2];
rz(-1.5593922) q[3];
sx q[3];
rz(-2.715761) q[3];
sx q[3];
rz(1.8886214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2389857) q[0];
sx q[0];
rz(-1.0317529) q[0];
sx q[0];
rz(1.2549988) q[0];
rz(-1.3764508) q[1];
sx q[1];
rz(-1.0245198) q[1];
sx q[1];
rz(2.4741516) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8534626) q[0];
sx q[0];
rz(-1.9024182) q[0];
sx q[0];
rz(2.1037071) q[0];
rz(-2.6450577) q[2];
sx q[2];
rz(-1.6081972) q[2];
sx q[2];
rz(-0.9031682) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.046958663) q[1];
sx q[1];
rz(-2.486645) q[1];
sx q[1];
rz(2.0597233) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7605704) q[3];
sx q[3];
rz(-1.8748611) q[3];
sx q[3];
rz(-1.1883433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5218375) q[2];
sx q[2];
rz(-0.73885584) q[2];
sx q[2];
rz(2.6832306) q[2];
rz(-1.5357664) q[3];
sx q[3];
rz(-2.063844) q[3];
sx q[3];
rz(0.88466907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2430724) q[0];
sx q[0];
rz(-2.5694818) q[0];
sx q[0];
rz(2.7761053) q[0];
rz(-2.1416523) q[1];
sx q[1];
rz(-1.702407) q[1];
sx q[1];
rz(-1.684729) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6479051) q[0];
sx q[0];
rz(-2.5468605) q[0];
sx q[0];
rz(-1.2059284) q[0];
x q[1];
rz(1.8387055) q[2];
sx q[2];
rz(-1.678907) q[2];
sx q[2];
rz(-1.8927758) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.748865) q[1];
sx q[1];
rz(-2.3727149) q[1];
sx q[1];
rz(-0.25090541) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69646466) q[3];
sx q[3];
rz(-1.2035771) q[3];
sx q[3];
rz(0.61456481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0008056) q[2];
sx q[2];
rz(-1.8203338) q[2];
sx q[2];
rz(-1.0055536) q[2];
rz(2.4908861) q[3];
sx q[3];
rz(-1.4395827) q[3];
sx q[3];
rz(2.2910291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-3.0477796) q[0];
sx q[0];
rz(-0.78421264) q[0];
sx q[0];
rz(-1.0017851) q[0];
rz(1.7306937) q[1];
sx q[1];
rz(-1.7041364) q[1];
sx q[1];
rz(-2.0028353) q[1];
rz(-1.9000353) q[2];
sx q[2];
rz(-1.6816947) q[2];
sx q[2];
rz(-1.1244863) q[2];
rz(1.6639573) q[3];
sx q[3];
rz(-1.4668844) q[3];
sx q[3];
rz(1.8528406) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
