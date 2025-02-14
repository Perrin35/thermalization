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
rz(-2.4969555) q[0];
sx q[0];
rz(-0.58965373) q[0];
sx q[0];
rz(2.0539505) q[0];
rz(-0.015406869) q[1];
sx q[1];
rz(3.5313731) q[1];
sx q[1];
rz(11.402147) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.312266) q[0];
sx q[0];
rz(-0.13008936) q[0];
sx q[0];
rz(1.6072558) q[0];
rz(0.6911152) q[2];
sx q[2];
rz(-1.8370312) q[2];
sx q[2];
rz(2.63111) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.30689209) q[1];
sx q[1];
rz(-1.0535208) q[1];
sx q[1];
rz(-0.24637988) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.075957493) q[3];
sx q[3];
rz(-1.4531368) q[3];
sx q[3];
rz(2.3880588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2414134) q[2];
sx q[2];
rz(-0.56168491) q[2];
sx q[2];
rz(-2.4955595) q[2];
rz(2.8863886) q[3];
sx q[3];
rz(-0.45707688) q[3];
sx q[3];
rz(1.7469762) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0274886) q[0];
sx q[0];
rz(-0.33559594) q[0];
sx q[0];
rz(0.80701971) q[0];
rz(-1.6837766) q[1];
sx q[1];
rz(-2.5648263) q[1];
sx q[1];
rz(1.3438276) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7248568) q[0];
sx q[0];
rz(-1.2829507) q[0];
sx q[0];
rz(-1.4132981) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0033002) q[2];
sx q[2];
rz(-1.2186745) q[2];
sx q[2];
rz(-2.9825236) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8763153) q[1];
sx q[1];
rz(-2.1723299) q[1];
sx q[1];
rz(-2.5391891) q[1];
rz(0.84193167) q[3];
sx q[3];
rz(-1.6042097) q[3];
sx q[3];
rz(2.9848826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.9606216) q[2];
sx q[2];
rz(-2.0714859) q[2];
sx q[2];
rz(-2.9478659) q[2];
rz(1.5242029) q[3];
sx q[3];
rz(-0.75810713) q[3];
sx q[3];
rz(2.962842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5209565) q[0];
sx q[0];
rz(-2.9036324) q[0];
sx q[0];
rz(-0.97610193) q[0];
rz(-0.8887662) q[1];
sx q[1];
rz(-2.3975394) q[1];
sx q[1];
rz(-0.1730473) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9219396) q[0];
sx q[0];
rz(-1.573607) q[0];
sx q[0];
rz(1.6083184) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6209774) q[2];
sx q[2];
rz(-1.4618502) q[2];
sx q[2];
rz(1.1586939) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5593619) q[1];
sx q[1];
rz(-1.7364864) q[1];
sx q[1];
rz(0.5104602) q[1];
rz(-pi) q[2];
rz(-1.4621696) q[3];
sx q[3];
rz(-1.8914239) q[3];
sx q[3];
rz(2.4214793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.139107) q[2];
sx q[2];
rz(-1.8128914) q[2];
sx q[2];
rz(2.4277182) q[2];
rz(2.930323) q[3];
sx q[3];
rz(-1.0372838) q[3];
sx q[3];
rz(-2.5762288) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36770377) q[0];
sx q[0];
rz(-1.9581032) q[0];
sx q[0];
rz(1.1327889) q[0];
rz(2.6702113) q[1];
sx q[1];
rz(-1.8476723) q[1];
sx q[1];
rz(-0.43101355) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9760808) q[0];
sx q[0];
rz(-1.3595194) q[0];
sx q[0];
rz(1.4552646) q[0];
x q[1];
rz(2.7796499) q[2];
sx q[2];
rz(-1.2133994) q[2];
sx q[2];
rz(1.5882284) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9412028) q[1];
sx q[1];
rz(-1.871849) q[1];
sx q[1];
rz(0.55834229) q[1];
rz(0.94915821) q[3];
sx q[3];
rz(-1.7931058) q[3];
sx q[3];
rz(-2.3883826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5203349) q[2];
sx q[2];
rz(-2.4020577) q[2];
sx q[2];
rz(0.4062824) q[2];
rz(-1.9351561) q[3];
sx q[3];
rz(-1.0932357) q[3];
sx q[3];
rz(-2.5549197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63221145) q[0];
sx q[0];
rz(-0.87566942) q[0];
sx q[0];
rz(-2.8152554) q[0];
rz(2.1874766) q[1];
sx q[1];
rz(-2.4931144) q[1];
sx q[1];
rz(3.1180678) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017880579) q[0];
sx q[0];
rz(-2.5766226) q[0];
sx q[0];
rz(-0.82161696) q[0];
x q[1];
rz(-0.40074375) q[2];
sx q[2];
rz(-0.6983499) q[2];
sx q[2];
rz(1.9633479) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.22241961) q[1];
sx q[1];
rz(-1.8559764) q[1];
sx q[1];
rz(1.213277) q[1];
x q[2];
rz(-0.29554772) q[3];
sx q[3];
rz(-1.7046527) q[3];
sx q[3];
rz(-0.16242151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.894459) q[2];
sx q[2];
rz(-2.6829312) q[2];
sx q[2];
rz(2.4776754) q[2];
rz(-0.89020056) q[3];
sx q[3];
rz(-1.592344) q[3];
sx q[3];
rz(0.012705407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7909872) q[0];
sx q[0];
rz(-0.43142879) q[0];
sx q[0];
rz(-2.8107693) q[0];
rz(0.36258969) q[1];
sx q[1];
rz(-1.9622842) q[1];
sx q[1];
rz(-2.6258452) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4048201) q[0];
sx q[0];
rz(-1.3546001) q[0];
sx q[0];
rz(-1.3275073) q[0];
x q[1];
rz(1.7605893) q[2];
sx q[2];
rz(-0.98514531) q[2];
sx q[2];
rz(0.3444258) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95553229) q[1];
sx q[1];
rz(-1.6481676) q[1];
sx q[1];
rz(-0.87826985) q[1];
rz(0.12966517) q[3];
sx q[3];
rz(-2.3414681) q[3];
sx q[3];
rz(-0.4365546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7911239) q[2];
sx q[2];
rz(-1.4501269) q[2];
sx q[2];
rz(-0.80751944) q[2];
rz(-3.0430072) q[3];
sx q[3];
rz(-2.2435296) q[3];
sx q[3];
rz(-1.0928104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7035141) q[0];
sx q[0];
rz(-0.57323891) q[0];
sx q[0];
rz(2.692063) q[0];
rz(0.685177) q[1];
sx q[1];
rz(-2.7339869) q[1];
sx q[1];
rz(2.5221241) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0280058) q[0];
sx q[0];
rz(-1.2613234) q[0];
sx q[0];
rz(1.600515) q[0];
rz(-3.0590942) q[2];
sx q[2];
rz(-0.52972016) q[2];
sx q[2];
rz(-2.8166383) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7799311) q[1];
sx q[1];
rz(-1.9059935) q[1];
sx q[1];
rz(0.72104309) q[1];
rz(0.18467091) q[3];
sx q[3];
rz(-1.369056) q[3];
sx q[3];
rz(1.7137485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5770136) q[2];
sx q[2];
rz(-1.6926293) q[2];
sx q[2];
rz(2.9435834) q[2];
rz(1.3624582) q[3];
sx q[3];
rz(-0.30006108) q[3];
sx q[3];
rz(0.70039606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.144416) q[0];
sx q[0];
rz(-2.5197025) q[0];
sx q[0];
rz(1.8560334) q[0];
rz(-0.27664912) q[1];
sx q[1];
rz(-1.3686907) q[1];
sx q[1];
rz(0.10861529) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1254107) q[0];
sx q[0];
rz(-1.4889532) q[0];
sx q[0];
rz(2.0731519) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1947961) q[2];
sx q[2];
rz(-1.5063926) q[2];
sx q[2];
rz(-1.4381806) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.020713652) q[1];
sx q[1];
rz(-0.87057897) q[1];
sx q[1];
rz(-0.28182272) q[1];
rz(-pi) q[2];
rz(1.8720988) q[3];
sx q[3];
rz(-0.40280324) q[3];
sx q[3];
rz(0.81217743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8331929) q[2];
sx q[2];
rz(-1.1966285) q[2];
sx q[2];
rz(1.2742554) q[2];
rz(-1.3537004) q[3];
sx q[3];
rz(-2.7022868) q[3];
sx q[3];
rz(0.86658365) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15247791) q[0];
sx q[0];
rz(-0.72565961) q[0];
sx q[0];
rz(3.0210378) q[0];
rz(-2.4001135) q[1];
sx q[1];
rz(-1.9375786) q[1];
sx q[1];
rz(-0.075686879) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78223373) q[0];
sx q[0];
rz(-1.4067211) q[0];
sx q[0];
rz(2.8316281) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7255177) q[2];
sx q[2];
rz(-0.46636367) q[2];
sx q[2];
rz(-1.5049962) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8676641) q[1];
sx q[1];
rz(-1.5455265) q[1];
sx q[1];
rz(-0.65622337) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.754154) q[3];
sx q[3];
rz(-0.9777189) q[3];
sx q[3];
rz(-1.0870522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0691444) q[2];
sx q[2];
rz(-0.88280237) q[2];
sx q[2];
rz(-0.34633386) q[2];
rz(2.196178) q[3];
sx q[3];
rz(-1.8190705) q[3];
sx q[3];
rz(2.1944428) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.866975) q[0];
sx q[0];
rz(-1.9434384) q[0];
sx q[0];
rz(-0.44788885) q[0];
rz(0.92944324) q[1];
sx q[1];
rz(-1.399469) q[1];
sx q[1];
rz(0.67451745) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4436627) q[0];
sx q[0];
rz(-0.81564689) q[0];
sx q[0];
rz(1.2944512) q[0];
rz(2.0777763) q[2];
sx q[2];
rz(-0.21179767) q[2];
sx q[2];
rz(-1.3174881) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.50215534) q[1];
sx q[1];
rz(-1.5719255) q[1];
sx q[1];
rz(-1.3974031) q[1];
rz(-3.0661656) q[3];
sx q[3];
rz(-1.1531246) q[3];
sx q[3];
rz(0.088069629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.931539) q[2];
sx q[2];
rz(-1.209582) q[2];
sx q[2];
rz(-0.26739576) q[2];
rz(2.0343272) q[3];
sx q[3];
rz(-0.78247726) q[3];
sx q[3];
rz(0.84899181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4103107) q[0];
sx q[0];
rz(-1.5189497) q[0];
sx q[0];
rz(2.0547163) q[0];
rz(-0.80401737) q[1];
sx q[1];
rz(-1.4878648) q[1];
sx q[1];
rz(1.0555242) q[1];
rz(1.4236535) q[2];
sx q[2];
rz(-1.6416807) q[2];
sx q[2];
rz(0.15105187) q[2];
rz(-1.6770398) q[3];
sx q[3];
rz(-0.67062702) q[3];
sx q[3];
rz(0.39556243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
