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
rz(0.1641195) q[0];
sx q[0];
rz(-2.1174705) q[0];
sx q[0];
rz(0.48409387) q[0];
rz(2.3789499) q[1];
sx q[1];
rz(4.7197309) q[1];
sx q[1];
rz(9.3698256) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23942102) q[0];
sx q[0];
rz(-2.4706984) q[0];
sx q[0];
rz(-1.5205199) q[0];
rz(0.058637549) q[2];
sx q[2];
rz(-1.3212886) q[2];
sx q[2];
rz(-0.80896689) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.311885) q[1];
sx q[1];
rz(-1.6685608) q[1];
sx q[1];
rz(0.13297653) q[1];
rz(-2.1922429) q[3];
sx q[3];
rz(-1.8481405) q[3];
sx q[3];
rz(-2.8327033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.83076465) q[2];
sx q[2];
rz(-0.97042933) q[2];
sx q[2];
rz(-0.27466276) q[2];
rz(2.2667387) q[3];
sx q[3];
rz(-2.0503876) q[3];
sx q[3];
rz(-2.8680475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65983588) q[0];
sx q[0];
rz(-2.4517224) q[0];
sx q[0];
rz(0.39875317) q[0];
rz(2.8975471) q[1];
sx q[1];
rz(-1.9052541) q[1];
sx q[1];
rz(1.1057378) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4717446) q[0];
sx q[0];
rz(-2.6378184) q[0];
sx q[0];
rz(-0.48724799) q[0];
rz(-pi) q[1];
rz(2.2695838) q[2];
sx q[2];
rz(-0.78561831) q[2];
sx q[2];
rz(0.18839041) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51197375) q[1];
sx q[1];
rz(-1.5675768) q[1];
sx q[1];
rz(-1.5681439) q[1];
rz(-1.8375754) q[3];
sx q[3];
rz(-1.6872129) q[3];
sx q[3];
rz(-2.4698225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4411321) q[2];
sx q[2];
rz(-1.9596142) q[2];
sx q[2];
rz(-2.0665118) q[2];
rz(-0.69617802) q[3];
sx q[3];
rz(-1.2735561) q[3];
sx q[3];
rz(-2.9785494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8657846) q[0];
sx q[0];
rz(-1.2514665) q[0];
sx q[0];
rz(-2.9248917) q[0];
rz(-1.7614583) q[1];
sx q[1];
rz(-2.7237027) q[1];
sx q[1];
rz(-0.61947852) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.769128) q[0];
sx q[0];
rz(-1.6004246) q[0];
sx q[0];
rz(1.6086786) q[0];
rz(-pi) q[1];
rz(3.0946629) q[2];
sx q[2];
rz(-1.6937733) q[2];
sx q[2];
rz(-2.8356981) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.463355) q[1];
sx q[1];
rz(-1.2690407) q[1];
sx q[1];
rz(-0.22308087) q[1];
rz(-pi) q[2];
rz(2.5314919) q[3];
sx q[3];
rz(-2.0629632) q[3];
sx q[3];
rz(2.6254693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2565903) q[2];
sx q[2];
rz(-1.2991178) q[2];
sx q[2];
rz(0.090506434) q[2];
rz(-2.6835486) q[3];
sx q[3];
rz(-0.48726714) q[3];
sx q[3];
rz(-2.0335782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3312382) q[0];
sx q[0];
rz(-2.2585223) q[0];
sx q[0];
rz(0.93165398) q[0];
rz(0.16904198) q[1];
sx q[1];
rz(-0.5642429) q[1];
sx q[1];
rz(1.0394675) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9477978) q[0];
sx q[0];
rz(-2.373824) q[0];
sx q[0];
rz(2.7916273) q[0];
rz(-1.3170856) q[2];
sx q[2];
rz(-1.3519962) q[2];
sx q[2];
rz(-1.4267061) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9586981) q[1];
sx q[1];
rz(-2.9258203) q[1];
sx q[1];
rz(-1.2879667) q[1];
rz(-pi) q[2];
rz(1.3576034) q[3];
sx q[3];
rz(-1.2925832) q[3];
sx q[3];
rz(-1.9813862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6197551) q[2];
sx q[2];
rz(-1.1014742) q[2];
sx q[2];
rz(-2.6452765) q[2];
rz(-0.79658341) q[3];
sx q[3];
rz(-0.28592548) q[3];
sx q[3];
rz(2.0485785) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26688823) q[0];
sx q[0];
rz(-0.58458352) q[0];
sx q[0];
rz(-0.97187483) q[0];
rz(3.019849) q[1];
sx q[1];
rz(-1.6106482) q[1];
sx q[1];
rz(-1.3294539) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15518256) q[0];
sx q[0];
rz(-1.2647795) q[0];
sx q[0];
rz(-3.0772692) q[0];
x q[1];
rz(0.8829863) q[2];
sx q[2];
rz(-0.60740439) q[2];
sx q[2];
rz(-2.4508053) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5962474) q[1];
sx q[1];
rz(-1.9831428) q[1];
sx q[1];
rz(2.418787) q[1];
rz(-1.6498783) q[3];
sx q[3];
rz(-1.5615765) q[3];
sx q[3];
rz(0.10445933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.10151265) q[2];
sx q[2];
rz(-1.2206581) q[2];
sx q[2];
rz(-0.15138781) q[2];
rz(-0.76465145) q[3];
sx q[3];
rz(-0.83573666) q[3];
sx q[3];
rz(1.5966655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0720035) q[0];
sx q[0];
rz(-1.66865) q[0];
sx q[0];
rz(-1.0714916) q[0];
rz(2.6103861) q[1];
sx q[1];
rz(-1.4899645) q[1];
sx q[1];
rz(-2.5154617) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75979489) q[0];
sx q[0];
rz(-3.1239428) q[0];
sx q[0];
rz(-1.2977029) q[0];
rz(-pi) q[1];
rz(-1.3254204) q[2];
sx q[2];
rz(-0.98773709) q[2];
sx q[2];
rz(-0.32223216) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5007361) q[1];
sx q[1];
rz(-1.2656414) q[1];
sx q[1];
rz(2.5049097) q[1];
x q[2];
rz(-0.71685426) q[3];
sx q[3];
rz(-0.62590137) q[3];
sx q[3];
rz(0.38960534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35817394) q[2];
sx q[2];
rz(-0.61966115) q[2];
sx q[2];
rz(2.1273071) q[2];
rz(1.8861534) q[3];
sx q[3];
rz(-1.3557326) q[3];
sx q[3];
rz(-2.0881418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96100539) q[0];
sx q[0];
rz(-2.2659681) q[0];
sx q[0];
rz(0.92051202) q[0];
rz(-3.0275184) q[1];
sx q[1];
rz(-1.4055077) q[1];
sx q[1];
rz(-1.1309518) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3217472) q[0];
sx q[0];
rz(-1.5764569) q[0];
sx q[0];
rz(-0.59771718) q[0];
x q[1];
rz(-0.99920706) q[2];
sx q[2];
rz(-1.0896557) q[2];
sx q[2];
rz(-2.4740117) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9662094) q[1];
sx q[1];
rz(-1.7031324) q[1];
sx q[1];
rz(2.2369409) q[1];
rz(0.5163631) q[3];
sx q[3];
rz(-1.6404021) q[3];
sx q[3];
rz(-2.7854837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34746927) q[2];
sx q[2];
rz(-2.4994714) q[2];
sx q[2];
rz(2.7327909) q[2];
rz(2.9583904) q[3];
sx q[3];
rz(-1.5902404) q[3];
sx q[3];
rz(2.0028152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0278397) q[0];
sx q[0];
rz(-3.0945859) q[0];
sx q[0];
rz(-0.29079944) q[0];
rz(-2.193702) q[1];
sx q[1];
rz(-0.46698505) q[1];
sx q[1];
rz(1.9416169) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3481846) q[0];
sx q[0];
rz(-0.68773341) q[0];
sx q[0];
rz(-2.4853112) q[0];
rz(-0.84917111) q[2];
sx q[2];
rz(-2.2361922) q[2];
sx q[2];
rz(-0.1553436) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.34741286) q[1];
sx q[1];
rz(-2.844226) q[1];
sx q[1];
rz(1.7296687) q[1];
x q[2];
rz(1.387595) q[3];
sx q[3];
rz(-1.6699413) q[3];
sx q[3];
rz(-2.9354036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36134186) q[2];
sx q[2];
rz(-1.6929408) q[2];
sx q[2];
rz(-1.348749) q[2];
rz(-1.6382943) q[3];
sx q[3];
rz(-0.88204757) q[3];
sx q[3];
rz(-2.9564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1143484) q[0];
sx q[0];
rz(-2.769727) q[0];
sx q[0];
rz(-1.0741023) q[0];
rz(0.4298003) q[1];
sx q[1];
rz(-2.0975515) q[1];
sx q[1];
rz(0.46357402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4922475) q[0];
sx q[0];
rz(-1.6657636) q[0];
sx q[0];
rz(-1.7596721) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26185449) q[2];
sx q[2];
rz(-1.7657649) q[2];
sx q[2];
rz(1.5443813) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8817655) q[1];
sx q[1];
rz(-0.9503839) q[1];
sx q[1];
rz(-0.46807351) q[1];
rz(-pi) q[2];
rz(-2.5829599) q[3];
sx q[3];
rz(-1.4945293) q[3];
sx q[3];
rz(2.7879451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8672436) q[2];
sx q[2];
rz(-2.5090019) q[2];
sx q[2];
rz(-0.80844936) q[2];
rz(-2.6570053) q[3];
sx q[3];
rz(-2.7633568) q[3];
sx q[3];
rz(2.7073879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.7782068) q[0];
sx q[0];
rz(-2.6860542) q[0];
sx q[0];
rz(-1.1325915) q[0];
rz(0.038381902) q[1];
sx q[1];
rz(-1.8033586) q[1];
sx q[1];
rz(-2.8448232) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0050186) q[0];
sx q[0];
rz(-1.8474192) q[0];
sx q[0];
rz(-2.121595) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1266668) q[2];
sx q[2];
rz(-1.9812366) q[2];
sx q[2];
rz(3.0164312) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6925507) q[1];
sx q[1];
rz(-1.3490648) q[1];
sx q[1];
rz(-1.6479083) q[1];
rz(-pi) q[2];
rz(2.2567883) q[3];
sx q[3];
rz(-1.8821041) q[3];
sx q[3];
rz(1.5112359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9639637) q[2];
sx q[2];
rz(-2.0202426) q[2];
sx q[2];
rz(-0.56619823) q[2];
rz(-1.1244134) q[3];
sx q[3];
rz(-3.0163613) q[3];
sx q[3];
rz(-1.9521149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88846702) q[0];
sx q[0];
rz(-2.4430226) q[0];
sx q[0];
rz(-3.0342614) q[0];
rz(-0.26168564) q[1];
sx q[1];
rz(-1.9410004) q[1];
sx q[1];
rz(0.95071361) q[1];
rz(0.28957257) q[2];
sx q[2];
rz(-0.085554335) q[2];
sx q[2];
rz(2.6680996) q[2];
rz(0.64601267) q[3];
sx q[3];
rz(-2.1727242) q[3];
sx q[3];
rz(0.77632191) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
