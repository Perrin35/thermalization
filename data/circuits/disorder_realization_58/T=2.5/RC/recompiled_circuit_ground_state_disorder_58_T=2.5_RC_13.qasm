OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6405606) q[0];
sx q[0];
rz(-0.80067331) q[0];
sx q[0];
rz(2.9963357) q[0];
rz(2.2486806) q[1];
sx q[1];
rz(-0.4141663) q[1];
sx q[1];
rz(-1.1069586) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1705748) q[0];
sx q[0];
rz(-0.98890328) q[0];
sx q[0];
rz(-2.2547743) q[0];
rz(-1.2453503) q[2];
sx q[2];
rz(-1.3017442) q[2];
sx q[2];
rz(-0.01435752) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.332079) q[1];
sx q[1];
rz(-2.4831536) q[1];
sx q[1];
rz(1.0926985) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1916394) q[3];
sx q[3];
rz(-1.3690476) q[3];
sx q[3];
rz(-2.4011322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5408111) q[2];
sx q[2];
rz(-1.3593707) q[2];
sx q[2];
rz(3.0493128) q[2];
rz(1.4394834) q[3];
sx q[3];
rz(-1.1441792) q[3];
sx q[3];
rz(0.93851844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1321024) q[0];
sx q[0];
rz(-1.9785896) q[0];
sx q[0];
rz(-2.5376885) q[0];
rz(1.5314792) q[1];
sx q[1];
rz(-1.3319301) q[1];
sx q[1];
rz(-1.6139222) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22900362) q[0];
sx q[0];
rz(-0.67978501) q[0];
sx q[0];
rz(2.0849821) q[0];
x q[1];
rz(2.1555734) q[2];
sx q[2];
rz(-0.9169609) q[2];
sx q[2];
rz(-2.3399835) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2774573) q[1];
sx q[1];
rz(-1.8445332) q[1];
sx q[1];
rz(-1.5828733) q[1];
rz(3.0433995) q[3];
sx q[3];
rz(-1.5714688) q[3];
sx q[3];
rz(2.0405586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4332726) q[2];
sx q[2];
rz(-1.4488139) q[2];
sx q[2];
rz(1.0478919) q[2];
rz(-0.6692872) q[3];
sx q[3];
rz(-0.71988121) q[3];
sx q[3];
rz(-1.9765114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.79353756) q[0];
sx q[0];
rz(-2.3568643) q[0];
sx q[0];
rz(1.0070356) q[0];
rz(-2.8496565) q[1];
sx q[1];
rz(-0.54436362) q[1];
sx q[1];
rz(0.67063355) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8438727) q[0];
sx q[0];
rz(-2.5150635) q[0];
sx q[0];
rz(-1.7202176) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2569855) q[2];
sx q[2];
rz(-0.90145196) q[2];
sx q[2];
rz(1.455292) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.933316) q[1];
sx q[1];
rz(-1.9066992) q[1];
sx q[1];
rz(-0.64888727) q[1];
rz(1.2420044) q[3];
sx q[3];
rz(-0.34600779) q[3];
sx q[3];
rz(3.0228928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2900419) q[2];
sx q[2];
rz(-2.2883577) q[2];
sx q[2];
rz(2.2500989) q[2];
rz(-2.0391035) q[3];
sx q[3];
rz(-2.1474371) q[3];
sx q[3];
rz(1.7360784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82146984) q[0];
sx q[0];
rz(-1.7550884) q[0];
sx q[0];
rz(1.0572877) q[0];
rz(-3.140246) q[1];
sx q[1];
rz(-2.3206382) q[1];
sx q[1];
rz(-0.79426208) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0214572) q[0];
sx q[0];
rz(-1.5524143) q[0];
sx q[0];
rz(-0.6619428) q[0];
x q[1];
rz(-2.7873331) q[2];
sx q[2];
rz(-1.2984167) q[2];
sx q[2];
rz(-0.59857063) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0483305) q[1];
sx q[1];
rz(-0.89846134) q[1];
sx q[1];
rz(0.59497084) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36093386) q[3];
sx q[3];
rz(-0.31692255) q[3];
sx q[3];
rz(1.8529056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0346251) q[2];
sx q[2];
rz(-2.5204973) q[2];
sx q[2];
rz(-2.1194439) q[2];
rz(-1.8978097) q[3];
sx q[3];
rz(-2.1918178) q[3];
sx q[3];
rz(-1.3589842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6720402) q[0];
sx q[0];
rz(-1.38009) q[0];
sx q[0];
rz(0.45355466) q[0];
rz(2.1039311) q[1];
sx q[1];
rz(-2.0062168) q[1];
sx q[1];
rz(-0.79197788) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4290473) q[0];
sx q[0];
rz(-0.28158108) q[0];
sx q[0];
rz(1.7292119) q[0];
rz(-pi) q[1];
rz(-0.99292314) q[2];
sx q[2];
rz(-1.2348334) q[2];
sx q[2];
rz(3.1022037) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5752069) q[1];
sx q[1];
rz(-2.7820911) q[1];
sx q[1];
rz(0.33728894) q[1];
x q[2];
rz(-1.4384934) q[3];
sx q[3];
rz(-3.0022394) q[3];
sx q[3];
rz(-0.15755586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8301293) q[2];
sx q[2];
rz(-2.4757803) q[2];
sx q[2];
rz(0.30544454) q[2];
rz(0.18181248) q[3];
sx q[3];
rz(-1.5287377) q[3];
sx q[3];
rz(2.4055433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55564725) q[0];
sx q[0];
rz(-0.82252994) q[0];
sx q[0];
rz(0.12538759) q[0];
rz(1.5665945) q[1];
sx q[1];
rz(-1.4488723) q[1];
sx q[1];
rz(-3.1094508) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9956995) q[0];
sx q[0];
rz(-0.75135485) q[0];
sx q[0];
rz(-2.8180647) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3022167) q[2];
sx q[2];
rz(-1.5245066) q[2];
sx q[2];
rz(2.3037825) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5405724) q[1];
sx q[1];
rz(-2.2443612) q[1];
sx q[1];
rz(-2.5190398) q[1];
x q[2];
rz(-1.4549082) q[3];
sx q[3];
rz(-0.82852302) q[3];
sx q[3];
rz(-2.7477086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0142168) q[2];
sx q[2];
rz(-1.2140423) q[2];
sx q[2];
rz(-2.0188913) q[2];
rz(0.073401062) q[3];
sx q[3];
rz(-0.95728907) q[3];
sx q[3];
rz(0.61298031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70925322) q[0];
sx q[0];
rz(-1.4839577) q[0];
sx q[0];
rz(0.51026979) q[0];
rz(2.7773652) q[1];
sx q[1];
rz(-2.7204456) q[1];
sx q[1];
rz(-1.4998923) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7059874) q[0];
sx q[0];
rz(-1.6765119) q[0];
sx q[0];
rz(-1.1054429) q[0];
rz(1.7287615) q[2];
sx q[2];
rz(-1.6876432) q[2];
sx q[2];
rz(0.37010461) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8270696) q[1];
sx q[1];
rz(-2.0644958) q[1];
sx q[1];
rz(-1.3160454) q[1];
rz(-pi) q[2];
rz(-2.4046201) q[3];
sx q[3];
rz(-1.342257) q[3];
sx q[3];
rz(1.4481973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4437272) q[2];
sx q[2];
rz(-1.5865734) q[2];
sx q[2];
rz(0.86722428) q[2];
rz(1.5813658) q[3];
sx q[3];
rz(-2.9428704) q[3];
sx q[3];
rz(-1.8038512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7241868) q[0];
sx q[0];
rz(-0.066411821) q[0];
sx q[0];
rz(-2.7253286) q[0];
rz(1.9789713) q[1];
sx q[1];
rz(-1.1888209) q[1];
sx q[1];
rz(-2.3847041) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74706057) q[0];
sx q[0];
rz(-1.9950657) q[0];
sx q[0];
rz(0.31655689) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95332884) q[2];
sx q[2];
rz(-2.4673415) q[2];
sx q[2];
rz(1.9566388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1637266) q[1];
sx q[1];
rz(-1.9975348) q[1];
sx q[1];
rz(-2.9336998) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0361597) q[3];
sx q[3];
rz(-0.98794395) q[3];
sx q[3];
rz(-2.7378156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4550712) q[2];
sx q[2];
rz(-1.7073809) q[2];
sx q[2];
rz(1.3580458) q[2];
rz(-0.81651917) q[3];
sx q[3];
rz(-1.2314545) q[3];
sx q[3];
rz(-2.9787298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4002976) q[0];
sx q[0];
rz(-1.6930027) q[0];
sx q[0];
rz(1.4073538) q[0];
rz(-0.74527144) q[1];
sx q[1];
rz(-0.73423568) q[1];
sx q[1];
rz(1.5819246) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3634264) q[0];
sx q[0];
rz(-0.72371783) q[0];
sx q[0];
rz(1.9947718) q[0];
x q[1];
rz(-1.5383155) q[2];
sx q[2];
rz(-0.92327061) q[2];
sx q[2];
rz(-0.10570174) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.386454) q[1];
sx q[1];
rz(-1.6999131) q[1];
sx q[1];
rz(2.0172202) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6635062) q[3];
sx q[3];
rz(-0.57766908) q[3];
sx q[3];
rz(-2.3232164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4173296) q[2];
sx q[2];
rz(-1.4549007) q[2];
sx q[2];
rz(2.0261185) q[2];
rz(-2.8880902) q[3];
sx q[3];
rz(-1.8799672) q[3];
sx q[3];
rz(-0.0089664627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1215006) q[0];
sx q[0];
rz(-1.2637063) q[0];
sx q[0];
rz(-0.76250917) q[0];
rz(0.94238845) q[1];
sx q[1];
rz(-1.0647048) q[1];
sx q[1];
rz(-1.9047838) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4651637) q[0];
sx q[0];
rz(-1.0753514) q[0];
sx q[0];
rz(-1.6406607) q[0];
rz(-pi) q[1];
rz(2.2549596) q[2];
sx q[2];
rz(-1.1012474) q[2];
sx q[2];
rz(-0.70882112) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33442228) q[1];
sx q[1];
rz(-1.5457193) q[1];
sx q[1];
rz(-1.5958171) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0701261) q[3];
sx q[3];
rz(-1.9744919) q[3];
sx q[3];
rz(1.5675327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2463871) q[2];
sx q[2];
rz(-0.095194101) q[2];
sx q[2];
rz(0.0072366317) q[2];
rz(0.17035189) q[3];
sx q[3];
rz(-1.4879613) q[3];
sx q[3];
rz(-2.4836704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.2687179) q[0];
sx q[0];
rz(-1.7418516) q[0];
sx q[0];
rz(-2.0280784) q[0];
rz(3.0524104) q[1];
sx q[1];
rz(-1.6249648) q[1];
sx q[1];
rz(-1.9393495) q[1];
rz(0.7128255) q[2];
sx q[2];
rz(-0.20607866) q[2];
sx q[2];
rz(-1.1417749) q[2];
rz(-0.30059697) q[3];
sx q[3];
rz(-1.3924122) q[3];
sx q[3];
rz(-2.198749) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
