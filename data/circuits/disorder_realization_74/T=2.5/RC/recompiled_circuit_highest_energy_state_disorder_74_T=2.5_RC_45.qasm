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
rz(0.81646252) q[0];
sx q[0];
rz(-3.0397968) q[0];
sx q[0];
rz(0.53959227) q[0];
rz(0.49630961) q[1];
sx q[1];
rz(-0.30975431) q[1];
sx q[1];
rz(-2.6024979) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70654725) q[0];
sx q[0];
rz(-0.4420949) q[0];
sx q[0];
rz(-3.0821783) q[0];
rz(-2.7560948) q[2];
sx q[2];
rz(-2.756167) q[2];
sx q[2];
rz(2.5138234) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8546363) q[1];
sx q[1];
rz(-2.7883734) q[1];
sx q[1];
rz(1.9853206) q[1];
x q[2];
rz(2.21978) q[3];
sx q[3];
rz(-1.6361437) q[3];
sx q[3];
rz(1.1999038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3794136) q[2];
sx q[2];
rz(-0.87612408) q[2];
sx q[2];
rz(-1.9525105) q[2];
rz(-1.9573697) q[3];
sx q[3];
rz(-2.2370179) q[3];
sx q[3];
rz(1.5466461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14520833) q[0];
sx q[0];
rz(-1.5518016) q[0];
sx q[0];
rz(-1.8183964) q[0];
rz(-0.48201758) q[1];
sx q[1];
rz(-0.91164416) q[1];
sx q[1];
rz(0.97420305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7071814) q[0];
sx q[0];
rz(-1.7493346) q[0];
sx q[0];
rz(-2.2212127) q[0];
rz(2.2902238) q[2];
sx q[2];
rz(-2.9915644) q[2];
sx q[2];
rz(2.0025557) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7108591) q[1];
sx q[1];
rz(-1.3361592) q[1];
sx q[1];
rz(0.55107848) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5752546) q[3];
sx q[3];
rz(-1.3674595) q[3];
sx q[3];
rz(-1.4161863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0049858967) q[2];
sx q[2];
rz(-0.58711457) q[2];
sx q[2];
rz(-0.74964398) q[2];
rz(-0.43198112) q[3];
sx q[3];
rz(-1.1155198) q[3];
sx q[3];
rz(-1.8384793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62464803) q[0];
sx q[0];
rz(-0.2722781) q[0];
sx q[0];
rz(-0.6063478) q[0];
rz(1.3336522) q[1];
sx q[1];
rz(-1.2140112) q[1];
sx q[1];
rz(2.8820754) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57588803) q[0];
sx q[0];
rz(-1.7262926) q[0];
sx q[0];
rz(-1.3208273) q[0];
rz(1.7735931) q[2];
sx q[2];
rz(-1.3992056) q[2];
sx q[2];
rz(-2.9362048) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1705679) q[1];
sx q[1];
rz(-1.9924506) q[1];
sx q[1];
rz(-2.7287911) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2058425) q[3];
sx q[3];
rz(-0.66422909) q[3];
sx q[3];
rz(2.2838044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8848662) q[2];
sx q[2];
rz(-1.7776411) q[2];
sx q[2];
rz(1.5941031) q[2];
rz(-0.78440845) q[3];
sx q[3];
rz(-1.6268077) q[3];
sx q[3];
rz(1.5756395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49482685) q[0];
sx q[0];
rz(-2.0443199) q[0];
sx q[0];
rz(-2.8821017) q[0];
rz(1.9208113) q[1];
sx q[1];
rz(-2.6916598) q[1];
sx q[1];
rz(-1.4422013) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3078559) q[0];
sx q[0];
rz(-0.74001946) q[0];
sx q[0];
rz(-0.83777512) q[0];
x q[1];
rz(-3.0109826) q[2];
sx q[2];
rz(-2.161986) q[2];
sx q[2];
rz(-1.3863877) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3825697) q[1];
sx q[1];
rz(-1.8174531) q[1];
sx q[1];
rz(-0.36212977) q[1];
rz(-pi) q[2];
rz(0.070017858) q[3];
sx q[3];
rz(-2.2517831) q[3];
sx q[3];
rz(2.3511166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0356902) q[2];
sx q[2];
rz(-1.8041958) q[2];
sx q[2];
rz(2.7316459) q[2];
rz(1.9299054) q[3];
sx q[3];
rz(-2.4544921) q[3];
sx q[3];
rz(-1.80779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.425151) q[0];
sx q[0];
rz(-2.4254159) q[0];
sx q[0];
rz(2.8675365) q[0];
rz(-2.7033499) q[1];
sx q[1];
rz(-1.2934877) q[1];
sx q[1];
rz(-0.73572198) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5490131) q[0];
sx q[0];
rz(-1.4837449) q[0];
sx q[0];
rz(-2.3866231) q[0];
x q[1];
rz(0.76916285) q[2];
sx q[2];
rz(-2.935098) q[2];
sx q[2];
rz(-2.1392876) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3408518) q[1];
sx q[1];
rz(-1.2781118) q[1];
sx q[1];
rz(2.8814948) q[1];
rz(2.208136) q[3];
sx q[3];
rz(-2.3846755) q[3];
sx q[3];
rz(0.33801038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.938544) q[2];
sx q[2];
rz(-1.6837348) q[2];
sx q[2];
rz(2.5277444) q[2];
rz(0.40890536) q[3];
sx q[3];
rz(-0.96858612) q[3];
sx q[3];
rz(-2.3992505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78855377) q[0];
sx q[0];
rz(-2.5034294) q[0];
sx q[0];
rz(1.9045389) q[0];
rz(1.6963814) q[1];
sx q[1];
rz(-0.45672363) q[1];
sx q[1];
rz(2.9663185) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1244558) q[0];
sx q[0];
rz(-1.5958028) q[0];
sx q[0];
rz(-1.7437922) q[0];
x q[1];
rz(-2.4779123) q[2];
sx q[2];
rz(-0.77518565) q[2];
sx q[2];
rz(-0.41942393) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8126909) q[1];
sx q[1];
rz(-1.1393424) q[1];
sx q[1];
rz(0.56568362) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95703362) q[3];
sx q[3];
rz(-2.0599457) q[3];
sx q[3];
rz(-3.1135984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1767629) q[2];
sx q[2];
rz(-1.8076597) q[2];
sx q[2];
rz(0.97243398) q[2];
rz(1.4771627) q[3];
sx q[3];
rz(-1.9921314) q[3];
sx q[3];
rz(-2.3798063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2870188) q[0];
sx q[0];
rz(-0.022495689) q[0];
sx q[0];
rz(0.35312411) q[0];
rz(1.0278206) q[1];
sx q[1];
rz(-1.9408344) q[1];
sx q[1];
rz(2.1305398) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9927427) q[0];
sx q[0];
rz(-1.5563838) q[0];
sx q[0];
rz(0.3573125) q[0];
rz(-pi) q[1];
rz(-0.6017466) q[2];
sx q[2];
rz(-2.1831552) q[2];
sx q[2];
rz(-1.5073206) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.331703) q[1];
sx q[1];
rz(-1.9192438) q[1];
sx q[1];
rz(2.9747559) q[1];
rz(2.3156741) q[3];
sx q[3];
rz(-2.7248235) q[3];
sx q[3];
rz(-0.59405223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.31356835) q[2];
sx q[2];
rz(-2.821533) q[2];
sx q[2];
rz(-1.3671406) q[2];
rz(1.2683055) q[3];
sx q[3];
rz(-1.6627848) q[3];
sx q[3];
rz(-0.84793276) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0276412) q[0];
sx q[0];
rz(-2.5340762) q[0];
sx q[0];
rz(1.129958) q[0];
rz(-2.9226411) q[1];
sx q[1];
rz(-1.7349225) q[1];
sx q[1];
rz(0.91046441) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8559056) q[0];
sx q[0];
rz(-2.160851) q[0];
sx q[0];
rz(-1.9737712) q[0];
x q[1];
rz(1.4850281) q[2];
sx q[2];
rz(-1.0940486) q[2];
sx q[2];
rz(-1.3439182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3748951) q[1];
sx q[1];
rz(-1.8378705) q[1];
sx q[1];
rz(-2.1007936) q[1];
x q[2];
rz(0.45121737) q[3];
sx q[3];
rz(-0.77773753) q[3];
sx q[3];
rz(1.3743708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8255446) q[2];
sx q[2];
rz(-2.3393708) q[2];
sx q[2];
rz(0.82553378) q[2];
rz(2.6801706) q[3];
sx q[3];
rz(-1.2666707) q[3];
sx q[3];
rz(-0.44574827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6398741) q[0];
sx q[0];
rz(-1.3383144) q[0];
sx q[0];
rz(-2.3097532) q[0];
rz(2.4141451) q[1];
sx q[1];
rz(-1.4201545) q[1];
sx q[1];
rz(3.0453392) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1883586) q[0];
sx q[0];
rz(-2.0287477) q[0];
sx q[0];
rz(-1.0012549) q[0];
rz(0.25419323) q[2];
sx q[2];
rz(-1.5113748) q[2];
sx q[2];
rz(2.3583902) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.61273328) q[1];
sx q[1];
rz(-1.7816356) q[1];
sx q[1];
rz(-0.092540578) q[1];
x q[2];
rz(-1.3745802) q[3];
sx q[3];
rz(-1.4150146) q[3];
sx q[3];
rz(2.2957612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7353797) q[2];
sx q[2];
rz(-1.7071743) q[2];
sx q[2];
rz(1.5926788) q[2];
rz(-0.63273543) q[3];
sx q[3];
rz(-1.9097208) q[3];
sx q[3];
rz(-0.24937853) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(11/(9*pi)) q[0];
sx q[0];
rz(-1.0405552) q[0];
sx q[0];
rz(-2.7885875) q[0];
rz(0.32265916) q[1];
sx q[1];
rz(-2.8162075) q[1];
sx q[1];
rz(-2.7696612) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7854561) q[0];
sx q[0];
rz(-1.8710941) q[0];
sx q[0];
rz(2.9343283) q[0];
rz(-pi) q[1];
rz(0.48666059) q[2];
sx q[2];
rz(-0.36937215) q[2];
sx q[2];
rz(-0.83115679) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1052211) q[1];
sx q[1];
rz(-2.7694355) q[1];
sx q[1];
rz(-0.45366617) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5266034) q[3];
sx q[3];
rz(-2.4165476) q[3];
sx q[3];
rz(2.8738662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.17452621) q[2];
sx q[2];
rz(-1.1897949) q[2];
sx q[2];
rz(1.2971499) q[2];
rz(2.2938812) q[3];
sx q[3];
rz(-1.1573557) q[3];
sx q[3];
rz(-2.1367836) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90947718) q[0];
sx q[0];
rz(-1.3790601) q[0];
sx q[0];
rz(-2.7098304) q[0];
rz(1.3879981) q[1];
sx q[1];
rz(-1.9276062) q[1];
sx q[1];
rz(3.066317) q[1];
rz(-2.3222011) q[2];
sx q[2];
rz(-2.1234305) q[2];
sx q[2];
rz(2.7101868) q[2];
rz(1.9907436) q[3];
sx q[3];
rz(-0.81880488) q[3];
sx q[3];
rz(2.8863751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
