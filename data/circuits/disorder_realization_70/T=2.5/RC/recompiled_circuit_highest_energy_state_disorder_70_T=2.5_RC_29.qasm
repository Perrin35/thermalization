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
rz(-0.47082585) q[0];
sx q[0];
rz(3.5916632) q[0];
sx q[0];
rz(9.8150742) q[0];
rz(-2.9974239) q[1];
sx q[1];
rz(-1.6427957) q[1];
sx q[1];
rz(1.0402476) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9269104) q[0];
sx q[0];
rz(-1.2217064) q[0];
sx q[0];
rz(0.30544282) q[0];
x q[1];
rz(-1.7261271) q[2];
sx q[2];
rz(-1.2110146) q[2];
sx q[2];
rz(1.8632165) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39506787) q[1];
sx q[1];
rz(-1.7774425) q[1];
sx q[1];
rz(2.2499491) q[1];
x q[2];
rz(3.0754117) q[3];
sx q[3];
rz(-2.3621763) q[3];
sx q[3];
rz(-0.73252892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1990004) q[2];
sx q[2];
rz(-0.60574836) q[2];
sx q[2];
rz(-1.595363) q[2];
rz(2.6664901) q[3];
sx q[3];
rz(-0.64517704) q[3];
sx q[3];
rz(-1.1554385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5763181) q[0];
sx q[0];
rz(-2.1529038) q[0];
sx q[0];
rz(-0.43462547) q[0];
rz(-2.077153) q[1];
sx q[1];
rz(-0.39355215) q[1];
sx q[1];
rz(0.89148608) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1284244) q[0];
sx q[0];
rz(-1.1718318) q[0];
sx q[0];
rz(-1.1815039) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72618809) q[2];
sx q[2];
rz(-1.1788158) q[2];
sx q[2];
rz(-1.7418944) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.75900092) q[1];
sx q[1];
rz(-2.0927974) q[1];
sx q[1];
rz(-2.0981789) q[1];
x q[2];
rz(0.029835506) q[3];
sx q[3];
rz(-0.91174698) q[3];
sx q[3];
rz(1.713879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6330304) q[2];
sx q[2];
rz(-0.56266251) q[2];
sx q[2];
rz(1.0594581) q[2];
rz(1.9153473) q[3];
sx q[3];
rz(-1.5303333) q[3];
sx q[3];
rz(0.1196158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0934963) q[0];
sx q[0];
rz(-0.39830783) q[0];
sx q[0];
rz(2.4523729) q[0];
rz(-2.7293909) q[1];
sx q[1];
rz(-2.0239794) q[1];
sx q[1];
rz(1.162792) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1173232) q[0];
sx q[0];
rz(-1.4179967) q[0];
sx q[0];
rz(-1.2627748) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77670375) q[2];
sx q[2];
rz(-1.113184) q[2];
sx q[2];
rz(2.6991778) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.38123576) q[1];
sx q[1];
rz(-2.3489526) q[1];
sx q[1];
rz(0.70602487) q[1];
x q[2];
rz(1.8602636) q[3];
sx q[3];
rz(-1.839387) q[3];
sx q[3];
rz(-0.1381607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.73843655) q[2];
sx q[2];
rz(-1.4102035) q[2];
sx q[2];
rz(2.9031244) q[2];
rz(-0.23369914) q[3];
sx q[3];
rz(-0.77478474) q[3];
sx q[3];
rz(0.15271798) q[3];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23976633) q[0];
sx q[0];
rz(-2.9775743) q[0];
sx q[0];
rz(0.30583403) q[0];
rz(0.26890525) q[1];
sx q[1];
rz(-2.3482359) q[1];
sx q[1];
rz(-2.7893524) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3650318) q[0];
sx q[0];
rz(-1.4370586) q[0];
sx q[0];
rz(-3.1015009) q[0];
rz(-pi) q[1];
rz(1.3046493) q[2];
sx q[2];
rz(-0.64322844) q[2];
sx q[2];
rz(0.51189724) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.70275408) q[1];
sx q[1];
rz(-1.9620336) q[1];
sx q[1];
rz(-2.4719098) q[1];
rz(-2.3533559) q[3];
sx q[3];
rz(-1.7468037) q[3];
sx q[3];
rz(0.22135425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8504146) q[2];
sx q[2];
rz(-1.6213657) q[2];
sx q[2];
rz(2.9746383) q[2];
rz(-0.11635612) q[3];
sx q[3];
rz(-2.6226624) q[3];
sx q[3];
rz(1.8714582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51036924) q[0];
sx q[0];
rz(-0.33741697) q[0];
sx q[0];
rz(-2.0262729) q[0];
rz(-1.2315617) q[1];
sx q[1];
rz(-1.7111338) q[1];
sx q[1];
rz(-3.0273052) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60576263) q[0];
sx q[0];
rz(-1.0530942) q[0];
sx q[0];
rz(0.96808956) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8631526) q[2];
sx q[2];
rz(-1.8634708) q[2];
sx q[2];
rz(-3.1029683) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.28262269) q[1];
sx q[1];
rz(-0.77853528) q[1];
sx q[1];
rz(0.60311879) q[1];
rz(-2.435569) q[3];
sx q[3];
rz(-1.5665652) q[3];
sx q[3];
rz(-0.26102456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.085122434) q[2];
sx q[2];
rz(-1.2135999) q[2];
sx q[2];
rz(-1.4465793) q[2];
rz(-2.1794686) q[3];
sx q[3];
rz(-2.646793) q[3];
sx q[3];
rz(1.8112633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1591448) q[0];
sx q[0];
rz(-1.5788989) q[0];
sx q[0];
rz(2.6084117) q[0];
rz(-1.606733) q[1];
sx q[1];
rz(-2.3695562) q[1];
sx q[1];
rz(-0.74347043) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6966382) q[0];
sx q[0];
rz(-1.6499203) q[0];
sx q[0];
rz(2.8487569) q[0];
x q[1];
rz(-1.1986046) q[2];
sx q[2];
rz(-0.57665885) q[2];
sx q[2];
rz(-1.3194336) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9333794) q[1];
sx q[1];
rz(-0.68720308) q[1];
sx q[1];
rz(-0.40172462) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0251158) q[3];
sx q[3];
rz(-0.66963306) q[3];
sx q[3];
rz(-2.7463934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.47034904) q[2];
sx q[2];
rz(-2.2942784) q[2];
sx q[2];
rz(0.49873763) q[2];
rz(0.70164743) q[3];
sx q[3];
rz(-2.4528153) q[3];
sx q[3];
rz(-0.65830314) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7094803) q[0];
sx q[0];
rz(-1.9146336) q[0];
sx q[0];
rz(3.0416601) q[0];
rz(-1.8607032) q[1];
sx q[1];
rz(-2.62968) q[1];
sx q[1];
rz(2.4023712) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4866132) q[0];
sx q[0];
rz(-0.57679048) q[0];
sx q[0];
rz(-0.77107112) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6373424) q[2];
sx q[2];
rz(-0.4302667) q[2];
sx q[2];
rz(0.86003785) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9770551) q[1];
sx q[1];
rz(-0.73325578) q[1];
sx q[1];
rz(-0.96340839) q[1];
x q[2];
rz(-1.5815063) q[3];
sx q[3];
rz(-2.3751343) q[3];
sx q[3];
rz(2.2720154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7382536) q[2];
sx q[2];
rz(-0.6311987) q[2];
sx q[2];
rz(-1.5367907) q[2];
rz(1.712045) q[3];
sx q[3];
rz(-1.4934544) q[3];
sx q[3];
rz(0.79609377) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1256063) q[0];
sx q[0];
rz(-0.37785372) q[0];
sx q[0];
rz(1.1146389) q[0];
rz(0.58427018) q[1];
sx q[1];
rz(-0.92002267) q[1];
sx q[1];
rz(0.11411962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7540365) q[0];
sx q[0];
rz(-0.68337959) q[0];
sx q[0];
rz(-2.6793733) q[0];
rz(-pi) q[1];
rz(0.80107208) q[2];
sx q[2];
rz(-1.1566887) q[2];
sx q[2];
rz(-1.0528477) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2456667) q[1];
sx q[1];
rz(-1.9787496) q[1];
sx q[1];
rz(2.308564) q[1];
rz(1.1646284) q[3];
sx q[3];
rz(-0.28985786) q[3];
sx q[3];
rz(1.690762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.70322651) q[2];
sx q[2];
rz(-0.72202903) q[2];
sx q[2];
rz(0.88877338) q[2];
rz(-1.599954) q[3];
sx q[3];
rz(-1.2069353) q[3];
sx q[3];
rz(2.3204939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0515161) q[0];
sx q[0];
rz(-0.30802825) q[0];
sx q[0];
rz(-2.8544881) q[0];
rz(1.2039394) q[1];
sx q[1];
rz(-0.86910373) q[1];
sx q[1];
rz(1.6285816) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8912761) q[0];
sx q[0];
rz(-0.34252942) q[0];
sx q[0];
rz(-0.22861679) q[0];
x q[1];
rz(-2.1356629) q[2];
sx q[2];
rz(-0.73137368) q[2];
sx q[2];
rz(-1.5767026) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.87344681) q[1];
sx q[1];
rz(-1.6153187) q[1];
sx q[1];
rz(2.4067307) q[1];
rz(-pi) q[2];
rz(1.318526) q[3];
sx q[3];
rz(-3.1384094) q[3];
sx q[3];
rz(-2.0306272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5055351) q[2];
sx q[2];
rz(-1.1670185) q[2];
sx q[2];
rz(2.8324845) q[2];
rz(-1.9348034) q[3];
sx q[3];
rz(-2.759582) q[3];
sx q[3];
rz(-0.25930723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.1768271) q[0];
sx q[0];
rz(-2.4137156) q[0];
sx q[0];
rz(-2.4859909) q[0];
rz(-0.62878311) q[1];
sx q[1];
rz(-1.4097593) q[1];
sx q[1];
rz(1.383673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38046471) q[0];
sx q[0];
rz(-1.5232067) q[0];
sx q[0];
rz(-1.6143198) q[0];
x q[1];
rz(-0.45510095) q[2];
sx q[2];
rz(-2.2544207) q[2];
sx q[2];
rz(0.065082642) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0023345) q[1];
sx q[1];
rz(-2.9911925) q[1];
sx q[1];
rz(0.12087442) q[1];
rz(-pi) q[2];
rz(-1.7114559) q[3];
sx q[3];
rz(-1.0753618) q[3];
sx q[3];
rz(-2.2666033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.76929602) q[2];
sx q[2];
rz(-1.6326222) q[2];
sx q[2];
rz(0.24202913) q[2];
rz(-2.5567143) q[3];
sx q[3];
rz(-2.4681028) q[3];
sx q[3];
rz(-0.52403319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.723421) q[0];
sx q[0];
rz(-2.5192498) q[0];
sx q[0];
rz(2.584516) q[0];
rz(-2.2611025) q[1];
sx q[1];
rz(-1.7387895) q[1];
sx q[1];
rz(-2.1008076) q[1];
rz(-1.3646094) q[2];
sx q[2];
rz(-1.577081) q[2];
sx q[2];
rz(-1.1834363) q[2];
rz(3.0467792) q[3];
sx q[3];
rz(-1.0066114) q[3];
sx q[3];
rz(-0.8361984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
