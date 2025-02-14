OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61152148) q[0];
sx q[0];
rz(-1.0728711) q[0];
sx q[0];
rz(1.7142417) q[0];
rz(0.99675769) q[1];
sx q[1];
rz(3.7550959) q[1];
sx q[1];
rz(9.4902314) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9491682) q[0];
sx q[0];
rz(-2.0473366) q[0];
sx q[0];
rz(1.7832827) q[0];
rz(-0.10721389) q[2];
sx q[2];
rz(-2.7709024) q[2];
sx q[2];
rz(3.0840741) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1004384) q[1];
sx q[1];
rz(-1.2693431) q[1];
sx q[1];
rz(-0.90561066) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64163107) q[3];
sx q[3];
rz(-2.4182582) q[3];
sx q[3];
rz(2.1154887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2685711) q[2];
sx q[2];
rz(-0.34653386) q[2];
sx q[2];
rz(-1.2968501) q[2];
rz(-2.0348564) q[3];
sx q[3];
rz(-2.1745671) q[3];
sx q[3];
rz(-1.4510252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035037128) q[0];
sx q[0];
rz(-2.4132001) q[0];
sx q[0];
rz(-1.534071) q[0];
rz(-2.2976177) q[1];
sx q[1];
rz(-1.5183828) q[1];
sx q[1];
rz(-0.27110505) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26626884) q[0];
sx q[0];
rz(-1.3014587) q[0];
sx q[0];
rz(-1.5311509) q[0];
x q[1];
rz(1.9366802) q[2];
sx q[2];
rz(-2.8717763) q[2];
sx q[2];
rz(-1.7797949) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18807377) q[1];
sx q[1];
rz(-1.155335) q[1];
sx q[1];
rz(-1.6759592) q[1];
rz(-0.45174349) q[3];
sx q[3];
rz(-0.16638936) q[3];
sx q[3];
rz(-1.1903669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8035182) q[2];
sx q[2];
rz(-0.98793554) q[2];
sx q[2];
rz(-0.79338497) q[2];
rz(-0.042898305) q[3];
sx q[3];
rz(-1.9148613) q[3];
sx q[3];
rz(0.66810098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0637829) q[0];
sx q[0];
rz(-2.6580647) q[0];
sx q[0];
rz(-0.10312816) q[0];
rz(0.25281301) q[1];
sx q[1];
rz(-1.4272855) q[1];
sx q[1];
rz(-0.3140744) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2677339) q[0];
sx q[0];
rz(-1.1815869) q[0];
sx q[0];
rz(0.53797526) q[0];
rz(-3.0843749) q[2];
sx q[2];
rz(-2.0524745) q[2];
sx q[2];
rz(1.526598) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.16437) q[1];
sx q[1];
rz(-0.80406351) q[1];
sx q[1];
rz(-2.3146446) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7400916) q[3];
sx q[3];
rz(-1.3162656) q[3];
sx q[3];
rz(-2.0930549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0146497) q[2];
sx q[2];
rz(-0.64121556) q[2];
sx q[2];
rz(0.81652299) q[2];
rz(2.4294295) q[3];
sx q[3];
rz(-1.4543337) q[3];
sx q[3];
rz(2.2848391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6461058) q[0];
sx q[0];
rz(-3.0480338) q[0];
sx q[0];
rz(-2.5529472) q[0];
rz(2.7694287) q[1];
sx q[1];
rz(-1.8446422) q[1];
sx q[1];
rz(-2.1968496) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78658453) q[0];
sx q[0];
rz(-1.2517126) q[0];
sx q[0];
rz(3.0727076) q[0];
x q[1];
rz(-2.0633374) q[2];
sx q[2];
rz(-1.8027935) q[2];
sx q[2];
rz(1.7556695) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.041145697) q[1];
sx q[1];
rz(-0.74029912) q[1];
sx q[1];
rz(-2.4802698) q[1];
x q[2];
rz(-2.0382763) q[3];
sx q[3];
rz(-0.40492461) q[3];
sx q[3];
rz(-1.5011464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9113691) q[2];
sx q[2];
rz(-1.5350124) q[2];
sx q[2];
rz(1.4851419) q[2];
rz(-0.39561513) q[3];
sx q[3];
rz(-1.6499358) q[3];
sx q[3];
rz(0.0609456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38451251) q[0];
sx q[0];
rz(-2.7482996) q[0];
sx q[0];
rz(-1.242189) q[0];
rz(3.0643265) q[1];
sx q[1];
rz(-1.2178414) q[1];
sx q[1];
rz(0.088931106) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0384232) q[0];
sx q[0];
rz(-0.22617243) q[0];
sx q[0];
rz(-2.3875176) q[0];
x q[1];
rz(-2.9776666) q[2];
sx q[2];
rz(-0.55183119) q[2];
sx q[2];
rz(2.0275807) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4684852) q[1];
sx q[1];
rz(-1.2947646) q[1];
sx q[1];
rz(2.0620146) q[1];
rz(1.4375646) q[3];
sx q[3];
rz(-1.4733088) q[3];
sx q[3];
rz(-2.9436802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58101216) q[2];
sx q[2];
rz(-1.3835013) q[2];
sx q[2];
rz(-1.2233454) q[2];
rz(1.2286202) q[3];
sx q[3];
rz(-2.2556428) q[3];
sx q[3];
rz(-2.3942053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0864047) q[0];
sx q[0];
rz(-0.815027) q[0];
sx q[0];
rz(-0.95274693) q[0];
rz(1.3505666) q[1];
sx q[1];
rz(-1.2341276) q[1];
sx q[1];
rz(-1.9780673) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02737795) q[0];
sx q[0];
rz(-2.3009536) q[0];
sx q[0];
rz(1.7213836) q[0];
rz(-pi) q[1];
rz(2.7609652) q[2];
sx q[2];
rz(-2.5109595) q[2];
sx q[2];
rz(2.6826114) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1310971) q[1];
sx q[1];
rz(-1.6225885) q[1];
sx q[1];
rz(1.3271922) q[1];
rz(-pi) q[2];
rz(-1.1436449) q[3];
sx q[3];
rz(-2.1329885) q[3];
sx q[3];
rz(-1.027299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20840883) q[2];
sx q[2];
rz(-1.2196502) q[2];
sx q[2];
rz(-1.5677933) q[2];
rz(-2.9257704) q[3];
sx q[3];
rz(-0.62041557) q[3];
sx q[3];
rz(-2.9847434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02483524) q[0];
sx q[0];
rz(-0.76501608) q[0];
sx q[0];
rz(-1.5997973) q[0];
rz(0.4190017) q[1];
sx q[1];
rz(-2.0332789) q[1];
sx q[1];
rz(-1.7760743) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9643819) q[0];
sx q[0];
rz(-2.1809077) q[0];
sx q[0];
rz(2.9047853) q[0];
rz(-pi) q[1];
rz(-2.1183242) q[2];
sx q[2];
rz(-2.27348) q[2];
sx q[2];
rz(1.7349402) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5225181) q[1];
sx q[1];
rz(-0.72942299) q[1];
sx q[1];
rz(-1.3213242) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0715874) q[3];
sx q[3];
rz(-2.4907101) q[3];
sx q[3];
rz(0.18903519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2447723) q[2];
sx q[2];
rz(-0.12910566) q[2];
sx q[2];
rz(1.4531892) q[2];
rz(1.3990654) q[3];
sx q[3];
rz(-1.1996256) q[3];
sx q[3];
rz(-0.80756584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6145265) q[0];
sx q[0];
rz(-0.87166059) q[0];
sx q[0];
rz(0.51323071) q[0];
rz(-0.97663438) q[1];
sx q[1];
rz(-0.66895023) q[1];
sx q[1];
rz(1.4248779) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0289405) q[0];
sx q[0];
rz(-1.7003577) q[0];
sx q[0];
rz(2.2718563) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4972673) q[2];
sx q[2];
rz(-2.1181501) q[2];
sx q[2];
rz(0.076736886) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.91814525) q[1];
sx q[1];
rz(-0.53358101) q[1];
sx q[1];
rz(0.97432889) q[1];
x q[2];
rz(-3.023671) q[3];
sx q[3];
rz(-1.7365082) q[3];
sx q[3];
rz(0.57512257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.49154115) q[2];
sx q[2];
rz(-1.0074793) q[2];
sx q[2];
rz(0.36637351) q[2];
rz(-0.63792396) q[3];
sx q[3];
rz(-1.3831235) q[3];
sx q[3];
rz(1.5095507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76846182) q[0];
sx q[0];
rz(-2.8167384) q[0];
sx q[0];
rz(1.984206) q[0];
rz(-2.6034082) q[1];
sx q[1];
rz(-1.3342074) q[1];
sx q[1];
rz(-0.023712637) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0483646) q[0];
sx q[0];
rz(-0.95010883) q[0];
sx q[0];
rz(-1.6439423) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7269283) q[2];
sx q[2];
rz(-0.86277308) q[2];
sx q[2];
rz(-1.9239349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7249604) q[1];
sx q[1];
rz(-1.9363931) q[1];
sx q[1];
rz(-0.54190062) q[1];
rz(-pi) q[2];
rz(3.0160041) q[3];
sx q[3];
rz(-1.7703345) q[3];
sx q[3];
rz(-0.34717595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9862765) q[2];
sx q[2];
rz(-1.8084869) q[2];
sx q[2];
rz(0.39815608) q[2];
rz(-0.46477535) q[3];
sx q[3];
rz(-2.0538752) q[3];
sx q[3];
rz(1.8006511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5584797) q[0];
sx q[0];
rz(-0.3322424) q[0];
sx q[0];
rz(-0.94859052) q[0];
rz(0.96725431) q[1];
sx q[1];
rz(-1.1525258) q[1];
sx q[1];
rz(-1.0088049) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44701496) q[0];
sx q[0];
rz(-0.80061808) q[0];
sx q[0];
rz(-3.1319764) q[0];
rz(-pi) q[1];
rz(-2.45935) q[2];
sx q[2];
rz(-2.3385404) q[2];
sx q[2];
rz(2.19953) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.12647835) q[1];
sx q[1];
rz(-0.51869828) q[1];
sx q[1];
rz(0.37094231) q[1];
x q[2];
rz(0.28078766) q[3];
sx q[3];
rz(-2.4641345) q[3];
sx q[3];
rz(-1.900713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.41445109) q[2];
sx q[2];
rz(-2.264617) q[2];
sx q[2];
rz(-3.1205175) q[2];
rz(-0.4840788) q[3];
sx q[3];
rz(-1.1295854) q[3];
sx q[3];
rz(-0.67897236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2912343) q[0];
sx q[0];
rz(-1.6073011) q[0];
sx q[0];
rz(-1.4144443) q[0];
rz(-2.5454632) q[1];
sx q[1];
rz(-2.3565751) q[1];
sx q[1];
rz(-0.48859488) q[1];
rz(2.1837744) q[2];
sx q[2];
rz(-0.25479813) q[2];
sx q[2];
rz(0.23164498) q[2];
rz(0.10808839) q[3];
sx q[3];
rz(-1.3664403) q[3];
sx q[3];
rz(-1.4667778) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
