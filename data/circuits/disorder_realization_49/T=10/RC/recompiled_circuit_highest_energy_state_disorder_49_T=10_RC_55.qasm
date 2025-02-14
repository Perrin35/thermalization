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
rz(-2.7558514) q[0];
sx q[0];
rz(-2.1585611) q[0];
sx q[0];
rz(-2.5842343) q[0];
rz(-1.9269257) q[1];
sx q[1];
rz(-2.2872556) q[1];
sx q[1];
rz(0.94205034) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5621126) q[0];
sx q[0];
rz(-1.3619553) q[0];
sx q[0];
rz(3.0203041) q[0];
rz(0.73150191) q[2];
sx q[2];
rz(-0.78353751) q[2];
sx q[2];
rz(-0.48043007) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5696733) q[1];
sx q[1];
rz(-1.5077356) q[1];
sx q[1];
rz(0.38625269) q[1];
rz(-0.73318847) q[3];
sx q[3];
rz(-2.3869143) q[3];
sx q[3];
rz(-1.4316991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4152834) q[2];
sx q[2];
rz(-1.3741263) q[2];
sx q[2];
rz(1.3709925) q[2];
rz(-2.3224984) q[3];
sx q[3];
rz(-0.84688014) q[3];
sx q[3];
rz(-3.0313361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6743728) q[0];
sx q[0];
rz(-1.233036) q[0];
sx q[0];
rz(-0.20587532) q[0];
rz(1.8241833) q[1];
sx q[1];
rz(-1.2818047) q[1];
sx q[1];
rz(0.46897108) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1519794) q[0];
sx q[0];
rz(-1.5474404) q[0];
sx q[0];
rz(-2.4863913) q[0];
rz(-pi) q[1];
rz(-1.1055752) q[2];
sx q[2];
rz(-1.376646) q[2];
sx q[2];
rz(2.126136) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.81509619) q[1];
sx q[1];
rz(-2.3760258) q[1];
sx q[1];
rz(0.86242843) q[1];
rz(-pi) q[2];
rz(-2.5303115) q[3];
sx q[3];
rz(-2.287902) q[3];
sx q[3];
rz(-2.1016013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5466902) q[2];
sx q[2];
rz(-2.3089843) q[2];
sx q[2];
rz(0.62758315) q[2];
rz(1.0707431) q[3];
sx q[3];
rz(-2.6379733) q[3];
sx q[3];
rz(-2.812775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.7367495) q[0];
sx q[0];
rz(-2.3312745) q[0];
sx q[0];
rz(-2.3531083) q[0];
rz(0.57111797) q[1];
sx q[1];
rz(-1.8126789) q[1];
sx q[1];
rz(-1.6956537) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3319237) q[0];
sx q[0];
rz(-1.8394711) q[0];
sx q[0];
rz(-0.013703811) q[0];
rz(1.532191) q[2];
sx q[2];
rz(-1.0330832) q[2];
sx q[2];
rz(0.78372389) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2647977) q[1];
sx q[1];
rz(-1.4709657) q[1];
sx q[1];
rz(0.2984751) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9218677) q[3];
sx q[3];
rz(-1.0255775) q[3];
sx q[3];
rz(-2.6355721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0995522) q[2];
sx q[2];
rz(-0.44241646) q[2];
sx q[2];
rz(0.37910795) q[2];
rz(-1.7673309) q[3];
sx q[3];
rz(-2.1702424) q[3];
sx q[3];
rz(2.1578535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9336201) q[0];
sx q[0];
rz(-0.75478983) q[0];
sx q[0];
rz(0.91243139) q[0];
rz(1.3937996) q[1];
sx q[1];
rz(-0.50519609) q[1];
sx q[1];
rz(2.8392653) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63589225) q[0];
sx q[0];
rz(-1.5109343) q[0];
sx q[0];
rz(-3.0144431) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2851587) q[2];
sx q[2];
rz(-1.8765479) q[2];
sx q[2];
rz(-0.46692525) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.15790882) q[1];
sx q[1];
rz(-1.527546) q[1];
sx q[1];
rz(1.3614348) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3477867) q[3];
sx q[3];
rz(-2.2798139) q[3];
sx q[3];
rz(0.10464917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6572774) q[2];
sx q[2];
rz(-0.70168287) q[2];
sx q[2];
rz(-0.62937984) q[2];
rz(1.4194007) q[3];
sx q[3];
rz(-1.8604167) q[3];
sx q[3];
rz(0.70845503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0486384) q[0];
sx q[0];
rz(-2.104367) q[0];
sx q[0];
rz(-1.8222437) q[0];
rz(-2.0535779) q[1];
sx q[1];
rz(-1.8698147) q[1];
sx q[1];
rz(-1.6768657) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78878337) q[0];
sx q[0];
rz(-1.6029198) q[0];
sx q[0];
rz(2.2706768) q[0];
rz(-pi) q[1];
rz(-1.6032463) q[2];
sx q[2];
rz(-2.2032149) q[2];
sx q[2];
rz(2.4995952) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8147827) q[1];
sx q[1];
rz(-1.2324573) q[1];
sx q[1];
rz(1.7717351) q[1];
rz(-pi) q[2];
rz(-1.7978908) q[3];
sx q[3];
rz(-2.164447) q[3];
sx q[3];
rz(2.32292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69270837) q[2];
sx q[2];
rz(-3.0227737) q[2];
sx q[2];
rz(0.23692712) q[2];
rz(2.1610625) q[3];
sx q[3];
rz(-1.7122995) q[3];
sx q[3];
rz(2.197649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15810814) q[0];
sx q[0];
rz(-3.0321002) q[0];
sx q[0];
rz(2.9964301) q[0];
rz(1.521184) q[1];
sx q[1];
rz(-2.3416134) q[1];
sx q[1];
rz(-3.1343585) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1370498) q[0];
sx q[0];
rz(-0.88456735) q[0];
sx q[0];
rz(2.8350825) q[0];
x q[1];
rz(2.1223567) q[2];
sx q[2];
rz(-0.92331112) q[2];
sx q[2];
rz(1.0873356) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5911094) q[1];
sx q[1];
rz(-1.8046011) q[1];
sx q[1];
rz(0.10661526) q[1];
x q[2];
rz(-2.9725644) q[3];
sx q[3];
rz(-2.9731186) q[3];
sx q[3];
rz(-3.1210312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89473692) q[2];
sx q[2];
rz(-1.1892908) q[2];
sx q[2];
rz(0.57752937) q[2];
rz(1.9891116) q[3];
sx q[3];
rz(-2.5566176) q[3];
sx q[3];
rz(-1.729689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57303992) q[0];
sx q[0];
rz(-0.19565208) q[0];
sx q[0];
rz(-2.0797119) q[0];
rz(1.3878239) q[1];
sx q[1];
rz(-1.2304708) q[1];
sx q[1];
rz(1.498163) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47344917) q[0];
sx q[0];
rz(-1.4381583) q[0];
sx q[0];
rz(0.42610618) q[0];
x q[1];
rz(0.42616828) q[2];
sx q[2];
rz(-2.4505601) q[2];
sx q[2];
rz(2.6994801) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0688144) q[1];
sx q[1];
rz(-2.1238951) q[1];
sx q[1];
rz(-2.8760853) q[1];
rz(-pi) q[2];
rz(2.5779547) q[3];
sx q[3];
rz(-2.5803714) q[3];
sx q[3];
rz(-1.7983675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7428703) q[2];
sx q[2];
rz(-3.0169432) q[2];
sx q[2];
rz(-2.2275662) q[2];
rz(2.3877609) q[3];
sx q[3];
rz(-1.2316615) q[3];
sx q[3];
rz(2.6851173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0085501) q[0];
sx q[0];
rz(-0.32519105) q[0];
sx q[0];
rz(-0.010566674) q[0];
rz(1.0039302) q[1];
sx q[1];
rz(-1.4164475) q[1];
sx q[1];
rz(-0.038987003) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8241725) q[0];
sx q[0];
rz(-1.7455818) q[0];
sx q[0];
rz(1.8932883) q[0];
rz(1.580367) q[2];
sx q[2];
rz(-1.2540069) q[2];
sx q[2];
rz(0.94632705) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5776313) q[1];
sx q[1];
rz(-2.3422514) q[1];
sx q[1];
rz(-0.54117898) q[1];
rz(2.4785068) q[3];
sx q[3];
rz(-1.6414101) q[3];
sx q[3];
rz(3.1200193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5799334) q[2];
sx q[2];
rz(-1.9114405) q[2];
sx q[2];
rz(3.0420493) q[2];
rz(1.8482515) q[3];
sx q[3];
rz(-1.8563396) q[3];
sx q[3];
rz(-2.7660363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67000166) q[0];
sx q[0];
rz(-1.3674068) q[0];
sx q[0];
rz(1.312183) q[0];
rz(-2.6147764) q[1];
sx q[1];
rz(-0.75378886) q[1];
sx q[1];
rz(0.83676338) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6569516) q[0];
sx q[0];
rz(-1.8645727) q[0];
sx q[0];
rz(-2.8696069) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9545045) q[2];
sx q[2];
rz(-1.109668) q[2];
sx q[2];
rz(-2.4949898) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6859582) q[1];
sx q[1];
rz(-2.1966329) q[1];
sx q[1];
rz(1.8954996) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1019245) q[3];
sx q[3];
rz(-1.9143081) q[3];
sx q[3];
rz(1.1706795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3706563) q[2];
sx q[2];
rz(-0.61345658) q[2];
sx q[2];
rz(-1.2072309) q[2];
rz(-1.8218482) q[3];
sx q[3];
rz(-1.5765669) q[3];
sx q[3];
rz(-0.79622644) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0470444) q[0];
sx q[0];
rz(-0.47484174) q[0];
sx q[0];
rz(-0.71665254) q[0];
rz(-1.563975) q[1];
sx q[1];
rz(-1.8378704) q[1];
sx q[1];
rz(-2.5242453) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7075609) q[0];
sx q[0];
rz(-2.4102978) q[0];
sx q[0];
rz(-1.3133658) q[0];
rz(2.9920299) q[2];
sx q[2];
rz(-1.1070651) q[2];
sx q[2];
rz(2.2123287) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58891728) q[1];
sx q[1];
rz(-2.008634) q[1];
sx q[1];
rz(-0.71029305) q[1];
rz(-pi) q[2];
rz(1.2003329) q[3];
sx q[3];
rz(-2.4672567) q[3];
sx q[3];
rz(0.5023027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8120332) q[2];
sx q[2];
rz(-0.38745189) q[2];
sx q[2];
rz(0.47920245) q[2];
rz(-1.6241578) q[3];
sx q[3];
rz(-1.4045709) q[3];
sx q[3];
rz(1.6421389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1351521) q[0];
sx q[0];
rz(-1.3545481) q[0];
sx q[0];
rz(-0.56726278) q[0];
rz(-0.80815036) q[1];
sx q[1];
rz(-1.6193401) q[1];
sx q[1];
rz(-1.5338939) q[1];
rz(-1.343518) q[2];
sx q[2];
rz(-0.53126104) q[2];
sx q[2];
rz(1.9781611) q[2];
rz(1.2841084) q[3];
sx q[3];
rz(-0.78249897) q[3];
sx q[3];
rz(1.8813871) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
