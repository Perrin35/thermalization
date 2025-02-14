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
rz(0.3857412) q[0];
sx q[0];
rz(2.1585611) q[0];
sx q[0];
rz(9.9821363) q[0];
rz(-1.9269257) q[1];
sx q[1];
rz(-2.2872556) q[1];
sx q[1];
rz(-2.1995423) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0165812) q[0];
sx q[0];
rz(-1.4521557) q[0];
sx q[0];
rz(-1.3604547) q[0];
x q[1];
rz(0.73150191) q[2];
sx q[2];
rz(-2.3580551) q[2];
sx q[2];
rz(-2.6611626) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1526133) q[1];
sx q[1];
rz(-2.7504814) q[1];
sx q[1];
rz(2.9755201) q[1];
rz(-pi) q[2];
rz(-0.60987925) q[3];
sx q[3];
rz(-1.0945265) q[3];
sx q[3];
rz(0.44157883) q[3];
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
rz(2.3224984) q[3];
sx q[3];
rz(-0.84688014) q[3];
sx q[3];
rz(-0.11025652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46721989) q[0];
sx q[0];
rz(-1.233036) q[0];
sx q[0];
rz(2.9357173) q[0];
rz(1.8241833) q[1];
sx q[1];
rz(-1.2818047) q[1];
sx q[1];
rz(0.46897108) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9896133) q[0];
sx q[0];
rz(-1.5474404) q[0];
sx q[0];
rz(0.65520136) q[0];
rz(-pi) q[1];
rz(-2.9250336) q[2];
sx q[2];
rz(-1.1149841) q[2];
sx q[2];
rz(2.6827981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3264965) q[1];
sx q[1];
rz(-0.76556682) q[1];
sx q[1];
rz(-2.2791642) q[1];
rz(-pi) q[2];
rz(-2.5303115) q[3];
sx q[3];
rz(-0.85369067) q[3];
sx q[3];
rz(-1.0399914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59490243) q[2];
sx q[2];
rz(-2.3089843) q[2];
sx q[2];
rz(-2.5140095) q[2];
rz(2.0708496) q[3];
sx q[3];
rz(-2.6379733) q[3];
sx q[3];
rz(-0.32881769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7367495) q[0];
sx q[0];
rz(-2.3312745) q[0];
sx q[0];
rz(2.3531083) q[0];
rz(-0.57111797) q[1];
sx q[1];
rz(-1.3289137) q[1];
sx q[1];
rz(-1.6956537) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3841032) q[0];
sx q[0];
rz(-1.5575842) q[0];
sx q[0];
rz(1.8394952) q[0];
rz(-pi) q[1];
rz(-0.064632434) q[2];
sx q[2];
rz(-2.6026313) q[2];
sx q[2];
rz(0.85899437) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4169374) q[1];
sx q[1];
rz(-1.8677399) q[1];
sx q[1];
rz(1.6752123) q[1];
rz(1.219725) q[3];
sx q[3];
rz(-2.1160151) q[3];
sx q[3];
rz(2.6355721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0995522) q[2];
sx q[2];
rz(-0.44241646) q[2];
sx q[2];
rz(-2.7624847) q[2];
rz(-1.3742617) q[3];
sx q[3];
rz(-2.1702424) q[3];
sx q[3];
rz(-2.1578535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2079726) q[0];
sx q[0];
rz(-0.75478983) q[0];
sx q[0];
rz(2.2291613) q[0];
rz(1.747793) q[1];
sx q[1];
rz(-2.6363966) q[1];
sx q[1];
rz(2.8392653) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5057004) q[0];
sx q[0];
rz(-1.5109343) q[0];
sx q[0];
rz(3.0144431) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85643391) q[2];
sx q[2];
rz(-1.2650448) q[2];
sx q[2];
rz(-2.6746674) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7378916) q[1];
sx q[1];
rz(-1.3616335) q[1];
sx q[1];
rz(3.0973781) q[1];
x q[2];
rz(-2.3477867) q[3];
sx q[3];
rz(-2.2798139) q[3];
sx q[3];
rz(-3.0369435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6572774) q[2];
sx q[2];
rz(-0.70168287) q[2];
sx q[2];
rz(-2.5122128) q[2];
rz(-1.7221919) q[3];
sx q[3];
rz(-1.281176) q[3];
sx q[3];
rz(2.4331376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0486384) q[0];
sx q[0];
rz(-1.0372256) q[0];
sx q[0];
rz(1.8222437) q[0];
rz(-1.0880148) q[1];
sx q[1];
rz(-1.8698147) q[1];
sx q[1];
rz(-1.4647269) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80905238) q[0];
sx q[0];
rz(-2.2702424) q[0];
sx q[0];
rz(-0.041985675) q[0];
rz(0.63266968) q[2];
sx q[2];
rz(-1.5969689) q[2];
sx q[2];
rz(-2.2319792) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22359554) q[1];
sx q[1];
rz(-0.3915266) q[1];
sx q[1];
rz(-2.6256204) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3220114) q[3];
sx q[3];
rz(-2.5109041) q[3];
sx q[3];
rz(1.2104152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4488843) q[2];
sx q[2];
rz(-0.11881891) q[2];
sx q[2];
rz(-0.23692712) q[2];
rz(-2.1610625) q[3];
sx q[3];
rz(-1.7122995) q[3];
sx q[3];
rz(-2.197649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9834845) q[0];
sx q[0];
rz(-3.0321002) q[0];
sx q[0];
rz(0.14516251) q[0];
rz(-1.6204087) q[1];
sx q[1];
rz(-2.3416134) q[1];
sx q[1];
rz(0.0072341166) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4677759) q[0];
sx q[0];
rz(-2.4002808) q[0];
sx q[0];
rz(-1.2178161) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0192359) q[2];
sx q[2];
rz(-2.2182815) q[2];
sx q[2];
rz(-1.0873356) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.09649) q[1];
sx q[1];
rz(-1.4670925) q[1];
sx q[1];
rz(-1.8058877) q[1];
rz(1.5421914) q[3];
sx q[3];
rz(-1.4047457) q[3];
sx q[3];
rz(-0.19197026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.89473692) q[2];
sx q[2];
rz(-1.9523018) q[2];
sx q[2];
rz(-2.5640633) q[2];
rz(1.1524811) q[3];
sx q[3];
rz(-2.5566176) q[3];
sx q[3];
rz(1.729689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
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
rz(-1.3878239) q[1];
sx q[1];
rz(-1.2304708) q[1];
sx q[1];
rz(1.6434297) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9842872) q[0];
sx q[0];
rz(-1.1486736) q[0];
sx q[0];
rz(-1.4253084) q[0];
rz(1.241356) q[2];
sx q[2];
rz(-2.1899275) q[2];
sx q[2];
rz(-3.0513024) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5504341) q[1];
sx q[1];
rz(-2.5341052) q[1];
sx q[1];
rz(1.1689069) q[1];
rz(2.6531336) q[3];
sx q[3];
rz(-1.2824713) q[3];
sx q[3];
rz(0.718887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3987223) q[2];
sx q[2];
rz(-3.0169432) q[2];
sx q[2];
rz(-0.91402641) q[2];
rz(0.75383178) q[3];
sx q[3];
rz(-1.9099312) q[3];
sx q[3];
rz(-0.45647538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.0085501) q[0];
sx q[0];
rz(-0.32519105) q[0];
sx q[0];
rz(3.131026) q[0];
rz(-1.0039302) q[1];
sx q[1];
rz(-1.4164475) q[1];
sx q[1];
rz(-3.1026057) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8241725) q[0];
sx q[0];
rz(-1.7455818) q[0];
sx q[0];
rz(-1.8932883) q[0];
rz(-pi) q[1];
rz(2.8247897) q[2];
sx q[2];
rz(-1.5798908) q[2];
sx q[2];
rz(0.62745079) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7379563) q[1];
sx q[1];
rz(-1.9490598) q[1];
sx q[1];
rz(2.4191816) q[1];
rz(-pi) q[2];
rz(-1.6603062) q[3];
sx q[3];
rz(-2.2319372) q[3];
sx q[3];
rz(-1.4941708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5799334) q[2];
sx q[2];
rz(-1.2301521) q[2];
sx q[2];
rz(0.099543355) q[2];
rz(1.2933412) q[3];
sx q[3];
rz(-1.8563396) q[3];
sx q[3];
rz(-0.37555638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67000166) q[0];
sx q[0];
rz(-1.7741859) q[0];
sx q[0];
rz(1.8294096) q[0];
rz(-0.52681628) q[1];
sx q[1];
rz(-0.75378886) q[1];
sx q[1];
rz(-0.83676338) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89032407) q[0];
sx q[0];
rz(-2.7439372) q[0];
sx q[0];
rz(-2.2969382) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2808321) q[2];
sx q[2];
rz(-2.3902811) q[2];
sx q[2];
rz(-1.6565135) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6859582) q[1];
sx q[1];
rz(-2.1966329) q[1];
sx q[1];
rz(1.8954996) q[1];
rz(-pi) q[2];
rz(-1.1019245) q[3];
sx q[3];
rz(-1.2272845) q[3];
sx q[3];
rz(1.1706795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3706563) q[2];
sx q[2];
rz(-0.61345658) q[2];
sx q[2];
rz(-1.2072309) q[2];
rz(1.3197445) q[3];
sx q[3];
rz(-1.5765669) q[3];
sx q[3];
rz(2.3453662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0470444) q[0];
sx q[0];
rz(-0.47484174) q[0];
sx q[0];
rz(-2.4249401) q[0];
rz(-1.563975) q[1];
sx q[1];
rz(-1.3037222) q[1];
sx q[1];
rz(-0.61734739) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43403175) q[0];
sx q[0];
rz(-0.73129485) q[0];
sx q[0];
rz(1.3133658) q[0];
rz(0.1495628) q[2];
sx q[2];
rz(-1.1070651) q[2];
sx q[2];
rz(0.92926393) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3315369) q[1];
sx q[1];
rz(-0.9390587) q[1];
sx q[1];
rz(1.017636) q[1];
rz(-pi) q[2];
rz(-0.93043296) q[3];
sx q[3];
rz(-1.798822) q[3];
sx q[3];
rz(-1.7785398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8120332) q[2];
sx q[2];
rz(-2.7541408) q[2];
sx q[2];
rz(2.6623902) q[2];
rz(-1.5174348) q[3];
sx q[3];
rz(-1.7370217) q[3];
sx q[3];
rz(-1.4994538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.0064405) q[0];
sx q[0];
rz(-1.3545481) q[0];
sx q[0];
rz(-0.56726278) q[0];
rz(-0.80815036) q[1];
sx q[1];
rz(-1.6193401) q[1];
sx q[1];
rz(-1.5338939) q[1];
rz(-1.7980747) q[2];
sx q[2];
rz(-2.6103316) q[2];
sx q[2];
rz(-1.1634315) q[2];
rz(-2.8675254) q[3];
sx q[3];
rz(-0.82809877) q[3];
sx q[3];
rz(-0.86622151) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
