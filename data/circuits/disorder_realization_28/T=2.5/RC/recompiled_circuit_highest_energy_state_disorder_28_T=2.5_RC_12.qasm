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
rz(0.75594354) q[0];
sx q[0];
rz(-2.8109901) q[0];
sx q[0];
rz(-0.27275738) q[0];
rz(2.7859712) q[1];
sx q[1];
rz(-2.1115117) q[1];
sx q[1];
rz(1.5707836) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8518675) q[0];
sx q[0];
rz(-3.0209732) q[0];
sx q[0];
rz(0.9319181) q[0];
rz(-0.46478125) q[2];
sx q[2];
rz(-1.2743093) q[2];
sx q[2];
rz(-3.1303984) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.570191) q[1];
sx q[1];
rz(-2.4610011) q[1];
sx q[1];
rz(-3.0742744) q[1];
rz(-0.36873534) q[3];
sx q[3];
rz(-1.4793385) q[3];
sx q[3];
rz(1.2386595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.97751272) q[2];
sx q[2];
rz(-1.7089607) q[2];
sx q[2];
rz(-0.88708964) q[2];
rz(1.2648434) q[3];
sx q[3];
rz(-2.0829468) q[3];
sx q[3];
rz(3.1177055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7510659) q[0];
sx q[0];
rz(-2.2845415) q[0];
sx q[0];
rz(-0.29399011) q[0];
rz(-1.3514306) q[1];
sx q[1];
rz(-2.0452812) q[1];
sx q[1];
rz(-1.9757804) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6751926) q[0];
sx q[0];
rz(-1.8972016) q[0];
sx q[0];
rz(2.9221888) q[0];
rz(-pi) q[1];
rz(1.4977156) q[2];
sx q[2];
rz(-1.5982619) q[2];
sx q[2];
rz(-2.0171201) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4519074) q[1];
sx q[1];
rz(-2.4396585) q[1];
sx q[1];
rz(-2.4268389) q[1];
rz(-2.5727859) q[3];
sx q[3];
rz(-1.3242464) q[3];
sx q[3];
rz(-0.55014551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4531648) q[2];
sx q[2];
rz(-2.5025949) q[2];
sx q[2];
rz(1.2028018) q[2];
rz(1.5419675) q[3];
sx q[3];
rz(-1.8034233) q[3];
sx q[3];
rz(0.28551027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0315392) q[0];
sx q[0];
rz(-3.0634614) q[0];
sx q[0];
rz(1.9670638) q[0];
rz(2.166676) q[1];
sx q[1];
rz(-0.87923032) q[1];
sx q[1];
rz(-0.45641315) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.31118) q[0];
sx q[0];
rz(-1.5390656) q[0];
sx q[0];
rz(-1.8422442) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1295616) q[2];
sx q[2];
rz(-0.44545275) q[2];
sx q[2];
rz(2.403806) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.76589903) q[1];
sx q[1];
rz(-1.0119034) q[1];
sx q[1];
rz(2.2431348) q[1];
rz(-pi) q[2];
rz(0.47243709) q[3];
sx q[3];
rz(-0.91210134) q[3];
sx q[3];
rz(2.6078826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5259033) q[2];
sx q[2];
rz(-0.38280767) q[2];
sx q[2];
rz(-2.3846386) q[2];
rz(-3.1225539) q[3];
sx q[3];
rz(-1.2913387) q[3];
sx q[3];
rz(2.3058057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86209908) q[0];
sx q[0];
rz(-2.4308496) q[0];
sx q[0];
rz(-0.94974649) q[0];
rz(0.21385916) q[1];
sx q[1];
rz(-2.2016134) q[1];
sx q[1];
rz(2.5423539) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99357641) q[0];
sx q[0];
rz(-0.92088078) q[0];
sx q[0];
rz(2.9666569) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2798115) q[2];
sx q[2];
rz(-2.6070234) q[2];
sx q[2];
rz(0.9415516) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.039700198) q[1];
sx q[1];
rz(-1.779161) q[1];
sx q[1];
rz(2.1096257) q[1];
rz(-pi) q[2];
rz(-1.4324801) q[3];
sx q[3];
rz(-0.99405046) q[3];
sx q[3];
rz(-2.6111297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.63138258) q[2];
sx q[2];
rz(-1.6333132) q[2];
sx q[2];
rz(-0.82376662) q[2];
rz(2.0508164) q[3];
sx q[3];
rz(-2.3407276) q[3];
sx q[3];
rz(-0.66528475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3983696) q[0];
sx q[0];
rz(-0.38341612) q[0];
sx q[0];
rz(-0.98633352) q[0];
rz(2.8979454) q[1];
sx q[1];
rz(-1.8758834) q[1];
sx q[1];
rz(-0.15959127) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8411113) q[0];
sx q[0];
rz(-0.87219014) q[0];
sx q[0];
rz(-0.41450702) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56964376) q[2];
sx q[2];
rz(-2.2692371) q[2];
sx q[2];
rz(1.4739571) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7437834) q[1];
sx q[1];
rz(-1.127217) q[1];
sx q[1];
rz(3.0438782) q[1];
x q[2];
rz(-2.2714204) q[3];
sx q[3];
rz(-1.3425832) q[3];
sx q[3];
rz(-2.1765577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7470982) q[2];
sx q[2];
rz(-1.8793841) q[2];
sx q[2];
rz(-2.0060284) q[2];
rz(-0.46843946) q[3];
sx q[3];
rz(-1.2218385) q[3];
sx q[3];
rz(2.2170317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0202476) q[0];
sx q[0];
rz(-1.0806885) q[0];
sx q[0];
rz(2.8998846) q[0];
rz(1.7860335) q[1];
sx q[1];
rz(-0.87239289) q[1];
sx q[1];
rz(2.2579069) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6304105) q[0];
sx q[0];
rz(-0.97994655) q[0];
sx q[0];
rz(-2.1686337) q[0];
x q[1];
rz(-0.31225754) q[2];
sx q[2];
rz(-1.3259282) q[2];
sx q[2];
rz(-1.6409525) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8203796) q[1];
sx q[1];
rz(-1.1324778) q[1];
sx q[1];
rz(-0.21578101) q[1];
rz(-1.571322) q[3];
sx q[3];
rz(-2.7306203) q[3];
sx q[3];
rz(-2.4961159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22150618) q[2];
sx q[2];
rz(-0.93044996) q[2];
sx q[2];
rz(-0.42631701) q[2];
rz(0.91655556) q[3];
sx q[3];
rz(-1.5519658) q[3];
sx q[3];
rz(-2.0385273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.351848) q[0];
sx q[0];
rz(-1.9639356) q[0];
sx q[0];
rz(-2.015693) q[0];
rz(-0.060983505) q[1];
sx q[1];
rz(-0.95044249) q[1];
sx q[1];
rz(-2.1263863) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39097745) q[0];
sx q[0];
rz(-2.1732507) q[0];
sx q[0];
rz(-1.0573122) q[0];
rz(-0.073748579) q[2];
sx q[2];
rz(-1.2332067) q[2];
sx q[2];
rz(0.13843564) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0132844) q[1];
sx q[1];
rz(-0.40402141) q[1];
sx q[1];
rz(1.4133417) q[1];
x q[2];
rz(-1.6105846) q[3];
sx q[3];
rz(-1.4574175) q[3];
sx q[3];
rz(1.6802406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.31275493) q[2];
sx q[2];
rz(-1.2042896) q[2];
sx q[2];
rz(-0.22605669) q[2];
rz(0.95212805) q[3];
sx q[3];
rz(-3.003037) q[3];
sx q[3];
rz(2.3329195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011618135) q[0];
sx q[0];
rz(-1.8267153) q[0];
sx q[0];
rz(-0.34039482) q[0];
rz(-3.0420692) q[1];
sx q[1];
rz(-2.0390022) q[1];
sx q[1];
rz(1.5460825) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6130644) q[0];
sx q[0];
rz(-0.97672594) q[0];
sx q[0];
rz(-1.1737002) q[0];
rz(0.056273862) q[2];
sx q[2];
rz(-1.2985126) q[2];
sx q[2];
rz(0.19993609) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4027962) q[1];
sx q[1];
rz(-1.5288884) q[1];
sx q[1];
rz(-0.099385029) q[1];
rz(-pi) q[2];
rz(0.45507064) q[3];
sx q[3];
rz(-1.651147) q[3];
sx q[3];
rz(-0.72380607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.16326363) q[2];
sx q[2];
rz(-1.1922057) q[2];
sx q[2];
rz(-2.1400129) q[2];
rz(-1.7892276) q[3];
sx q[3];
rz(-0.47836256) q[3];
sx q[3];
rz(-1.1866239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1773862) q[0];
sx q[0];
rz(-1.0564251) q[0];
sx q[0];
rz(-0.00038432234) q[0];
rz(0.43801019) q[1];
sx q[1];
rz(-1.6338467) q[1];
sx q[1];
rz(2.6843574) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23822902) q[0];
sx q[0];
rz(-1.6203124) q[0];
sx q[0];
rz(1.7286506) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3136131) q[2];
sx q[2];
rz(-2.3739907) q[2];
sx q[2];
rz(2.0274251) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.25343382) q[1];
sx q[1];
rz(-2.7256084) q[1];
sx q[1];
rz(0.57749727) q[1];
rz(-2.5874279) q[3];
sx q[3];
rz(-2.5636052) q[3];
sx q[3];
rz(-2.223034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3823201) q[2];
sx q[2];
rz(-0.98733941) q[2];
sx q[2];
rz(-1.5116723) q[2];
rz(0.090653732) q[3];
sx q[3];
rz(-1.0414618) q[3];
sx q[3];
rz(2.4618728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35336211) q[0];
sx q[0];
rz(-1.7857977) q[0];
sx q[0];
rz(-0.28025383) q[0];
rz(1.7578112) q[1];
sx q[1];
rz(-0.46418142) q[1];
sx q[1];
rz(1.9162477) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2458151) q[0];
sx q[0];
rz(-1.4207463) q[0];
sx q[0];
rz(1.4729985) q[0];
rz(-pi) q[1];
rz(-2.7011048) q[2];
sx q[2];
rz(-2.2058528) q[2];
sx q[2];
rz(1.6556513) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6514488) q[1];
sx q[1];
rz(-1.3281392) q[1];
sx q[1];
rz(-0.76702422) q[1];
x q[2];
rz(-0.7627714) q[3];
sx q[3];
rz(-2.7342334) q[3];
sx q[3];
rz(-0.4688022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.19259351) q[2];
sx q[2];
rz(-1.5263564) q[2];
sx q[2];
rz(0.68142778) q[2];
rz(2.0415908) q[3];
sx q[3];
rz(-0.93183485) q[3];
sx q[3];
rz(-1.3070235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6901907) q[0];
sx q[0];
rz(-0.9342397) q[0];
sx q[0];
rz(-1.1697212) q[0];
rz(-0.38047094) q[1];
sx q[1];
rz(-0.66743757) q[1];
sx q[1];
rz(0.22421511) q[1];
rz(-1.8809594) q[2];
sx q[2];
rz(-0.72523965) q[2];
sx q[2];
rz(0.88015866) q[2];
rz(-2.7093665) q[3];
sx q[3];
rz(-2.1189011) q[3];
sx q[3];
rz(-0.44174474) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
