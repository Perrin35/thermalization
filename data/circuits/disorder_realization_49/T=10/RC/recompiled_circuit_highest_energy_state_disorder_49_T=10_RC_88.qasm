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
rz(-0.98303151) q[0];
sx q[0];
rz(-0.55735832) q[0];
rz(-1.9269257) q[1];
sx q[1];
rz(7.1375224) q[1];
sx q[1];
rz(8.4827276) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57948008) q[0];
sx q[0];
rz(-1.3619553) q[0];
sx q[0];
rz(3.0203041) q[0];
rz(-pi) q[1];
rz(2.1579956) q[2];
sx q[2];
rz(-1.0178121) q[2];
sx q[2];
rz(-1.7584973) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1683386) q[1];
sx q[1];
rz(-1.9562408) q[1];
sx q[1];
rz(-1.5027352) q[1];
x q[2];
rz(2.4084042) q[3];
sx q[3];
rz(-0.7546784) q[3];
sx q[3];
rz(1.4316991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4152834) q[2];
sx q[2];
rz(-1.3741263) q[2];
sx q[2];
rz(-1.3709925) q[2];
rz(-0.81909424) q[3];
sx q[3];
rz(-2.2947125) q[3];
sx q[3];
rz(0.11025652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46721989) q[0];
sx q[0];
rz(-1.233036) q[0];
sx q[0];
rz(0.20587532) q[0];
rz(-1.3174093) q[1];
sx q[1];
rz(-1.2818047) q[1];
sx q[1];
rz(-2.6726216) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9896133) q[0];
sx q[0];
rz(-1.5941522) q[0];
sx q[0];
rz(-2.4863913) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9838721) q[2];
sx q[2];
rz(-2.6402355) q[2];
sx q[2];
rz(0.92228466) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9392021) q[1];
sx q[1];
rz(-1.1031045) q[1];
sx q[1];
rz(-2.2012996) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5303115) q[3];
sx q[3];
rz(-0.85369067) q[3];
sx q[3];
rz(-1.0399914) q[3];
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
rz(-2.5140095) q[2];
rz(2.0708496) q[3];
sx q[3];
rz(-0.50361931) q[3];
sx q[3];
rz(-2.812775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7367495) q[0];
sx q[0];
rz(-2.3312745) q[0];
sx q[0];
rz(-0.78848439) q[0];
rz(0.57111797) q[1];
sx q[1];
rz(-1.8126789) q[1];
sx q[1];
rz(1.4459389) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86125042) q[0];
sx q[0];
rz(-0.26901562) q[0];
sx q[0];
rz(1.6205257) q[0];
x q[1];
rz(-0.064632434) q[2];
sx q[2];
rz(-0.53896133) q[2];
sx q[2];
rz(-0.85899437) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72465529) q[1];
sx q[1];
rz(-1.2738527) q[1];
sx q[1];
rz(1.4663803) q[1];
rz(0.57351968) q[3];
sx q[3];
rz(-1.2723426) q[3];
sx q[3];
rz(-0.87707601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0995522) q[2];
sx q[2];
rz(-0.44241646) q[2];
sx q[2];
rz(-2.7624847) q[2];
rz(-1.7673309) q[3];
sx q[3];
rz(-2.1702424) q[3];
sx q[3];
rz(2.1578535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2079726) q[0];
sx q[0];
rz(-2.3868028) q[0];
sx q[0];
rz(2.2291613) q[0];
rz(-1.747793) q[1];
sx q[1];
rz(-2.6363966) q[1];
sx q[1];
rz(-2.8392653) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1990406) q[0];
sx q[0];
rz(-1.4438757) q[0];
sx q[0];
rz(1.5104483) q[0];
x q[1];
rz(2.0197844) q[2];
sx q[2];
rz(-0.76631472) q[2];
sx q[2];
rz(-1.7036071) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7378916) q[1];
sx q[1];
rz(-1.7799592) q[1];
sx q[1];
rz(-3.0973781) q[1];
rz(-pi) q[2];
rz(-2.3477867) q[3];
sx q[3];
rz(-2.2798139) q[3];
sx q[3];
rz(-3.0369435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48431524) q[2];
sx q[2];
rz(-0.70168287) q[2];
sx q[2];
rz(-2.5122128) q[2];
rz(1.4194007) q[3];
sx q[3];
rz(-1.281176) q[3];
sx q[3];
rz(2.4331376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092954271) q[0];
sx q[0];
rz(-1.0372256) q[0];
sx q[0];
rz(1.8222437) q[0];
rz(2.0535779) q[1];
sx q[1];
rz(-1.8698147) q[1];
sx q[1];
rz(1.6768657) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7438904) q[0];
sx q[0];
rz(-0.7004929) q[0];
sx q[0];
rz(-1.6206436) q[0];
rz(-pi) q[1];
rz(-1.6032463) q[2];
sx q[2];
rz(-2.2032149) q[2];
sx q[2];
rz(2.4995952) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32680997) q[1];
sx q[1];
rz(-1.9091354) q[1];
sx q[1];
rz(1.7717351) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60576906) q[3];
sx q[3];
rz(-1.7585227) q[3];
sx q[3];
rz(2.5180205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4488843) q[2];
sx q[2];
rz(-3.0227737) q[2];
sx q[2];
rz(-2.9046655) q[2];
rz(0.98053011) q[3];
sx q[3];
rz(-1.4292932) q[3];
sx q[3];
rz(-0.94394365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9834845) q[0];
sx q[0];
rz(-3.0321002) q[0];
sx q[0];
rz(0.14516251) q[0];
rz(1.6204087) q[1];
sx q[1];
rz(-0.79997921) q[1];
sx q[1];
rz(0.0072341166) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0045429) q[0];
sx q[0];
rz(-0.88456735) q[0];
sx q[0];
rz(-2.8350825) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4154948) q[2];
sx q[2];
rz(-1.1396004) q[2];
sx q[2];
rz(-0.83881718) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.09649) q[1];
sx q[1];
rz(-1.4670925) q[1];
sx q[1];
rz(-1.3357049) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9725644) q[3];
sx q[3];
rz(-2.9731186) q[3];
sx q[3];
rz(-3.1210312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2468557) q[2];
sx q[2];
rz(-1.9523018) q[2];
sx q[2];
rz(2.5640633) q[2];
rz(1.9891116) q[3];
sx q[3];
rz(-2.5566176) q[3];
sx q[3];
rz(1.4119036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57303992) q[0];
sx q[0];
rz(-2.9459406) q[0];
sx q[0];
rz(-2.0797119) q[0];
rz(1.3878239) q[1];
sx q[1];
rz(-1.9111218) q[1];
sx q[1];
rz(-1.498163) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6681435) q[0];
sx q[0];
rz(-1.7034344) q[0];
sx q[0];
rz(0.42610618) q[0];
rz(2.4961124) q[2];
sx q[2];
rz(-1.8374075) q[2];
sx q[2];
rz(-1.4652166) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0688144) q[1];
sx q[1];
rz(-2.1238951) q[1];
sx q[1];
rz(2.8760853) q[1];
rz(-pi) q[2];
rz(1.2467674) q[3];
sx q[3];
rz(-1.1041485) q[3];
sx q[3];
rz(-2.4396536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7428703) q[2];
sx q[2];
rz(-3.0169432) q[2];
sx q[2];
rz(2.2275662) q[2];
rz(2.3877609) q[3];
sx q[3];
rz(-1.9099312) q[3];
sx q[3];
rz(0.45647538) q[3];
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
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1330426) q[0];
sx q[0];
rz(-0.32519105) q[0];
sx q[0];
rz(0.010566674) q[0];
rz(1.0039302) q[1];
sx q[1];
rz(-1.4164475) q[1];
sx q[1];
rz(3.1026057) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31742014) q[0];
sx q[0];
rz(-1.3960109) q[0];
sx q[0];
rz(1.8932883) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5612256) q[2];
sx q[2];
rz(-1.2540069) q[2];
sx q[2];
rz(-2.1952656) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5639613) q[1];
sx q[1];
rz(-0.79934123) q[1];
sx q[1];
rz(-0.54117898) q[1];
rz(-1.6603062) q[3];
sx q[3];
rz(-0.90965547) q[3];
sx q[3];
rz(-1.6474219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5799334) q[2];
sx q[2];
rz(-1.9114405) q[2];
sx q[2];
rz(-0.099543355) q[2];
rz(1.8482515) q[3];
sx q[3];
rz(-1.8563396) q[3];
sx q[3];
rz(0.37555638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.471591) q[0];
sx q[0];
rz(-1.3674068) q[0];
sx q[0];
rz(-1.8294096) q[0];
rz(0.52681628) q[1];
sx q[1];
rz(-2.3878038) q[1];
sx q[1];
rz(2.3048293) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1360224) q[0];
sx q[0];
rz(-1.3107398) q[0];
sx q[0];
rz(-1.8751161) q[0];
rz(2.1870881) q[2];
sx q[2];
rz(-1.109668) q[2];
sx q[2];
rz(2.4949898) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4556344) q[1];
sx q[1];
rz(-0.94495979) q[1];
sx q[1];
rz(-1.246093) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2403735) q[3];
sx q[3];
rz(-0.57358426) q[3];
sx q[3];
rz(0.18665126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7709363) q[2];
sx q[2];
rz(-0.61345658) q[2];
sx q[2];
rz(1.2072309) q[2];
rz(1.3197445) q[3];
sx q[3];
rz(-1.5650257) q[3];
sx q[3];
rz(0.79622644) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094548263) q[0];
sx q[0];
rz(-2.6667509) q[0];
sx q[0];
rz(2.4249401) q[0];
rz(1.563975) q[1];
sx q[1];
rz(-1.8378704) q[1];
sx q[1];
rz(2.5242453) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43403175) q[0];
sx q[0];
rz(-2.4102978) q[0];
sx q[0];
rz(1.3133658) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2812219) q[2];
sx q[2];
rz(-2.6560142) q[2];
sx q[2];
rz(-2.537279) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6177298) q[1];
sx q[1];
rz(-0.81392127) q[1];
sx q[1];
rz(2.518923) q[1];
rz(-pi) q[2];
rz(0.93043296) q[3];
sx q[3];
rz(-1.3427707) q[3];
sx q[3];
rz(1.3630529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8120332) q[2];
sx q[2];
rz(-2.7541408) q[2];
sx q[2];
rz(-0.47920245) q[2];
rz(-1.5174348) q[3];
sx q[3];
rz(-1.7370217) q[3];
sx q[3];
rz(1.6421389) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1351521) q[0];
sx q[0];
rz(-1.3545481) q[0];
sx q[0];
rz(-0.56726278) q[0];
rz(-2.3334423) q[1];
sx q[1];
rz(-1.5222526) q[1];
sx q[1];
rz(1.6076988) q[1];
rz(2.0907503) q[2];
sx q[2];
rz(-1.4563917) q[2];
sx q[2];
rz(0.2105486) q[2];
rz(-1.8574842) q[3];
sx q[3];
rz(-0.78249897) q[3];
sx q[3];
rz(1.8813871) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
