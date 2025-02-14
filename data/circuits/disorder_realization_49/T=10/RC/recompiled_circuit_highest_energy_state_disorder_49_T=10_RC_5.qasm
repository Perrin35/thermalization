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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5621126) q[0];
sx q[0];
rz(-1.7796374) q[0];
sx q[0];
rz(-3.0203041) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5036132) q[2];
sx q[2];
rz(-2.061741) q[2];
sx q[2];
rz(-0.52396905) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97325402) q[1];
sx q[1];
rz(-1.9562408) q[1];
sx q[1];
rz(-1.5027352) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60987925) q[3];
sx q[3];
rz(-1.0945265) q[3];
sx q[3];
rz(-0.44157883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72630924) q[2];
sx q[2];
rz(-1.3741263) q[2];
sx q[2];
rz(-1.7706002) q[2];
rz(-2.3224984) q[3];
sx q[3];
rz(-2.2947125) q[3];
sx q[3];
rz(3.0313361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6743728) q[0];
sx q[0];
rz(-1.233036) q[0];
sx q[0];
rz(-2.9357173) q[0];
rz(-1.8241833) q[1];
sx q[1];
rz(-1.2818047) q[1];
sx q[1];
rz(2.6726216) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1519794) q[0];
sx q[0];
rz(-1.5474404) q[0];
sx q[0];
rz(-2.4863913) q[0];
rz(-pi) q[1];
rz(-1.1577206) q[2];
sx q[2];
rz(-2.6402355) q[2];
sx q[2];
rz(2.219308) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6862985) q[1];
sx q[1];
rz(-1.0166234) q[1];
sx q[1];
rz(-0.55880736) q[1];
rz(-0.75403611) q[3];
sx q[3];
rz(-1.1234576) q[3];
sx q[3];
rz(-0.099180982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5466902) q[2];
sx q[2];
rz(-2.3089843) q[2];
sx q[2];
rz(-0.62758315) q[2];
rz(-2.0708496) q[3];
sx q[3];
rz(-0.50361931) q[3];
sx q[3];
rz(2.812775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(1.4459389) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86125042) q[0];
sx q[0];
rz(-2.872577) q[0];
sx q[0];
rz(1.6205257) q[0];
rz(-pi) q[1];
rz(-3.0769602) q[2];
sx q[2];
rz(-0.53896133) q[2];
sx q[2];
rz(-2.2825983) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2647977) q[1];
sx q[1];
rz(-1.4709657) q[1];
sx q[1];
rz(-2.8431176) q[1];
rz(-pi) q[2];
rz(0.51579185) q[3];
sx q[3];
rz(-2.502877) q[3];
sx q[3];
rz(2.0206919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0420405) q[2];
sx q[2];
rz(-2.6991762) q[2];
sx q[2];
rz(-0.37910795) q[2];
rz(1.7673309) q[3];
sx q[3];
rz(-2.1702424) q[3];
sx q[3];
rz(-2.1578535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.2079726) q[0];
sx q[0];
rz(-2.3868028) q[0];
sx q[0];
rz(2.2291613) q[0];
rz(-1.747793) q[1];
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
rz(-0.63589225) q[0];
sx q[0];
rz(-1.5109343) q[0];
sx q[0];
rz(-0.12714956) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2851587) q[2];
sx q[2];
rz(-1.2650448) q[2];
sx q[2];
rz(2.6746674) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5279549) q[1];
sx q[1];
rz(-0.21371811) q[1];
sx q[1];
rz(1.3655013) q[1];
rz(-0.87731464) q[3];
sx q[3];
rz(-2.1319509) q[3];
sx q[3];
rz(-1.1059685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.48431524) q[2];
sx q[2];
rz(-2.4399098) q[2];
sx q[2];
rz(0.62937984) q[2];
rz(-1.4194007) q[3];
sx q[3];
rz(-1.281176) q[3];
sx q[3];
rz(-2.4331376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092954271) q[0];
sx q[0];
rz(-2.104367) q[0];
sx q[0];
rz(1.3193489) q[0];
rz(2.0535779) q[1];
sx q[1];
rz(-1.8698147) q[1];
sx q[1];
rz(-1.4647269) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80905238) q[0];
sx q[0];
rz(-2.2702424) q[0];
sx q[0];
rz(0.041985675) q[0];
x q[1];
rz(-3.0973487) q[2];
sx q[2];
rz(-2.508456) q[2];
sx q[2];
rz(2.4447332) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3114918) q[1];
sx q[1];
rz(-1.7602008) q[1];
sx q[1];
rz(-2.7968391) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7978908) q[3];
sx q[3];
rz(-0.97714564) q[3];
sx q[3];
rz(0.81867262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.69270837) q[2];
sx q[2];
rz(-0.11881891) q[2];
sx q[2];
rz(-0.23692712) q[2];
rz(-2.1610625) q[3];
sx q[3];
rz(-1.4292932) q[3];
sx q[3];
rz(-0.94394365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15810814) q[0];
sx q[0];
rz(-3.0321002) q[0];
sx q[0];
rz(-0.14516251) q[0];
rz(1.521184) q[1];
sx q[1];
rz(-0.79997921) q[1];
sx q[1];
rz(3.1343585) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0045429) q[0];
sx q[0];
rz(-2.2570253) q[0];
sx q[0];
rz(-0.30651019) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72609781) q[2];
sx q[2];
rz(-2.0019922) q[2];
sx q[2];
rz(-0.83881718) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0238432) q[1];
sx q[1];
rz(-2.8850318) q[1];
sx q[1];
rz(-1.1505949) q[1];
rz(-2.9754753) q[3];
sx q[3];
rz(-1.5990077) q[3];
sx q[3];
rz(-1.3835554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2468557) q[2];
sx q[2];
rz(-1.9523018) q[2];
sx q[2];
rz(-2.5640633) q[2];
rz(1.9891116) q[3];
sx q[3];
rz(-0.58497506) q[3];
sx q[3];
rz(1.729689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5685527) q[0];
sx q[0];
rz(-2.9459406) q[0];
sx q[0];
rz(-2.0797119) q[0];
rz(1.3878239) q[1];
sx q[1];
rz(-1.9111218) q[1];
sx q[1];
rz(1.6434297) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6681435) q[0];
sx q[0];
rz(-1.7034344) q[0];
sx q[0];
rz(-2.7154865) q[0];
rz(-pi) q[1];
x q[1];
rz(1.241356) q[2];
sx q[2];
rz(-2.1899275) q[2];
sx q[2];
rz(-3.0513024) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.59115852) q[1];
sx q[1];
rz(-0.60748749) q[1];
sx q[1];
rz(1.1689069) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5779547) q[3];
sx q[3];
rz(-2.5803714) q[3];
sx q[3];
rz(-1.3432251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3987223) q[2];
sx q[2];
rz(-0.12464945) q[2];
sx q[2];
rz(-0.91402641) q[2];
rz(-2.3877609) q[3];
sx q[3];
rz(-1.9099312) q[3];
sx q[3];
rz(-0.45647538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0085501) q[0];
sx q[0];
rz(-0.32519105) q[0];
sx q[0];
rz(3.131026) q[0];
rz(-2.1376624) q[1];
sx q[1];
rz(-1.4164475) q[1];
sx q[1];
rz(3.1026057) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4083791) q[0];
sx q[0];
rz(-0.36535242) q[0];
sx q[0];
rz(1.062458) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.580367) q[2];
sx q[2];
rz(-1.8875857) q[2];
sx q[2];
rz(0.94632705) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7379563) q[1];
sx q[1];
rz(-1.9490598) q[1];
sx q[1];
rz(0.72241108) q[1];
rz(-0.66308588) q[3];
sx q[3];
rz(-1.6414101) q[3];
sx q[3];
rz(-0.021573349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5799334) q[2];
sx q[2];
rz(-1.9114405) q[2];
sx q[2];
rz(-0.099543355) q[2];
rz(1.8482515) q[3];
sx q[3];
rz(-1.2852531) q[3];
sx q[3];
rz(-0.37555638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.471591) q[0];
sx q[0];
rz(-1.7741859) q[0];
sx q[0];
rz(-1.8294096) q[0];
rz(0.52681628) q[1];
sx q[1];
rz(-0.75378886) q[1];
sx q[1];
rz(-2.3048293) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1360224) q[0];
sx q[0];
rz(-1.3107398) q[0];
sx q[0];
rz(-1.2664766) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86076059) q[2];
sx q[2];
rz(-2.3902811) q[2];
sx q[2];
rz(-1.4850791) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2075314) q[1];
sx q[1];
rz(-2.4467) q[1];
sx q[1];
rz(2.7259105) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38133867) q[3];
sx q[3];
rz(-1.1313038) q[3];
sx q[3];
rz(0.56909305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3706563) q[2];
sx q[2];
rz(-0.61345658) q[2];
sx q[2];
rz(-1.9343617) q[2];
rz(1.8218482) q[3];
sx q[3];
rz(-1.5765669) q[3];
sx q[3];
rz(-2.3453662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.094548263) q[0];
sx q[0];
rz(-0.47484174) q[0];
sx q[0];
rz(-0.71665254) q[0];
rz(1.5776177) q[1];
sx q[1];
rz(-1.8378704) q[1];
sx q[1];
rz(-2.5242453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43403175) q[0];
sx q[0];
rz(-0.73129485) q[0];
sx q[0];
rz(-1.8282268) q[0];
rz(-pi) q[1];
rz(2.9920299) q[2];
sx q[2];
rz(-2.0345275) q[2];
sx q[2];
rz(0.92926393) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.58891728) q[1];
sx q[1];
rz(-2.008634) q[1];
sx q[1];
rz(-2.4312996) q[1];
rz(-pi) q[2];
rz(-2.8598911) q[3];
sx q[3];
rz(-2.1920077) q[3];
sx q[3];
rz(-3.1007183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3295595) q[2];
sx q[2];
rz(-0.38745189) q[2];
sx q[2];
rz(-0.47920245) q[2];
rz(1.5174348) q[3];
sx q[3];
rz(-1.7370217) q[3];
sx q[3];
rz(1.4994538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
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
rz(-2.3334423) q[1];
sx q[1];
rz(-1.5222526) q[1];
sx q[1];
rz(1.6076988) q[1];
rz(1.343518) q[2];
sx q[2];
rz(-2.6103316) q[2];
sx q[2];
rz(-1.1634315) q[2];
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
