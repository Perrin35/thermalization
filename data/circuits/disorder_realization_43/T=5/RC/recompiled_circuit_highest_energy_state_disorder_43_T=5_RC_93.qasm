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
rz(-0.22696683) q[0];
sx q[0];
rz(-2.5300955) q[0];
sx q[0];
rz(-0.55846941) q[0];
rz(0.59977579) q[1];
sx q[1];
rz(1.37473) q[1];
sx q[1];
rz(9.3461499) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087705165) q[0];
sx q[0];
rz(-1.0358521) q[0];
sx q[0];
rz(-2.22552) q[0];
rz(-pi) q[1];
rz(2.9869674) q[2];
sx q[2];
rz(-0.68462709) q[2];
sx q[2];
rz(-1.3473617) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.47657644) q[1];
sx q[1];
rz(-1.652583) q[1];
sx q[1];
rz(0.25340124) q[1];
x q[2];
rz(2.4047002) q[3];
sx q[3];
rz(-0.46023638) q[3];
sx q[3];
rz(2.5207507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.89995304) q[2];
sx q[2];
rz(-3.0934379) q[2];
sx q[2];
rz(1.0345577) q[2];
rz(-3.0311846) q[3];
sx q[3];
rz(-1.4805099) q[3];
sx q[3];
rz(0.30606562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13650525) q[0];
sx q[0];
rz(-1.6965447) q[0];
sx q[0];
rz(1.1862296) q[0];
rz(1.8738481) q[1];
sx q[1];
rz(-0.45919752) q[1];
sx q[1];
rz(1.7523821) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25042501) q[0];
sx q[0];
rz(-2.1019263) q[0];
sx q[0];
rz(2.7851339) q[0];
x q[1];
rz(0.48149646) q[2];
sx q[2];
rz(-1.4717558) q[2];
sx q[2];
rz(1.5426829) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4965044) q[1];
sx q[1];
rz(-1.2821272) q[1];
sx q[1];
rz(1.9157857) q[1];
rz(-pi) q[2];
rz(1.4382382) q[3];
sx q[3];
rz(-2.5131559) q[3];
sx q[3];
rz(-0.6836764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.17239751) q[2];
sx q[2];
rz(-1.1149422) q[2];
sx q[2];
rz(-1.4364852) q[2];
rz(-1.8819594) q[3];
sx q[3];
rz(-0.45537046) q[3];
sx q[3];
rz(0.026738515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36419511) q[0];
sx q[0];
rz(-1.3854249) q[0];
sx q[0];
rz(-1.2170323) q[0];
rz(-2.9265535) q[1];
sx q[1];
rz(-1.2478849) q[1];
sx q[1];
rz(1.7020114) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0501409) q[0];
sx q[0];
rz(-0.9147075) q[0];
sx q[0];
rz(-1.4138282) q[0];
rz(-pi) q[1];
rz(0.093482253) q[2];
sx q[2];
rz(-2.0801525) q[2];
sx q[2];
rz(-2.2619154) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7172723) q[1];
sx q[1];
rz(-1.3417146) q[1];
sx q[1];
rz(-3.0113264) q[1];
rz(0.21814367) q[3];
sx q[3];
rz(-1.3036332) q[3];
sx q[3];
rz(1.657079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3415459) q[2];
sx q[2];
rz(-1.505082) q[2];
sx q[2];
rz(-0.81614196) q[2];
rz(0.0084361313) q[3];
sx q[3];
rz(-0.71788994) q[3];
sx q[3];
rz(-1.5661904) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7208045) q[0];
sx q[0];
rz(-1.1710465) q[0];
sx q[0];
rz(2.8724443) q[0];
rz(-1.7587657) q[1];
sx q[1];
rz(-2.5803284) q[1];
sx q[1];
rz(-0.47163481) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1203493) q[0];
sx q[0];
rz(-2.2902852) q[0];
sx q[0];
rz(0.18622744) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.083375562) q[2];
sx q[2];
rz(-2.7365587) q[2];
sx q[2];
rz(2.0138559) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.562512) q[1];
sx q[1];
rz(-1.0483841) q[1];
sx q[1];
rz(1.9396616) q[1];
rz(-2.7876623) q[3];
sx q[3];
rz(-0.9359064) q[3];
sx q[3];
rz(0.77087444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.65583324) q[2];
sx q[2];
rz(-1.4071608) q[2];
sx q[2];
rz(-1.887623) q[2];
rz(-0.29821011) q[3];
sx q[3];
rz(-1.8658274) q[3];
sx q[3];
rz(1.5378753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0373847) q[0];
sx q[0];
rz(-0.91008121) q[0];
sx q[0];
rz(0.75697672) q[0];
rz(2.0825999) q[1];
sx q[1];
rz(-1.4451566) q[1];
sx q[1];
rz(2.7991378) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.379102) q[0];
sx q[0];
rz(-0.44371334) q[0];
sx q[0];
rz(1.1136454) q[0];
rz(-pi) q[1];
x q[1];
rz(1.191869) q[2];
sx q[2];
rz(-2.4034326) q[2];
sx q[2];
rz(-1.2305232) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40552822) q[1];
sx q[1];
rz(-1.9384512) q[1];
sx q[1];
rz(2.6339732) q[1];
x q[2];
rz(1.2737688) q[3];
sx q[3];
rz(-1.1714463) q[3];
sx q[3];
rz(0.12055971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4464438) q[2];
sx q[2];
rz(-1.6764287) q[2];
sx q[2];
rz(-2.5214419) q[2];
rz(-0.54689637) q[3];
sx q[3];
rz(-0.20709012) q[3];
sx q[3];
rz(1.0774405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3243489) q[0];
sx q[0];
rz(-2.9381848) q[0];
sx q[0];
rz(1.2212344) q[0];
rz(1.1034032) q[1];
sx q[1];
rz(-1.1627448) q[1];
sx q[1];
rz(2.3308636) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2632836) q[0];
sx q[0];
rz(-1.7297525) q[0];
sx q[0];
rz(-2.8829888) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76580183) q[2];
sx q[2];
rz(-2.3277936) q[2];
sx q[2];
rz(1.9794996) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83119694) q[1];
sx q[1];
rz(-1.4321021) q[1];
sx q[1];
rz(3.0866361) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91355904) q[3];
sx q[3];
rz(-2.676042) q[3];
sx q[3];
rz(-0.54009932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3409884) q[2];
sx q[2];
rz(-1.1977414) q[2];
sx q[2];
rz(1.1654589) q[2];
rz(0.10945877) q[3];
sx q[3];
rz(-2.1200924) q[3];
sx q[3];
rz(-0.060733184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.3133746) q[0];
sx q[0];
rz(-3.1309541) q[0];
sx q[0];
rz(-2.4995372) q[0];
rz(2.8984046) q[1];
sx q[1];
rz(-2.7311192) q[1];
sx q[1];
rz(2.5004255) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12578861) q[0];
sx q[0];
rz(-2.1282853) q[0];
sx q[0];
rz(-0.95161887) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2095318) q[2];
sx q[2];
rz(-0.85339499) q[2];
sx q[2];
rz(-0.016005767) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23785628) q[1];
sx q[1];
rz(-2.3163539) q[1];
sx q[1];
rz(-1.115487) q[1];
rz(-pi) q[2];
rz(1.0585467) q[3];
sx q[3];
rz(-0.44971684) q[3];
sx q[3];
rz(2.0459156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13333653) q[2];
sx q[2];
rz(-1.9828321) q[2];
sx q[2];
rz(-2.1978417) q[2];
rz(-2.4933695) q[3];
sx q[3];
rz(-0.74130487) q[3];
sx q[3];
rz(-0.6664204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7159395) q[0];
sx q[0];
rz(-2.3305927) q[0];
sx q[0];
rz(0.42914036) q[0];
rz(-2.7375713) q[1];
sx q[1];
rz(-0.41530135) q[1];
sx q[1];
rz(-2.2228352) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0149834) q[0];
sx q[0];
rz(-2.1169239) q[0];
sx q[0];
rz(-1.2864497) q[0];
rz(1.1303004) q[2];
sx q[2];
rz(-2.8545032) q[2];
sx q[2];
rz(-0.31781351) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.017581878) q[1];
sx q[1];
rz(-0.38330829) q[1];
sx q[1];
rz(-0.95335754) q[1];
rz(-pi) q[2];
rz(1.7946967) q[3];
sx q[3];
rz(-0.32156518) q[3];
sx q[3];
rz(2.0617776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1227485) q[2];
sx q[2];
rz(-2.4825373) q[2];
sx q[2];
rz(-1.3502236) q[2];
rz(-0.33251479) q[3];
sx q[3];
rz(-0.88034383) q[3];
sx q[3];
rz(-1.3688709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74828446) q[0];
sx q[0];
rz(-0.64993334) q[0];
sx q[0];
rz(-2.6848324) q[0];
rz(1.7204174) q[1];
sx q[1];
rz(-1.8616118) q[1];
sx q[1];
rz(-3.0866887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0676014) q[0];
sx q[0];
rz(-1.5419863) q[0];
sx q[0];
rz(2.3037698) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1403528) q[2];
sx q[2];
rz(-1.1932185) q[2];
sx q[2];
rz(2.0098639) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.12704472) q[1];
sx q[1];
rz(-2.603122) q[1];
sx q[1];
rz(-2.0035067) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53590448) q[3];
sx q[3];
rz(-1.2782989) q[3];
sx q[3];
rz(-1.3939987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2533337) q[2];
sx q[2];
rz(-2.4312225) q[2];
sx q[2];
rz(-0.19676512) q[2];
rz(-2.0678068) q[3];
sx q[3];
rz(-1.1353227) q[3];
sx q[3];
rz(-0.73608583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1573023) q[0];
sx q[0];
rz(-2.2881916) q[0];
sx q[0];
rz(-0.89944696) q[0];
rz(-0.31117123) q[1];
sx q[1];
rz(-1.9322461) q[1];
sx q[1];
rz(-1.7412294) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6104855) q[0];
sx q[0];
rz(-0.60043469) q[0];
sx q[0];
rz(-1.9470293) q[0];
rz(-pi) q[1];
rz(-0.24263361) q[2];
sx q[2];
rz(-1.0972692) q[2];
sx q[2];
rz(2.1514016) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3653917) q[1];
sx q[1];
rz(-1.1836056) q[1];
sx q[1];
rz(-2.4939818) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66026779) q[3];
sx q[3];
rz(-0.83383027) q[3];
sx q[3];
rz(1.6136606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0805936) q[2];
sx q[2];
rz(-0.62005764) q[2];
sx q[2];
rz(1.9798123) q[2];
rz(-1.0787841) q[3];
sx q[3];
rz(-0.99083841) q[3];
sx q[3];
rz(-2.259528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5921191) q[0];
sx q[0];
rz(-1.5220806) q[0];
sx q[0];
rz(1.4694389) q[0];
rz(-3.0081765) q[1];
sx q[1];
rz(-1.7620371) q[1];
sx q[1];
rz(-0.98099991) q[1];
rz(1.5518804) q[2];
sx q[2];
rz(-0.053413548) q[2];
sx q[2];
rz(1.5656768) q[2];
rz(2.6399163) q[3];
sx q[3];
rz(-1.6335084) q[3];
sx q[3];
rz(-0.15214534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
