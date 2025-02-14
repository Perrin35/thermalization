OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1898243) q[0];
sx q[0];
rz(-0.70280743) q[0];
sx q[0];
rz(1.2220609) q[0];
rz(-3.1932073) q[1];
sx q[1];
rz(1.5049223) q[1];
sx q[1];
rz(7.7636889) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5334789) q[0];
sx q[0];
rz(-2.2256269) q[0];
sx q[0];
rz(2.9346908) q[0];
rz(-pi) q[1];
rz(1.0009345) q[2];
sx q[2];
rz(-2.9409932) q[2];
sx q[2];
rz(-0.72884411) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6340629) q[1];
sx q[1];
rz(-0.46137091) q[1];
sx q[1];
rz(-0.4005691) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2961404) q[3];
sx q[3];
rz(-0.61224711) q[3];
sx q[3];
rz(-2.1795764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1808971) q[2];
sx q[2];
rz(-0.81520671) q[2];
sx q[2];
rz(-1.0308456) q[2];
rz(-0.24813063) q[3];
sx q[3];
rz(-1.7957567) q[3];
sx q[3];
rz(-2.8743751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88028518) q[0];
sx q[0];
rz(-1.945865) q[0];
sx q[0];
rz(-1.1241359) q[0];
rz(-1.0719489) q[1];
sx q[1];
rz(-1.5771461) q[1];
sx q[1];
rz(-0.42676485) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.916722) q[0];
sx q[0];
rz(-1.568075) q[0];
sx q[0];
rz(0.82758521) q[0];
rz(-pi) q[1];
rz(0.63964455) q[2];
sx q[2];
rz(-0.42176127) q[2];
sx q[2];
rz(2.1157921) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2093231) q[1];
sx q[1];
rz(-0.70581064) q[1];
sx q[1];
rz(2.9747444) q[1];
x q[2];
rz(-1.9598448) q[3];
sx q[3];
rz(-1.3836244) q[3];
sx q[3];
rz(2.1120918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.65853226) q[2];
sx q[2];
rz(-2.6250562) q[2];
sx q[2];
rz(0.10022441) q[2];
rz(2.0641067) q[3];
sx q[3];
rz(-1.5651549) q[3];
sx q[3];
rz(-1.235435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60233068) q[0];
sx q[0];
rz(-2.1644008) q[0];
sx q[0];
rz(1.7072898) q[0];
rz(-0.82646787) q[1];
sx q[1];
rz(-1.4974599) q[1];
sx q[1];
rz(2.7109587) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5204091) q[0];
sx q[0];
rz(-0.31776515) q[0];
sx q[0];
rz(0.54917021) q[0];
x q[1];
rz(-1.1330092) q[2];
sx q[2];
rz(-0.29972813) q[2];
sx q[2];
rz(2.9309811) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0890177) q[1];
sx q[1];
rz(-2.3529426) q[1];
sx q[1];
rz(2.1441205) q[1];
x q[2];
rz(-2.1149733) q[3];
sx q[3];
rz(-1.9154356) q[3];
sx q[3];
rz(1.0349864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2677801) q[2];
sx q[2];
rz(-1.1649106) q[2];
sx q[2];
rz(2.5679892) q[2];
rz(2.7576647) q[3];
sx q[3];
rz(-1.826518) q[3];
sx q[3];
rz(-0.65281502) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7225994) q[0];
sx q[0];
rz(-2.3498131) q[0];
sx q[0];
rz(2.0844039) q[0];
rz(-2.3253697) q[1];
sx q[1];
rz(-2.1482601) q[1];
sx q[1];
rz(2.9073471) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5200096) q[0];
sx q[0];
rz(-2.710973) q[0];
sx q[0];
rz(-1.8355811) q[0];
rz(-pi) q[1];
rz(-2.7170638) q[2];
sx q[2];
rz(-2.0574951) q[2];
sx q[2];
rz(-0.27482143) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2709759) q[1];
sx q[1];
rz(-0.91667316) q[1];
sx q[1];
rz(0.65947094) q[1];
x q[2];
rz(0.81984249) q[3];
sx q[3];
rz(-2.4944135) q[3];
sx q[3];
rz(2.1539713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.28698841) q[2];
sx q[2];
rz(-1.333326) q[2];
sx q[2];
rz(2.2387779) q[2];
rz(1.758894) q[3];
sx q[3];
rz(-0.32302502) q[3];
sx q[3];
rz(2.2949016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6803902) q[0];
sx q[0];
rz(-1.2795376) q[0];
sx q[0];
rz(0.70988208) q[0];
rz(-0.4711802) q[1];
sx q[1];
rz(-1.2268365) q[1];
sx q[1];
rz(3.0772298) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.775255) q[0];
sx q[0];
rz(-1.4357872) q[0];
sx q[0];
rz(2.058821) q[0];
x q[1];
rz(2.0662599) q[2];
sx q[2];
rz(-1.3643985) q[2];
sx q[2];
rz(-1.0546233) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.086972039) q[1];
sx q[1];
rz(-0.68212648) q[1];
sx q[1];
rz(2.4404018) q[1];
x q[2];
rz(-3.0547906) q[3];
sx q[3];
rz(-1.5802812) q[3];
sx q[3];
rz(-0.61244165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.71089661) q[2];
sx q[2];
rz(-2.5888093) q[2];
sx q[2];
rz(2.1882679) q[2];
rz(2.583875) q[3];
sx q[3];
rz(-1.6744813) q[3];
sx q[3];
rz(-2.6127156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8271269) q[0];
sx q[0];
rz(-1.0449266) q[0];
sx q[0];
rz(-0.99660981) q[0];
rz(3.0820471) q[1];
sx q[1];
rz(-0.98809067) q[1];
sx q[1];
rz(-1.3522118) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23350643) q[0];
sx q[0];
rz(-0.75859374) q[0];
sx q[0];
rz(-0.80484606) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3699721) q[2];
sx q[2];
rz(-1.3305656) q[2];
sx q[2];
rz(0.26773237) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.775573) q[1];
sx q[1];
rz(-1.8692353) q[1];
sx q[1];
rz(-1.8587023) q[1];
rz(2.7963403) q[3];
sx q[3];
rz(-2.1998029) q[3];
sx q[3];
rz(-2.0382413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.754564) q[2];
sx q[2];
rz(-1.2754385) q[2];
sx q[2];
rz(-0.55956364) q[2];
rz(2.3061421) q[3];
sx q[3];
rz(-2.7159034) q[3];
sx q[3];
rz(1.3028418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4556722) q[0];
sx q[0];
rz(-2.4113825) q[0];
sx q[0];
rz(0.36668229) q[0];
rz(-2.2600251) q[1];
sx q[1];
rz(-2.5036948) q[1];
sx q[1];
rz(1.4311904) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87302215) q[0];
sx q[0];
rz(-2.2013469) q[0];
sx q[0];
rz(2.0591169) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9164247) q[2];
sx q[2];
rz(-2.2273438) q[2];
sx q[2];
rz(-0.57166568) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.28364946) q[1];
sx q[1];
rz(-1.4542631) q[1];
sx q[1];
rz(1.5416508) q[1];
rz(-pi) q[2];
rz(-0.65815417) q[3];
sx q[3];
rz(-1.3835211) q[3];
sx q[3];
rz(1.8204456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35393474) q[2];
sx q[2];
rz(-2.6255609) q[2];
sx q[2];
rz(2.3616974) q[2];
rz(-2.6025313) q[3];
sx q[3];
rz(-1.6293679) q[3];
sx q[3];
rz(-2.5299634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.014538177) q[0];
sx q[0];
rz(-2.4712565) q[0];
sx q[0];
rz(-0.77914733) q[0];
rz(-0.51891333) q[1];
sx q[1];
rz(-1.8592368) q[1];
sx q[1];
rz(-0.40294495) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3765744) q[0];
sx q[0];
rz(-0.85581644) q[0];
sx q[0];
rz(-0.9299703) q[0];
rz(-pi) q[1];
rz(2.05415) q[2];
sx q[2];
rz(-1.8515808) q[2];
sx q[2];
rz(-1.3705685) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1379304) q[1];
sx q[1];
rz(-2.3786585) q[1];
sx q[1];
rz(2.5502045) q[1];
rz(-pi) q[2];
rz(-1.9698079) q[3];
sx q[3];
rz(-2.5305989) q[3];
sx q[3];
rz(-0.64507592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.11401033) q[2];
sx q[2];
rz(-1.9738013) q[2];
sx q[2];
rz(-2.0243417) q[2];
rz(2.4452325) q[3];
sx q[3];
rz(-2.0317234) q[3];
sx q[3];
rz(-1.2874359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6663412) q[0];
sx q[0];
rz(-0.23463686) q[0];
sx q[0];
rz(-0.41467211) q[0];
rz(2.0798202) q[1];
sx q[1];
rz(-0.88645187) q[1];
sx q[1];
rz(-0.83962238) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79796004) q[0];
sx q[0];
rz(-1.1251483) q[0];
sx q[0];
rz(-2.6889422) q[0];
x q[1];
rz(-2.1346852) q[2];
sx q[2];
rz(-2.208935) q[2];
sx q[2];
rz(-1.9319122) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.190358) q[1];
sx q[1];
rz(-0.91427416) q[1];
sx q[1];
rz(1.7503529) q[1];
rz(1.8168338) q[3];
sx q[3];
rz(-1.2241505) q[3];
sx q[3];
rz(2.3510166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7825369) q[2];
sx q[2];
rz(-1.3691207) q[2];
sx q[2];
rz(-3.0699406) q[2];
rz(3.0960633) q[3];
sx q[3];
rz(-1.1979016) q[3];
sx q[3];
rz(0.45904747) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59938207) q[0];
sx q[0];
rz(-2.5686503) q[0];
sx q[0];
rz(-2.086916) q[0];
rz(-0.032546267) q[1];
sx q[1];
rz(-0.97144214) q[1];
sx q[1];
rz(-0.5221101) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3536026) q[0];
sx q[0];
rz(-0.66594571) q[0];
sx q[0];
rz(2.6329726) q[0];
x q[1];
rz(2.1565151) q[2];
sx q[2];
rz(-2.217948) q[2];
sx q[2];
rz(-3.0708608) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.13455277) q[1];
sx q[1];
rz(-1.1940446) q[1];
sx q[1];
rz(1.1621848) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2100092) q[3];
sx q[3];
rz(-0.57769247) q[3];
sx q[3];
rz(-3.011774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1471499) q[2];
sx q[2];
rz(-2.9394737) q[2];
sx q[2];
rz(1.7333376) q[2];
rz(-1.5857006) q[3];
sx q[3];
rz(-0.2318016) q[3];
sx q[3];
rz(-2.2304992) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0178575) q[0];
sx q[0];
rz(-1.8386848) q[0];
sx q[0];
rz(0.89717502) q[0];
rz(-0.97902117) q[1];
sx q[1];
rz(-1.9985825) q[1];
sx q[1];
rz(-0.40252007) q[1];
rz(-1.2382837) q[2];
sx q[2];
rz(-1.7560183) q[2];
sx q[2];
rz(2.8011372) q[2];
rz(-1.3835945) q[3];
sx q[3];
rz(-0.91619195) q[3];
sx q[3];
rz(-0.34294101) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
