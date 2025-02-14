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
rz(2.9146258) q[0];
sx q[0];
rz(-0.6114971) q[0];
sx q[0];
rz(0.55846941) q[0];
rz(0.59977579) q[1];
sx q[1];
rz(-1.7668626) q[1];
sx q[1];
rz(0.078628063) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0315379) q[0];
sx q[0];
rz(-1.0194091) q[0];
sx q[0];
rz(0.64161513) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67875454) q[2];
sx q[2];
rz(-1.6683443) q[2];
sx q[2];
rz(0.3435979) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0262194) q[1];
sx q[1];
rz(-1.3182606) q[1];
sx q[1];
rz(1.6552684) q[1];
rz(1.2492248) q[3];
sx q[3];
rz(-1.235629) q[3];
sx q[3];
rz(-0.17091076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2416396) q[2];
sx q[2];
rz(-3.0934379) q[2];
sx q[2];
rz(1.0345577) q[2];
rz(-0.11040802) q[3];
sx q[3];
rz(-1.4805099) q[3];
sx q[3];
rz(-0.30606562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0050874) q[0];
sx q[0];
rz(-1.445048) q[0];
sx q[0];
rz(-1.1862296) q[0];
rz(1.8738481) q[1];
sx q[1];
rz(-0.45919752) q[1];
sx q[1];
rz(1.7523821) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8911676) q[0];
sx q[0];
rz(-2.1019263) q[0];
sx q[0];
rz(-2.7851339) q[0];
rz(-2.6600962) q[2];
sx q[2];
rz(-1.6698368) q[2];
sx q[2];
rz(-1.5426829) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4965044) q[1];
sx q[1];
rz(-1.2821272) q[1];
sx q[1];
rz(-1.225807) q[1];
rz(-pi) q[2];
rz(1.4382382) q[3];
sx q[3];
rz(-2.5131559) q[3];
sx q[3];
rz(-0.6836764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17239751) q[2];
sx q[2];
rz(-2.0266504) q[2];
sx q[2];
rz(1.7051075) q[2];
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
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36419511) q[0];
sx q[0];
rz(-1.7561678) q[0];
sx q[0];
rz(1.9245603) q[0];
rz(0.2150391) q[1];
sx q[1];
rz(-1.8937078) q[1];
sx q[1];
rz(-1.7020114) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091451784) q[0];
sx q[0];
rz(-0.9147075) q[0];
sx q[0];
rz(1.4138282) q[0];
rz(-0.093482253) q[2];
sx q[2];
rz(-2.0801525) q[2];
sx q[2];
rz(-0.87967726) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9653757) q[1];
sx q[1];
rz(-1.6976408) q[1];
sx q[1];
rz(-1.3398257) q[1];
rz(2.2398364) q[3];
sx q[3];
rz(-2.7983411) q[3];
sx q[3];
rz(0.95860976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.8000468) q[2];
sx q[2];
rz(-1.6365106) q[2];
sx q[2];
rz(-0.81614196) q[2];
rz(-3.1331565) q[3];
sx q[3];
rz(-0.71788994) q[3];
sx q[3];
rz(1.5754023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7208045) q[0];
sx q[0];
rz(-1.9705462) q[0];
sx q[0];
rz(-0.26914832) q[0];
rz(1.7587657) q[1];
sx q[1];
rz(-2.5803284) q[1];
sx q[1];
rz(-2.6699578) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29971805) q[0];
sx q[0];
rz(-0.73900676) q[0];
sx q[0];
rz(-1.7790545) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5351065) q[2];
sx q[2];
rz(-1.9743414) q[2];
sx q[2];
rz(-1.2184136) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5790807) q[1];
sx q[1];
rz(-1.0483841) q[1];
sx q[1];
rz(1.201931) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7876623) q[3];
sx q[3];
rz(-0.9359064) q[3];
sx q[3];
rz(-2.3707182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.65583324) q[2];
sx q[2];
rz(-1.4071608) q[2];
sx q[2];
rz(-1.2539697) q[2];
rz(2.8433825) q[3];
sx q[3];
rz(-1.2757653) q[3];
sx q[3];
rz(1.6037174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1042079) q[0];
sx q[0];
rz(-2.2315114) q[0];
sx q[0];
rz(2.3846159) q[0];
rz(-2.0825999) q[1];
sx q[1];
rz(-1.4451566) q[1];
sx q[1];
rz(0.34245488) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9152074) q[0];
sx q[0];
rz(-1.7614375) q[0];
sx q[0];
rz(1.9739499) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2724897) q[2];
sx q[2];
rz(-1.822374) q[2];
sx q[2];
rz(3.0877496) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7389512) q[1];
sx q[1];
rz(-0.61719751) q[1];
sx q[1];
rz(0.67006858) q[1];
x q[2];
rz(0.41564055) q[3];
sx q[3];
rz(-1.2977674) q[3];
sx q[3];
rz(1.8098097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4464438) q[2];
sx q[2];
rz(-1.4651639) q[2];
sx q[2];
rz(-2.5214419) q[2];
rz(2.5946963) q[3];
sx q[3];
rz(-2.9345025) q[3];
sx q[3];
rz(-1.0774405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81724375) q[0];
sx q[0];
rz(-2.9381848) q[0];
sx q[0];
rz(1.2212344) q[0];
rz(-1.1034032) q[1];
sx q[1];
rz(-1.1627448) q[1];
sx q[1];
rz(0.81072909) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2632836) q[0];
sx q[0];
rz(-1.4118402) q[0];
sx q[0];
rz(-2.8829888) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3757908) q[2];
sx q[2];
rz(-0.81379902) q[2];
sx q[2];
rz(1.162093) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6891046) q[1];
sx q[1];
rz(-0.14912046) q[1];
sx q[1];
rz(1.945687) q[1];
x q[2];
rz(1.1922513) q[3];
sx q[3];
rz(-1.2929799) q[3];
sx q[3];
rz(-1.5072106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3409884) q[2];
sx q[2];
rz(-1.9438513) q[2];
sx q[2];
rz(1.1654589) q[2];
rz(0.10945877) q[3];
sx q[3];
rz(-2.1200924) q[3];
sx q[3];
rz(3.0808595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8282181) q[0];
sx q[0];
rz(-3.1309541) q[0];
sx q[0];
rz(-2.4995372) q[0];
rz(-2.8984046) q[1];
sx q[1];
rz(-0.41047341) q[1];
sx q[1];
rz(-0.64116716) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0571447) q[0];
sx q[0];
rz(-2.085745) q[0];
sx q[0];
rz(-0.65339974) q[0];
x q[1];
rz(1.9320609) q[2];
sx q[2];
rz(-2.2881977) q[2];
sx q[2];
rz(0.016005767) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7542759) q[1];
sx q[1];
rz(-2.2914304) q[1];
sx q[1];
rz(0.44447036) q[1];
rz(-1.9690597) q[3];
sx q[3];
rz(-1.7855111) q[3];
sx q[3];
rz(2.1977148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.13333653) q[2];
sx q[2];
rz(-1.9828321) q[2];
sx q[2];
rz(0.94375098) q[2];
rz(-0.64822316) q[3];
sx q[3];
rz(-0.74130487) q[3];
sx q[3];
rz(-2.4751723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42565313) q[0];
sx q[0];
rz(-2.3305927) q[0];
sx q[0];
rz(0.42914036) q[0];
rz(-2.7375713) q[1];
sx q[1];
rz(-2.7262913) q[1];
sx q[1];
rz(2.2228352) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63916535) q[0];
sx q[0];
rz(-0.60898655) q[0];
sx q[0];
rz(-0.43242411) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8317675) q[2];
sx q[2];
rz(-1.6918285) q[2];
sx q[2];
rz(0.82838917) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1707216) q[1];
sx q[1];
rz(-1.7890463) q[1];
sx q[1];
rz(1.2531325) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0677614) q[3];
sx q[3];
rz(-1.8840577) q[3];
sx q[3];
rz(1.3153803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1227485) q[2];
sx q[2];
rz(-2.4825373) q[2];
sx q[2];
rz(-1.3502236) q[2];
rz(0.33251479) q[3];
sx q[3];
rz(-2.2612488) q[3];
sx q[3];
rz(-1.3688709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3933082) q[0];
sx q[0];
rz(-0.64993334) q[0];
sx q[0];
rz(0.45676029) q[0];
rz(1.4211753) q[1];
sx q[1];
rz(-1.8616118) q[1];
sx q[1];
rz(-0.054904003) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0676014) q[0];
sx q[0];
rz(-1.5419863) q[0];
sx q[0];
rz(0.83782283) q[0];
rz(-pi) q[1];
x q[1];
rz(2.730092) q[2];
sx q[2];
rz(-1.1724645) q[2];
sx q[2];
rz(0.27137941) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5209321) q[1];
sx q[1];
rz(-2.0550613) q[1];
sx q[1];
rz(2.8961532) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.233912) q[3];
sx q[3];
rz(-2.0816605) q[3];
sx q[3];
rz(0.0072016933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.88825893) q[2];
sx q[2];
rz(-2.4312225) q[2];
sx q[2];
rz(-0.19676512) q[2];
rz(1.0737859) q[3];
sx q[3];
rz(-1.1353227) q[3];
sx q[3];
rz(-0.73608583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1573023) q[0];
sx q[0];
rz(-2.2881916) q[0];
sx q[0];
rz(2.2421457) q[0];
rz(-2.8304214) q[1];
sx q[1];
rz(-1.9322461) q[1];
sx q[1];
rz(-1.4003632) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.78681) q[0];
sx q[0];
rz(-1.3616832) q[0];
sx q[0];
rz(2.1378921) q[0];
rz(-pi) q[1];
rz(-0.24263361) q[2];
sx q[2];
rz(-1.0972692) q[2];
sx q[2];
rz(-0.99019105) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33161641) q[1];
sx q[1];
rz(-0.73996937) q[1];
sx q[1];
rz(-0.59438594) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71615852) q[3];
sx q[3];
rz(-1.0993488) q[3];
sx q[3];
rz(2.6177419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.060999) q[2];
sx q[2];
rz(-0.62005764) q[2];
sx q[2];
rz(-1.9798123) q[2];
rz(-2.0628085) q[3];
sx q[3];
rz(-2.1507542) q[3];
sx q[3];
rz(0.88206464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5921191) q[0];
sx q[0];
rz(-1.5220806) q[0];
sx q[0];
rz(1.4694389) q[0];
rz(-0.1334162) q[1];
sx q[1];
rz(-1.3795556) q[1];
sx q[1];
rz(2.1605927) q[1];
rz(1.5897122) q[2];
sx q[2];
rz(-3.0881791) q[2];
sx q[2];
rz(-1.5759158) q[2];
rz(3.0117494) q[3];
sx q[3];
rz(-2.6363439) q[3];
sx q[3];
rz(1.5324203) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
