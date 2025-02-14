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
rz(-2.5831232) q[0];
rz(0.59977579) q[1];
sx q[1];
rz(-1.7668626) q[1];
sx q[1];
rz(0.078628063) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0692908) q[0];
sx q[0];
rz(-2.3218394) q[0];
sx q[0];
rz(-0.79901521) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6958649) q[2];
sx q[2];
rz(-2.2457221) q[2];
sx q[2];
rz(-1.9927911) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3529571) q[1];
sx q[1];
rz(-2.87559) q[1];
sx q[1];
rz(2.8255844) q[1];
rz(-pi) q[2];
rz(0.73689245) q[3];
sx q[3];
rz(-0.46023638) q[3];
sx q[3];
rz(-2.5207507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2416396) q[2];
sx q[2];
rz(-3.0934379) q[2];
sx q[2];
rz(2.107035) q[2];
rz(-0.11040802) q[3];
sx q[3];
rz(-1.6610828) q[3];
sx q[3];
rz(-2.835527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13650525) q[0];
sx q[0];
rz(-1.445048) q[0];
sx q[0];
rz(1.1862296) q[0];
rz(-1.8738481) q[1];
sx q[1];
rz(-0.45919752) q[1];
sx q[1];
rz(1.3892106) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1339552) q[0];
sx q[0];
rz(-1.2651772) q[0];
sx q[0];
rz(-1.0108749) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6824426) q[2];
sx q[2];
rz(-2.0497339) q[2];
sx q[2];
rz(-0.079733036) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8858365) q[1];
sx q[1];
rz(-0.44609082) q[1];
sx q[1];
rz(0.85020937) q[1];
rz(-pi) q[2];
rz(-1.4382382) q[3];
sx q[3];
rz(-0.62843674) q[3];
sx q[3];
rz(2.4579163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9691951) q[2];
sx q[2];
rz(-1.1149422) q[2];
sx q[2];
rz(-1.4364852) q[2];
rz(-1.8819594) q[3];
sx q[3];
rz(-2.6862222) q[3];
sx q[3];
rz(-0.026738515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36419511) q[0];
sx q[0];
rz(-1.3854249) q[0];
sx q[0];
rz(-1.2170323) q[0];
rz(-0.2150391) q[1];
sx q[1];
rz(-1.8937078) q[1];
sx q[1];
rz(1.7020114) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3453044) q[0];
sx q[0];
rz(-0.6718967) q[0];
sx q[0];
rz(-0.20033269) q[0];
x q[1];
rz(1.7363988) q[2];
sx q[2];
rz(-0.51711997) q[2];
sx q[2];
rz(1.0696326) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9653757) q[1];
sx q[1];
rz(-1.6976408) q[1];
sx q[1];
rz(1.3398257) q[1];
x q[2];
rz(-2.2398364) q[3];
sx q[3];
rz(-0.34325156) q[3];
sx q[3];
rz(-2.1829829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3415459) q[2];
sx q[2];
rz(-1.505082) q[2];
sx q[2];
rz(-2.3254507) q[2];
rz(-3.1331565) q[3];
sx q[3];
rz(-0.71788994) q[3];
sx q[3];
rz(1.5754023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.7208045) q[0];
sx q[0];
rz(-1.1710465) q[0];
sx q[0];
rz(0.26914832) q[0];
rz(1.3828269) q[1];
sx q[1];
rz(-2.5803284) q[1];
sx q[1];
rz(-0.47163481) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4260226) q[0];
sx q[0];
rz(-1.4310808) q[0];
sx q[0];
rz(0.84265291) q[0];
rz(-pi) q[1];
x q[1];
rz(0.083375562) q[2];
sx q[2];
rz(-2.7365587) q[2];
sx q[2];
rz(-2.0138559) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.1822575) q[1];
sx q[1];
rz(-1.253009) q[1];
sx q[1];
rz(0.55303581) q[1];
rz(-pi) q[2];
rz(-1.1310234) q[3];
sx q[3];
rz(-0.71478715) q[3];
sx q[3];
rz(-0.21372114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4857594) q[2];
sx q[2];
rz(-1.7344319) q[2];
sx q[2];
rz(-1.2539697) q[2];
rz(0.29821011) q[3];
sx q[3];
rz(-1.2757653) q[3];
sx q[3];
rz(1.5378753) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1042079) q[0];
sx q[0];
rz(-2.2315114) q[0];
sx q[0];
rz(0.75697672) q[0];
rz(1.0589927) q[1];
sx q[1];
rz(-1.4451566) q[1];
sx q[1];
rz(-2.7991378) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22638527) q[0];
sx q[0];
rz(-1.7614375) q[0];
sx q[0];
rz(-1.9739499) q[0];
rz(-2.8169698) q[2];
sx q[2];
rz(-0.89541905) q[2];
sx q[2];
rz(1.7243232) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7389512) q[1];
sx q[1];
rz(-2.5243951) q[1];
sx q[1];
rz(0.67006858) q[1];
rz(-pi) q[2];
rz(-2.7259521) q[3];
sx q[3];
rz(-1.8438253) q[3];
sx q[3];
rz(1.331783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69514889) q[2];
sx q[2];
rz(-1.4651639) q[2];
sx q[2];
rz(0.62015074) q[2];
rz(-2.5946963) q[3];
sx q[3];
rz(-0.20709012) q[3];
sx q[3];
rz(2.0641522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81724375) q[0];
sx q[0];
rz(-0.20340782) q[0];
sx q[0];
rz(1.9203583) q[0];
rz(2.0381894) q[1];
sx q[1];
rz(-1.1627448) q[1];
sx q[1];
rz(0.81072909) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65064036) q[0];
sx q[0];
rz(-1.3155259) q[0];
sx q[0];
rz(1.7351236) q[0];
rz(2.3757908) q[2];
sx q[2];
rz(-0.81379902) q[2];
sx q[2];
rz(1.162093) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6891046) q[1];
sx q[1];
rz(-2.9924722) q[1];
sx q[1];
rz(1.945687) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2280336) q[3];
sx q[3];
rz(-2.676042) q[3];
sx q[3];
rz(-2.6014933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3409884) q[2];
sx q[2];
rz(-1.9438513) q[2];
sx q[2];
rz(-1.1654589) q[2];
rz(3.0321339) q[3];
sx q[3];
rz(-1.0215003) q[3];
sx q[3];
rz(-0.060733184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.3133746) q[0];
sx q[0];
rz(-0.010638588) q[0];
sx q[0];
rz(2.4995372) q[0];
rz(-2.8984046) q[1];
sx q[1];
rz(-0.41047341) q[1];
sx q[1];
rz(-0.64116716) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.084448) q[0];
sx q[0];
rz(-2.085745) q[0];
sx q[0];
rz(0.65339974) q[0];
rz(-pi) q[1];
rz(-1.2095318) q[2];
sx q[2];
rz(-2.2881977) q[2];
sx q[2];
rz(0.016005767) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6536246) q[1];
sx q[1];
rz(-1.89978) q[1];
sx q[1];
rz(-0.79923652) q[1];
rz(-2.9092714) q[3];
sx q[3];
rz(-1.9594176) q[3];
sx q[3];
rz(0.53750932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13333653) q[2];
sx q[2];
rz(-1.1587605) q[2];
sx q[2];
rz(-0.94375098) q[2];
rz(2.4933695) q[3];
sx q[3];
rz(-0.74130487) q[3];
sx q[3];
rz(-2.4751723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42565313) q[0];
sx q[0];
rz(-0.81099993) q[0];
sx q[0];
rz(0.42914036) q[0];
rz(-0.40402135) q[1];
sx q[1];
rz(-2.7262913) q[1];
sx q[1];
rz(0.9187575) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8480534) q[0];
sx q[0];
rz(-1.812879) q[0];
sx q[0];
rz(2.5771067) q[0];
rz(-1.1303004) q[2];
sx q[2];
rz(-2.8545032) q[2];
sx q[2];
rz(0.31781351) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4705953) q[1];
sx q[1];
rz(-1.2609224) q[1];
sx q[1];
rz(-2.9122374) q[1];
x q[2];
rz(-0.073831255) q[3];
sx q[3];
rz(-1.8840577) q[3];
sx q[3];
rz(1.3153803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.018844133) q[2];
sx q[2];
rz(-0.65905535) q[2];
sx q[2];
rz(-1.3502236) q[2];
rz(-2.8090779) q[3];
sx q[3];
rz(-2.2612488) q[3];
sx q[3];
rz(-1.3688709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3933082) q[0];
sx q[0];
rz(-0.64993334) q[0];
sx q[0];
rz(0.45676029) q[0];
rz(-1.7204174) q[1];
sx q[1];
rz(-1.8616118) q[1];
sx q[1];
rz(-0.054904003) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0739912) q[0];
sx q[0];
rz(-1.5996063) q[0];
sx q[0];
rz(2.3037698) q[0];
x q[1];
rz(-2.0012399) q[2];
sx q[2];
rz(-1.1932185) q[2];
sx q[2];
rz(1.1317288) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.12704472) q[1];
sx q[1];
rz(-2.603122) q[1];
sx q[1];
rz(1.1380859) q[1];
rz(0.53284455) q[3];
sx q[3];
rz(-0.60361629) q[3];
sx q[3];
rz(0.6288022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.88825893) q[2];
sx q[2];
rz(-2.4312225) q[2];
sx q[2];
rz(2.9448275) q[2];
rz(-2.0678068) q[3];
sx q[3];
rz(-2.00627) q[3];
sx q[3];
rz(-2.4055068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
rz(2.2421457) q[0];
rz(2.8304214) q[1];
sx q[1];
rz(-1.2093465) q[1];
sx q[1];
rz(1.7412294) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084564805) q[0];
sx q[0];
rz(-1.0175144) q[0];
sx q[0];
rz(2.8951089) q[0];
rz(-0.24263361) q[2];
sx q[2];
rz(-1.0972692) q[2];
sx q[2];
rz(2.1514016) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3653917) q[1];
sx q[1];
rz(-1.9579871) q[1];
sx q[1];
rz(0.64761083) q[1];
x q[2];
rz(-2.1651044) q[3];
sx q[3];
rz(-2.1955238) q[3];
sx q[3];
rz(2.4710666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.060999) q[2];
sx q[2];
rz(-0.62005764) q[2];
sx q[2];
rz(-1.1617804) q[2];
rz(-1.0787841) q[3];
sx q[3];
rz(-2.1507542) q[3];
sx q[3];
rz(2.259528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5494736) q[0];
sx q[0];
rz(-1.5220806) q[0];
sx q[0];
rz(1.4694389) q[0];
rz(0.1334162) q[1];
sx q[1];
rz(-1.7620371) q[1];
sx q[1];
rz(-0.98099991) q[1];
rz(-1.5173923) q[2];
sx q[2];
rz(-1.5718062) q[2];
sx q[2];
rz(0.013769416) q[2];
rz(-1.4992989) q[3];
sx q[3];
rz(-2.0713948) q[3];
sx q[3];
rz(-1.757302) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
