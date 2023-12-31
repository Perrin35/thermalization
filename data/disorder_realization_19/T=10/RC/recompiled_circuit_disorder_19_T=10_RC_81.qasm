OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1146381) q[0];
sx q[0];
rz(-1.4517598) q[0];
sx q[0];
rz(2.4960158) q[0];
rz(0.37880701) q[1];
sx q[1];
rz(4.9103476) q[1];
sx q[1];
rz(11.068439) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71292114) q[0];
sx q[0];
rz(-0.38654583) q[0];
sx q[0];
rz(-1.1791694) q[0];
rz(-2.2013118) q[2];
sx q[2];
rz(-1.6271546) q[2];
sx q[2];
rz(-1.8688569) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.77817569) q[1];
sx q[1];
rz(-2.2925804) q[1];
sx q[1];
rz(-1.282882) q[1];
rz(0.50819355) q[3];
sx q[3];
rz(-0.71532202) q[3];
sx q[3];
rz(3.1014266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2215185) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(2.8011838) q[2];
rz(2.3085964) q[3];
sx q[3];
rz(-2.4383128) q[3];
sx q[3];
rz(-1.7849281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4689363) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(-0.54498589) q[0];
rz(2.2333721) q[1];
sx q[1];
rz(-2.814099) q[1];
sx q[1];
rz(2.3166336) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93281125) q[0];
sx q[0];
rz(-2.6926846) q[0];
sx q[0];
rz(-2.0138028) q[0];
rz(-pi) q[1];
rz(1.5812133) q[2];
sx q[2];
rz(-1.9736819) q[2];
sx q[2];
rz(1.252623) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7823239) q[1];
sx q[1];
rz(-1.6869079) q[1];
sx q[1];
rz(-0.5639204) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0042079) q[3];
sx q[3];
rz(-1.7363747) q[3];
sx q[3];
rz(2.2752938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5543582) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(2.9651802) q[2];
rz(-0.13088626) q[3];
sx q[3];
rz(-1.1754879) q[3];
sx q[3];
rz(3.1089354) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762887) q[0];
sx q[0];
rz(-0.58350199) q[0];
sx q[0];
rz(2.6718455) q[0];
rz(-1.5247955) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(3.1030531) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24667106) q[0];
sx q[0];
rz(-1.5630432) q[0];
sx q[0];
rz(0.00064108032) q[0];
rz(0.23080319) q[2];
sx q[2];
rz(-1.7204086) q[2];
sx q[2];
rz(-0.19300592) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.027027834) q[1];
sx q[1];
rz(-1.6323043) q[1];
sx q[1];
rz(-1.0807178) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1707118) q[3];
sx q[3];
rz(-1.0661134) q[3];
sx q[3];
rz(-1.2952627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.58275756) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(2.2290686) q[2];
rz(1.3085261) q[3];
sx q[3];
rz(-2.1378744) q[3];
sx q[3];
rz(-1.7002038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38347605) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(2.7918949) q[0];
rz(-1.8967459) q[1];
sx q[1];
rz(-2.8379776) q[1];
sx q[1];
rz(-2.9464338) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0450268) q[0];
sx q[0];
rz(-1.4782527) q[0];
sx q[0];
rz(1.4499614) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85907016) q[2];
sx q[2];
rz(-2.4675183) q[2];
sx q[2];
rz(0.40875834) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3217722) q[1];
sx q[1];
rz(-1.1204549) q[1];
sx q[1];
rz(-0.26178534) q[1];
rz(0.021899453) q[3];
sx q[3];
rz(-2.5023513) q[3];
sx q[3];
rz(0.99614776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.23094709) q[2];
sx q[2];
rz(-0.4898943) q[2];
sx q[2];
rz(-1.267743) q[2];
rz(-2.0729444) q[3];
sx q[3];
rz(-2.1290776) q[3];
sx q[3];
rz(-0.8403362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1068263) q[0];
sx q[0];
rz(-1.4522469) q[0];
sx q[0];
rz(3.0019794) q[0];
rz(-1.0768249) q[1];
sx q[1];
rz(-1.040841) q[1];
sx q[1];
rz(-2.7635014) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17079167) q[0];
sx q[0];
rz(-2.1149201) q[0];
sx q[0];
rz(-2.8118954) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81850448) q[2];
sx q[2];
rz(-0.92220014) q[2];
sx q[2];
rz(1.4865781) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7996019) q[1];
sx q[1];
rz(-1.2050036) q[1];
sx q[1];
rz(-1.1127383) q[1];
rz(-2.7630745) q[3];
sx q[3];
rz(-2.6444728) q[3];
sx q[3];
rz(-2.5147223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2698764) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(0.88796973) q[2];
rz(0.97638431) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(-2.1358657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3865005) q[0];
sx q[0];
rz(-0.077682406) q[0];
sx q[0];
rz(0.37762541) q[0];
rz(-0.32304421) q[1];
sx q[1];
rz(-0.66345614) q[1];
sx q[1];
rz(2.2264218) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5765502) q[0];
sx q[0];
rz(-0.33320198) q[0];
sx q[0];
rz(-2.8898426) q[0];
rz(-pi) q[1];
rz(2.633811) q[2];
sx q[2];
rz(-2.0785329) q[2];
sx q[2];
rz(-2.7024262) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.311971) q[1];
sx q[1];
rz(-1.258007) q[1];
sx q[1];
rz(-1.642986) q[1];
rz(-0.69453199) q[3];
sx q[3];
rz(-1.6581222) q[3];
sx q[3];
rz(-0.15184034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6043828) q[2];
sx q[2];
rz(-2.9769124) q[2];
sx q[2];
rz(-0.79052314) q[2];
rz(0.28997713) q[3];
sx q[3];
rz(-2.4089549) q[3];
sx q[3];
rz(-2.524232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4054366) q[0];
sx q[0];
rz(-1.0671395) q[0];
sx q[0];
rz(-2.7745568) q[0];
rz(1.908318) q[1];
sx q[1];
rz(-1.9676625) q[1];
sx q[1];
rz(-1.7162011) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0910949) q[0];
sx q[0];
rz(-2.0031345) q[0];
sx q[0];
rz(-1.7763441) q[0];
x q[1];
rz(-0.44600365) q[2];
sx q[2];
rz(-1.5127231) q[2];
sx q[2];
rz(2.8446978) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.094147625) q[1];
sx q[1];
rz(-0.74456753) q[1];
sx q[1];
rz(0.57569699) q[1];
rz(-1.5464877) q[3];
sx q[3];
rz(-2.6377502) q[3];
sx q[3];
rz(-2.3712096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1620862) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(-1.210775) q[2];
rz(3.1397505) q[3];
sx q[3];
rz(-2.3760934) q[3];
sx q[3];
rz(0.65729284) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84412557) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(-0.076518245) q[0];
rz(2.466295) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(-0.75235596) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9936179) q[0];
sx q[0];
rz(-0.60950845) q[0];
sx q[0];
rz(-2.5698635) q[0];
rz(-pi) q[1];
rz(-0.1444507) q[2];
sx q[2];
rz(-0.92717294) q[2];
sx q[2];
rz(0.87063906) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89531089) q[1];
sx q[1];
rz(-1.6591676) q[1];
sx q[1];
rz(-1.6842711) q[1];
rz(-0.30267834) q[3];
sx q[3];
rz(-2.5303826) q[3];
sx q[3];
rz(-1.6906893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.9383119) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(0.00014649815) q[2];
rz(-1.1095307) q[3];
sx q[3];
rz(-2.2347982) q[3];
sx q[3];
rz(2.8439567) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14934854) q[0];
sx q[0];
rz(-2.902817) q[0];
sx q[0];
rz(-2.1355656) q[0];
rz(-0.26793119) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(-1.8267652) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49004236) q[0];
sx q[0];
rz(-2.393912) q[0];
sx q[0];
rz(2.7689395) q[0];
x q[1];
rz(2.7776412) q[2];
sx q[2];
rz(-1.6953354) q[2];
sx q[2];
rz(0.044791128) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6445551) q[1];
sx q[1];
rz(-0.92004787) q[1];
sx q[1];
rz(-0.52849309) q[1];
x q[2];
rz(-2.3141765) q[3];
sx q[3];
rz(-1.5762868) q[3];
sx q[3];
rz(-1.0178125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7018147) q[2];
sx q[2];
rz(-2.5952227) q[2];
sx q[2];
rz(2.8961199) q[2];
rz(-0.43073511) q[3];
sx q[3];
rz(-1.0916748) q[3];
sx q[3];
rz(2.664264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7219287) q[0];
sx q[0];
rz(-2.1491282) q[0];
sx q[0];
rz(0.28668177) q[0];
rz(0.87896705) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(-0.39189664) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.058404) q[0];
sx q[0];
rz(-1.579111) q[0];
sx q[0];
rz(-2.1283098) q[0];
x q[1];
rz(3.0060843) q[2];
sx q[2];
rz(-0.91943179) q[2];
sx q[2];
rz(1.1609921) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74589409) q[1];
sx q[1];
rz(-0.77334009) q[1];
sx q[1];
rz(-1.3032773) q[1];
x q[2];
rz(-2.5819671) q[3];
sx q[3];
rz(-1.2979753) q[3];
sx q[3];
rz(0.19459693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.60951704) q[2];
sx q[2];
rz(-2.1485907) q[2];
sx q[2];
rz(-1.9840476) q[2];
rz(2.5907497) q[3];
sx q[3];
rz(-1.563787) q[3];
sx q[3];
rz(1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3863603) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(1.339284) q[1];
sx q[1];
rz(-2.5279999) q[1];
sx q[1];
rz(-2.7816714) q[1];
rz(-2.3618868) q[2];
sx q[2];
rz(-1.4479326) q[2];
sx q[2];
rz(-0.89276199) q[2];
rz(-2.7950381) q[3];
sx q[3];
rz(-1.0481491) q[3];
sx q[3];
rz(-0.92878503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
