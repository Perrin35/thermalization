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
rz(-1.4323956) q[0];
sx q[0];
rz(-0.30269912) q[0];
sx q[0];
rz(1.0330638) q[0];
rz(-1.2248224) q[1];
sx q[1];
rz(-1.4900102) q[1];
sx q[1];
rz(2.9472247) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5759566) q[0];
sx q[0];
rz(-2.4969184) q[0];
sx q[0];
rz(-2.8863532) q[0];
x q[1];
rz(-0.69673522) q[2];
sx q[2];
rz(-1.8039304) q[2];
sx q[2];
rz(-1.8270884) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7558507) q[1];
sx q[1];
rz(-3.0839018) q[1];
sx q[1];
rz(2.2123446) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1367873) q[3];
sx q[3];
rz(-2.1510486) q[3];
sx q[3];
rz(-0.88742076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4091461) q[2];
sx q[2];
rz(-2.0398085) q[2];
sx q[2];
rz(2.5771602) q[2];
rz(1.268528) q[3];
sx q[3];
rz(-0.0081491834) q[3];
sx q[3];
rz(2.234999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059435189) q[0];
sx q[0];
rz(-3.1340288) q[0];
sx q[0];
rz(2.2285158) q[0];
rz(-0.68436855) q[1];
sx q[1];
rz(-3.1412536) q[1];
sx q[1];
rz(-0.93904644) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53227753) q[0];
sx q[0];
rz(-1.0160722) q[0];
sx q[0];
rz(2.0164665) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46397658) q[2];
sx q[2];
rz(-1.5590406) q[2];
sx q[2];
rz(-1.8092312) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1176743) q[1];
sx q[1];
rz(-0.54421762) q[1];
sx q[1];
rz(-0.88440374) q[1];
rz(-pi) q[2];
rz(-1.7250118) q[3];
sx q[3];
rz(-0.41114488) q[3];
sx q[3];
rz(-1.4158885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82274503) q[2];
sx q[2];
rz(-3.1343967) q[2];
sx q[2];
rz(-0.89246559) q[2];
rz(-1.578963) q[3];
sx q[3];
rz(-3.1171418) q[3];
sx q[3];
rz(1.310937) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40176323) q[0];
sx q[0];
rz(-0.50245291) q[0];
sx q[0];
rz(-0.40973642) q[0];
rz(-3.1306664) q[1];
sx q[1];
rz(-2.9210563) q[1];
sx q[1];
rz(1.797537) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6861658) q[0];
sx q[0];
rz(-2.6302105) q[0];
sx q[0];
rz(-1.5288985) q[0];
rz(-1.686269) q[2];
sx q[2];
rz(-1.5073379) q[2];
sx q[2];
rz(-2.7977914) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2933767) q[1];
sx q[1];
rz(-0.266222) q[1];
sx q[1];
rz(1.3220399) q[1];
rz(-pi) q[2];
rz(-1.4555172) q[3];
sx q[3];
rz(-1.5805827) q[3];
sx q[3];
rz(-0.20231314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.672013) q[2];
sx q[2];
rz(-1.6894222) q[2];
sx q[2];
rz(-3.0984042) q[2];
rz(2.0016661) q[3];
sx q[3];
rz(-2.98525) q[3];
sx q[3];
rz(-1.9388916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4129055) q[0];
sx q[0];
rz(-2.7396956) q[0];
sx q[0];
rz(-1.3380949) q[0];
rz(-0.69522229) q[1];
sx q[1];
rz(-0.1278563) q[1];
sx q[1];
rz(2.7241838) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.886717) q[0];
sx q[0];
rz(-2.4435197) q[0];
sx q[0];
rz(-1.8398192) q[0];
x q[1];
rz(1.2222671) q[2];
sx q[2];
rz(-0.67643316) q[2];
sx q[2];
rz(-1.2471022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7898832) q[1];
sx q[1];
rz(-1.369006) q[1];
sx q[1];
rz(2.2191597) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.999971) q[3];
sx q[3];
rz(-1.9998261) q[3];
sx q[3];
rz(-2.1528368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7030316) q[2];
sx q[2];
rz(-0.87623864) q[2];
sx q[2];
rz(-1.7523127) q[2];
rz(-2.321068) q[3];
sx q[3];
rz(-1.5503784) q[3];
sx q[3];
rz(-2.4337721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97838068) q[0];
sx q[0];
rz(-0.037540171) q[0];
sx q[0];
rz(0.96330825) q[0];
rz(-2.9523201) q[1];
sx q[1];
rz(-0.015733868) q[1];
sx q[1];
rz(2.9761369) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4625797) q[0];
sx q[0];
rz(-1.5272584) q[0];
sx q[0];
rz(-1.5971558) q[0];
rz(-pi) q[1];
rz(0.29954977) q[2];
sx q[2];
rz(-0.89847696) q[2];
sx q[2];
rz(-1.7393665) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.206659) q[1];
sx q[1];
rz(-2.3689785) q[1];
sx q[1];
rz(1.4040703) q[1];
rz(-pi) q[2];
rz(0.9521552) q[3];
sx q[3];
rz(-2.6252271) q[3];
sx q[3];
rz(-1.0960033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8339707) q[2];
sx q[2];
rz(-2.622719) q[2];
sx q[2];
rz(-1.0597672) q[2];
rz(-0.69093949) q[3];
sx q[3];
rz(-2.2773404) q[3];
sx q[3];
rz(-1.8701657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6178013) q[0];
sx q[0];
rz(-0.093136223) q[0];
sx q[0];
rz(-1.5366489) q[0];
rz(2.6887584) q[1];
sx q[1];
rz(-0.0080527877) q[1];
sx q[1];
rz(1.7013928) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9337517) q[0];
sx q[0];
rz(-1.7329204) q[0];
sx q[0];
rz(-0.16926814) q[0];
rz(-pi) q[1];
rz(2.9199706) q[2];
sx q[2];
rz(-1.3563507) q[2];
sx q[2];
rz(-0.947244) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78919856) q[1];
sx q[1];
rz(-1.4161613) q[1];
sx q[1];
rz(1.5782389) q[1];
rz(-pi) q[2];
rz(-0.099248107) q[3];
sx q[3];
rz(-1.7407047) q[3];
sx q[3];
rz(-2.3836294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77233934) q[2];
sx q[2];
rz(-0.99268308) q[2];
sx q[2];
rz(-0.079027979) q[2];
rz(0.8655656) q[3];
sx q[3];
rz(-0.95439684) q[3];
sx q[3];
rz(2.1178093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9158151) q[0];
sx q[0];
rz(-3.1362035) q[0];
sx q[0];
rz(2.916577) q[0];
rz(0.30613884) q[1];
sx q[1];
rz(-3.1248416) q[1];
sx q[1];
rz(0.87047815) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8057833) q[0];
sx q[0];
rz(-3.0558944) q[0];
sx q[0];
rz(0.58914574) q[0];
rz(-0.6323496) q[2];
sx q[2];
rz(-0.82201695) q[2];
sx q[2];
rz(-2.1123304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6806769) q[1];
sx q[1];
rz(-2.3391238) q[1];
sx q[1];
rz(1.6548619) q[1];
rz(0.41272687) q[3];
sx q[3];
rz(-1.1201292) q[3];
sx q[3];
rz(3.009575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6656782) q[2];
sx q[2];
rz(-0.12683882) q[2];
sx q[2];
rz(1.0352943) q[2];
rz(-2.3546442) q[3];
sx q[3];
rz(-3.0988155) q[3];
sx q[3];
rz(2.8662445) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03054522) q[0];
sx q[0];
rz(-0.011140911) q[0];
sx q[0];
rz(0.025644843) q[0];
rz(-1.8772839) q[1];
sx q[1];
rz(-0.023921078) q[1];
sx q[1];
rz(0.69902507) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2866488) q[0];
sx q[0];
rz(-0.35084769) q[0];
sx q[0];
rz(-0.69706608) q[0];
x q[1];
rz(-0.46866663) q[2];
sx q[2];
rz(-1.4184409) q[2];
sx q[2];
rz(3.0909757) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.069534171) q[1];
sx q[1];
rz(-1.1873108) q[1];
sx q[1];
rz(1.6504565) q[1];
x q[2];
rz(2.5949536) q[3];
sx q[3];
rz(-2.1501503) q[3];
sx q[3];
rz(-2.1648615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.25253025) q[2];
sx q[2];
rz(-2.3928596) q[2];
sx q[2];
rz(0.32386455) q[2];
rz(-0.1570541) q[3];
sx q[3];
rz(-1.2374977) q[3];
sx q[3];
rz(2.9878555) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5667628) q[0];
sx q[0];
rz(-3.1269508) q[0];
sx q[0];
rz(0.59004849) q[0];
rz(-0.7198965) q[1];
sx q[1];
rz(-3.0820334) q[1];
sx q[1];
rz(-2.2743026) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.86751) q[0];
sx q[0];
rz(-1.6757586) q[0];
sx q[0];
rz(0.63469736) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9407104) q[2];
sx q[2];
rz(-0.6354161) q[2];
sx q[2];
rz(-0.70178343) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1220951) q[1];
sx q[1];
rz(-1.6830793) q[1];
sx q[1];
rz(-0.068536802) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9525824) q[3];
sx q[3];
rz(-2.2750686) q[3];
sx q[3];
rz(-0.49956043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1934293) q[2];
sx q[2];
rz(-0.89274222) q[2];
sx q[2];
rz(0.11290045) q[2];
rz(-1.0490949) q[3];
sx q[3];
rz(-1.5821404) q[3];
sx q[3];
rz(0.73672867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0650487) q[0];
sx q[0];
rz(-1.5815409) q[0];
sx q[0];
rz(-1.6381868) q[0];
rz(-0.63649559) q[1];
sx q[1];
rz(-0.46000767) q[1];
sx q[1];
rz(1.5719302) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8419452) q[0];
sx q[0];
rz(-1.6498066) q[0];
sx q[0];
rz(1.6539198) q[0];
rz(-pi) q[1];
rz(-0.0031329069) q[2];
sx q[2];
rz(-1.571448) q[2];
sx q[2];
rz(0.80564431) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1366982) q[1];
sx q[1];
rz(-1.5678798) q[1];
sx q[1];
rz(1.5705622) q[1];
rz(0.28255372) q[3];
sx q[3];
rz(-0.5324978) q[3];
sx q[3];
rz(3.0498216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.87127176) q[2];
sx q[2];
rz(-1.511829) q[2];
sx q[2];
rz(2.9809269) q[2];
rz(1.9130982) q[3];
sx q[3];
rz(-0.03385032) q[3];
sx q[3];
rz(-0.20139774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0495618) q[0];
sx q[0];
rz(-1.1783538) q[0];
sx q[0];
rz(0.023068064) q[0];
rz(1.5184215) q[1];
sx q[1];
rz(-2.7802614) q[1];
sx q[1];
rz(-2.8526715) q[1];
rz(1.6466181) q[2];
sx q[2];
rz(-2.8479912) q[2];
sx q[2];
rz(-2.8206024) q[2];
rz(-1.8801943) q[3];
sx q[3];
rz(-2.6831476) q[3];
sx q[3];
rz(0.89098709) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
