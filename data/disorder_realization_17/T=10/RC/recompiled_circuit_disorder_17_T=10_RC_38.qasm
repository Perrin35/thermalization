OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2330115) q[0];
sx q[0];
rz(-1.1865948) q[0];
sx q[0];
rz(-1.9957805) q[0];
rz(-2.1687578) q[1];
sx q[1];
rz(-1.6701148) q[1];
sx q[1];
rz(2.8491128) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7115294) q[0];
sx q[0];
rz(-1.9646514) q[0];
sx q[0];
rz(-3.0815365) q[0];
rz(-pi) q[1];
rz(-0.77941676) q[2];
sx q[2];
rz(-1.5343108) q[2];
sx q[2];
rz(1.905029) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9521192) q[1];
sx q[1];
rz(-1.6166302) q[1];
sx q[1];
rz(-0.97656753) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1396477) q[3];
sx q[3];
rz(-2.0446152) q[3];
sx q[3];
rz(-0.67809826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4108489) q[2];
sx q[2];
rz(-2.0480672) q[2];
sx q[2];
rz(-2.6921819) q[2];
rz(-0.64569008) q[3];
sx q[3];
rz(-2.460545) q[3];
sx q[3];
rz(-1.2908363) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21355024) q[0];
sx q[0];
rz(-0.88590652) q[0];
sx q[0];
rz(-0.74203062) q[0];
rz(1.4713326) q[1];
sx q[1];
rz(-2.4447618) q[1];
sx q[1];
rz(-1.7785633) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5308967) q[0];
sx q[0];
rz(-2.1801821) q[0];
sx q[0];
rz(-2.0531528) q[0];
rz(1.3554553) q[2];
sx q[2];
rz(-2.0305579) q[2];
sx q[2];
rz(-0.41972566) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.14582536) q[1];
sx q[1];
rz(-1.2312504) q[1];
sx q[1];
rz(1.5911566) q[1];
x q[2];
rz(-2.9491049) q[3];
sx q[3];
rz(-1.3093595) q[3];
sx q[3];
rz(-2.1646433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8403975) q[2];
sx q[2];
rz(-2.6448554) q[2];
sx q[2];
rz(-1.2724686) q[2];
rz(-0.33723351) q[3];
sx q[3];
rz(-1.1703346) q[3];
sx q[3];
rz(-0.01005323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067588016) q[0];
sx q[0];
rz(-2.8019866) q[0];
sx q[0];
rz(-0.2628251) q[0];
rz(2.7858531) q[1];
sx q[1];
rz(-1.6825312) q[1];
sx q[1];
rz(-1.9062818) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4609769) q[0];
sx q[0];
rz(-1.5605698) q[0];
sx q[0];
rz(0.015120487) q[0];
x q[1];
rz(-1.8112196) q[2];
sx q[2];
rz(-1.5449107) q[2];
sx q[2];
rz(-0.38088574) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0318109) q[1];
sx q[1];
rz(-0.33966741) q[1];
sx q[1];
rz(2.8207645) q[1];
x q[2];
rz(0.87576207) q[3];
sx q[3];
rz(-1.2998298) q[3];
sx q[3];
rz(0.7489487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1245023) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(-0.18033218) q[2];
rz(-0.63594937) q[3];
sx q[3];
rz(-1.9790117) q[3];
sx q[3];
rz(-2.7699871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76698774) q[0];
sx q[0];
rz(-0.42414442) q[0];
sx q[0];
rz(2.4225127) q[0];
rz(-2.6285697) q[1];
sx q[1];
rz(-1.2027556) q[1];
sx q[1];
rz(1.0346574) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7906404) q[0];
sx q[0];
rz(-1.6296903) q[0];
sx q[0];
rz(-1.5509997) q[0];
rz(-pi) q[1];
rz(-1.7059533) q[2];
sx q[2];
rz(-1.5237336) q[2];
sx q[2];
rz(-2.1152209) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0532916) q[1];
sx q[1];
rz(-2.3846855) q[1];
sx q[1];
rz(-2.7234368) q[1];
x q[2];
rz(1.3961117) q[3];
sx q[3];
rz(-2.0377199) q[3];
sx q[3];
rz(0.13388453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.93306142) q[2];
sx q[2];
rz(-0.26991093) q[2];
sx q[2];
rz(-0.2362403) q[2];
rz(1.158372) q[3];
sx q[3];
rz(-2.1393496) q[3];
sx q[3];
rz(-0.85197824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029595705) q[0];
sx q[0];
rz(-2.1313666) q[0];
sx q[0];
rz(0.90233666) q[0];
rz(-0.92102712) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(-0.62228084) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3233933) q[0];
sx q[0];
rz(-2.8022031) q[0];
sx q[0];
rz(-2.1470977) q[0];
rz(3.0201966) q[2];
sx q[2];
rz(-2.4897794) q[2];
sx q[2];
rz(-2.7898942) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2942549) q[1];
sx q[1];
rz(-1.6604742) q[1];
sx q[1];
rz(1.7404106) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6617674) q[3];
sx q[3];
rz(-1.2180426) q[3];
sx q[3];
rz(1.477581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92419147) q[2];
sx q[2];
rz(-1.8175586) q[2];
sx q[2];
rz(0.021281555) q[2];
rz(-1.4422669) q[3];
sx q[3];
rz(-0.71443021) q[3];
sx q[3];
rz(0.91059476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4845881) q[0];
sx q[0];
rz(-2.5452884) q[0];
sx q[0];
rz(0.044629991) q[0];
rz(2.1394829) q[1];
sx q[1];
rz(-1.0005181) q[1];
sx q[1];
rz(0.034428509) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7444367) q[0];
sx q[0];
rz(-1.7560871) q[0];
sx q[0];
rz(0.013884355) q[0];
x q[1];
rz(-0.12182932) q[2];
sx q[2];
rz(-1.0708772) q[2];
sx q[2];
rz(3.1140529) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.62528246) q[1];
sx q[1];
rz(-2.1045661) q[1];
sx q[1];
rz(2.7620402) q[1];
x q[2];
rz(0.19272007) q[3];
sx q[3];
rz(-2.4270714) q[3];
sx q[3];
rz(-0.89380985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3520711) q[2];
sx q[2];
rz(-1.2723203) q[2];
sx q[2];
rz(-2.1748523) q[2];
rz(-0.012185193) q[3];
sx q[3];
rz(-0.92958486) q[3];
sx q[3];
rz(-3.0325082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0063342) q[0];
sx q[0];
rz(-0.2671347) q[0];
sx q[0];
rz(-2.1475041) q[0];
rz(0.055123568) q[1];
sx q[1];
rz(-1.2722641) q[1];
sx q[1];
rz(-0.73928839) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71541876) q[0];
sx q[0];
rz(-2.9746375) q[0];
sx q[0];
rz(-0.963361) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0544858) q[2];
sx q[2];
rz(-2.0927883) q[2];
sx q[2];
rz(-2.1453478) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2412487) q[1];
sx q[1];
rz(-2.5785891) q[1];
sx q[1];
rz(2.7808933) q[1];
rz(0.39412734) q[3];
sx q[3];
rz(-2.4836802) q[3];
sx q[3];
rz(3.1042299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58549515) q[2];
sx q[2];
rz(-1.408351) q[2];
sx q[2];
rz(-0.56376702) q[2];
rz(0.24027763) q[3];
sx q[3];
rz(-2.6085745) q[3];
sx q[3];
rz(1.7668004) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4450842) q[0];
sx q[0];
rz(-2.962528) q[0];
sx q[0];
rz(-2.5575496) q[0];
rz(-0.42075992) q[1];
sx q[1];
rz(-1.0412443) q[1];
sx q[1];
rz(0.27841321) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1906492) q[0];
sx q[0];
rz(-1.6961432) q[0];
sx q[0];
rz(2.1423253) q[0];
rz(-pi) q[1];
rz(-2.7114696) q[2];
sx q[2];
rz(-2.7134036) q[2];
sx q[2];
rz(0.45454121) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1380996) q[1];
sx q[1];
rz(-1.3521191) q[1];
sx q[1];
rz(1.2117282) q[1];
rz(-1.0233381) q[3];
sx q[3];
rz(-1.3609144) q[3];
sx q[3];
rz(1.4048214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.86928308) q[2];
sx q[2];
rz(-2.3217106) q[2];
sx q[2];
rz(2.753624) q[2];
rz(-1.7913473) q[3];
sx q[3];
rz(-2.6781121) q[3];
sx q[3];
rz(-2.8665682) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85912722) q[0];
sx q[0];
rz(-2.9282741) q[0];
sx q[0];
rz(2.4097089) q[0];
rz(-2.9751119) q[1];
sx q[1];
rz(-0.44547588) q[1];
sx q[1];
rz(-3.0678715) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6442935) q[0];
sx q[0];
rz(-1.4696346) q[0];
sx q[0];
rz(1.1343207) q[0];
rz(-0.080967112) q[2];
sx q[2];
rz(-1.4514187) q[2];
sx q[2];
rz(-2.2031684) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8158055) q[1];
sx q[1];
rz(-0.98762883) q[1];
sx q[1];
rz(-1.5259685) q[1];
rz(1.681626) q[3];
sx q[3];
rz(-1.0004527) q[3];
sx q[3];
rz(-0.092872083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59447294) q[2];
sx q[2];
rz(-0.23590817) q[2];
sx q[2];
rz(0.42868844) q[2];
rz(1.3171014) q[3];
sx q[3];
rz(-1.3573815) q[3];
sx q[3];
rz(-1.930442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73228943) q[0];
sx q[0];
rz(-1.6560873) q[0];
sx q[0];
rz(2.0987341) q[0];
rz(1.6304784) q[1];
sx q[1];
rz(-1.7621367) q[1];
sx q[1];
rz(1.3483378) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3139075) q[0];
sx q[0];
rz(-1.5450918) q[0];
sx q[0];
rz(-1.4575507) q[0];
rz(-pi) q[1];
rz(0.83421631) q[2];
sx q[2];
rz(-1.5110821) q[2];
sx q[2];
rz(2.5028253) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9650324) q[1];
sx q[1];
rz(-1.5000048) q[1];
sx q[1];
rz(-1.2464347) q[1];
rz(-pi) q[2];
rz(-1.7616368) q[3];
sx q[3];
rz(-1.3340501) q[3];
sx q[3];
rz(2.5773347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3412791) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(0.094816118) q[2];
rz(-1.9684277) q[3];
sx q[3];
rz(-1.7404107) q[3];
sx q[3];
rz(-2.9828984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6288347) q[0];
sx q[0];
rz(-2.6678968) q[0];
sx q[0];
rz(-2.0976023) q[0];
rz(-1.6013153) q[1];
sx q[1];
rz(-1.5581144) q[1];
sx q[1];
rz(1.5099572) q[1];
rz(-1.9376471) q[2];
sx q[2];
rz(-2.9137712) q[2];
sx q[2];
rz(-3.055345) q[2];
rz(0.81090609) q[3];
sx q[3];
rz(-1.1699642) q[3];
sx q[3];
rz(-1.5311833) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
