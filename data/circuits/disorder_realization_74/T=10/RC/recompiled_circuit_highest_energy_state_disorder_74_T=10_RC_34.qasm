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
rz(1.9407152) q[0];
sx q[0];
rz(-2.241029) q[0];
sx q[0];
rz(-2.9086034) q[0];
rz(-1.5268582) q[1];
sx q[1];
rz(-2.0695217) q[1];
sx q[1];
rz(1.1195247) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4716063) q[0];
sx q[0];
rz(-1.8052088) q[0];
sx q[0];
rz(-3.1283911) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1231261) q[2];
sx q[2];
rz(-0.56150836) q[2];
sx q[2];
rz(0.36040053) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0927375) q[1];
sx q[1];
rz(-1.6462781) q[1];
sx q[1];
rz(-1.2176179) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9078211) q[3];
sx q[3];
rz(-0.98857461) q[3];
sx q[3];
rz(-1.9523417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.69922525) q[2];
sx q[2];
rz(-1.0590326) q[2];
sx q[2];
rz(-0.11288682) q[2];
rz(2.8804307) q[3];
sx q[3];
rz(-1.3449202) q[3];
sx q[3];
rz(-1.2507218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3498822) q[0];
sx q[0];
rz(-0.078502027) q[0];
sx q[0];
rz(-3.0872524) q[0];
rz(-2.9229274) q[1];
sx q[1];
rz(-1.668914) q[1];
sx q[1];
rz(-0.36453077) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35848356) q[0];
sx q[0];
rz(-2.3158584) q[0];
sx q[0];
rz(-2.4381105) q[0];
rz(-pi) q[1];
rz(-2.7414315) q[2];
sx q[2];
rz(-2.1605943) q[2];
sx q[2];
rz(-1.9831744) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5082701) q[1];
sx q[1];
rz(-0.62707096) q[1];
sx q[1];
rz(-0.12992628) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.909799) q[3];
sx q[3];
rz(-1.399013) q[3];
sx q[3];
rz(1.671227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.24078044) q[2];
sx q[2];
rz(-1.4419) q[2];
sx q[2];
rz(1.1085294) q[2];
rz(2.945914) q[3];
sx q[3];
rz(-1.207573) q[3];
sx q[3];
rz(1.4669363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7473258) q[0];
sx q[0];
rz(-2.1013923) q[0];
sx q[0];
rz(-2.105383) q[0];
rz(-0.58468435) q[1];
sx q[1];
rz(-1.4671289) q[1];
sx q[1];
rz(2.5856957) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5824048) q[0];
sx q[0];
rz(-0.83600658) q[0];
sx q[0];
rz(0.57484267) q[0];
rz(-pi) q[1];
rz(2.0228407) q[2];
sx q[2];
rz(-3.0556745) q[2];
sx q[2];
rz(0.8762067) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30196721) q[1];
sx q[1];
rz(-2.4263067) q[1];
sx q[1];
rz(1.7710745) q[1];
x q[2];
rz(2.6922052) q[3];
sx q[3];
rz(-1.3419749) q[3];
sx q[3];
rz(2.902361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3802152) q[2];
sx q[2];
rz(-0.20647241) q[2];
sx q[2];
rz(-1.9492487) q[2];
rz(0.10382593) q[3];
sx q[3];
rz(-1.1622279) q[3];
sx q[3];
rz(-1.9942358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035606774) q[0];
sx q[0];
rz(-2.5150531) q[0];
sx q[0];
rz(-1.4242127) q[0];
rz(2.4183938) q[1];
sx q[1];
rz(-2.9056748) q[1];
sx q[1];
rz(2.9211488) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.854326) q[0];
sx q[0];
rz(-1.1136855) q[0];
sx q[0];
rz(-0.4099538) q[0];
rz(-pi) q[1];
rz(-1.1197243) q[2];
sx q[2];
rz(-0.60904087) q[2];
sx q[2];
rz(-0.90897564) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7539636) q[1];
sx q[1];
rz(-2.7668608) q[1];
sx q[1];
rz(-1.8829569) q[1];
rz(-pi) q[2];
rz(-1.1573767) q[3];
sx q[3];
rz(-1.2261021) q[3];
sx q[3];
rz(2.7935296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7580938) q[2];
sx q[2];
rz(-1.8269962) q[2];
sx q[2];
rz(-2.626075) q[2];
rz(-0.085144194) q[3];
sx q[3];
rz(-0.65535039) q[3];
sx q[3];
rz(2.3084124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67007095) q[0];
sx q[0];
rz(-2.9586198) q[0];
sx q[0];
rz(-2.4772189) q[0];
rz(0.5270671) q[1];
sx q[1];
rz(-1.138569) q[1];
sx q[1];
rz(1.9648431) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7648286) q[0];
sx q[0];
rz(-1.6559147) q[0];
sx q[0];
rz(3.1133988) q[0];
x q[1];
rz(-1.0254622) q[2];
sx q[2];
rz(-2.3787913) q[2];
sx q[2];
rz(-0.57833407) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0461573) q[1];
sx q[1];
rz(-1.5901788) q[1];
sx q[1];
rz(0.23829112) q[1];
rz(-pi) q[2];
x q[2];
rz(2.142435) q[3];
sx q[3];
rz(-1.7222341) q[3];
sx q[3];
rz(-0.12584036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.027017) q[2];
sx q[2];
rz(-0.21275529) q[2];
sx q[2];
rz(-0.93811718) q[2];
rz(2.0959057) q[3];
sx q[3];
rz(-2.9595879) q[3];
sx q[3];
rz(2.3312881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3351347) q[0];
sx q[0];
rz(-1.9430176) q[0];
sx q[0];
rz(-1.2499811) q[0];
rz(-0.20420034) q[1];
sx q[1];
rz(-0.67664346) q[1];
sx q[1];
rz(2.2883889) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6224504) q[0];
sx q[0];
rz(-1.2850437) q[0];
sx q[0];
rz(-0.35730548) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61978842) q[2];
sx q[2];
rz(-2.0896122) q[2];
sx q[2];
rz(-2.3960115) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0895929) q[1];
sx q[1];
rz(-1.1530563) q[1];
sx q[1];
rz(2.3485648) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9167658) q[3];
sx q[3];
rz(-0.36133535) q[3];
sx q[3];
rz(2.0193577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.044067232) q[2];
sx q[2];
rz(-1.7996457) q[2];
sx q[2];
rz(-1.9419144) q[2];
rz(2.1357644) q[3];
sx q[3];
rz(-0.89322105) q[3];
sx q[3];
rz(-0.83542663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47578874) q[0];
sx q[0];
rz(-0.26837334) q[0];
sx q[0];
rz(0.98261181) q[0];
rz(-1.3081374) q[1];
sx q[1];
rz(-1.9255226) q[1];
sx q[1];
rz(1.7024202) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7076232) q[0];
sx q[0];
rz(-2.6585007) q[0];
sx q[0];
rz(0.77204319) q[0];
rz(-pi) q[1];
rz(-1.7050883) q[2];
sx q[2];
rz(-2.6156313) q[2];
sx q[2];
rz(-0.44912042) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8285259) q[1];
sx q[1];
rz(-1.8983574) q[1];
sx q[1];
rz(0.86159535) q[1];
rz(-2.6643326) q[3];
sx q[3];
rz(-1.4443732) q[3];
sx q[3];
rz(1.8074769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16284379) q[2];
sx q[2];
rz(-0.89357251) q[2];
sx q[2];
rz(0.54235512) q[2];
rz(0.98313037) q[3];
sx q[3];
rz(-1.1636846) q[3];
sx q[3];
rz(-2.2014309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5777609) q[0];
sx q[0];
rz(-1.4262119) q[0];
sx q[0];
rz(0.78052178) q[0];
rz(-1.966194) q[1];
sx q[1];
rz(-0.54882097) q[1];
sx q[1];
rz(2.1072047) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8261125) q[0];
sx q[0];
rz(-1.2010935) q[0];
sx q[0];
rz(-1.3913054) q[0];
x q[1];
rz(2.0883191) q[2];
sx q[2];
rz(-0.1642326) q[2];
sx q[2];
rz(-2.9688778) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3745034) q[1];
sx q[1];
rz(-1.4562291) q[1];
sx q[1];
rz(2.7388045) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7195549) q[3];
sx q[3];
rz(-0.58941089) q[3];
sx q[3];
rz(2.5333719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9913651) q[2];
sx q[2];
rz(-1.4069724) q[2];
sx q[2];
rz(0.11037174) q[2];
rz(1.2982093) q[3];
sx q[3];
rz(-1.7896264) q[3];
sx q[3];
rz(2.8492294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2503535) q[0];
sx q[0];
rz(-2.1284916) q[0];
sx q[0];
rz(-1.2497586) q[0];
rz(-0.091014422) q[1];
sx q[1];
rz(-2.4200771) q[1];
sx q[1];
rz(-0.7116085) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7859416) q[0];
sx q[0];
rz(-2.3632746) q[0];
sx q[0];
rz(-3.1319736) q[0];
x q[1];
rz(-2.3046012) q[2];
sx q[2];
rz(-0.92304936) q[2];
sx q[2];
rz(2.456712) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8597128) q[1];
sx q[1];
rz(-1.2285474) q[1];
sx q[1];
rz(0.48515353) q[1];
rz(-pi) q[2];
rz(-1.2290088) q[3];
sx q[3];
rz(-2.2967489) q[3];
sx q[3];
rz(-2.1330698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8083501) q[2];
sx q[2];
rz(-1.658172) q[2];
sx q[2];
rz(-1.9688155) q[2];
rz(1.5277537) q[3];
sx q[3];
rz(-2.0634191) q[3];
sx q[3];
rz(3.0465904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.461819) q[0];
sx q[0];
rz(-0.12400308) q[0];
sx q[0];
rz(-1.5661731) q[0];
rz(-0.35789403) q[1];
sx q[1];
rz(-2.0707668) q[1];
sx q[1];
rz(2.2391052) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6704412) q[0];
sx q[0];
rz(-0.67865279) q[0];
sx q[0];
rz(1.9891692) q[0];
rz(-pi) q[1];
rz(-0.64401354) q[2];
sx q[2];
rz(-1.5531544) q[2];
sx q[2];
rz(-2.0669075) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.984362) q[1];
sx q[1];
rz(-0.8547201) q[1];
sx q[1];
rz(2.6916719) q[1];
rz(1.020458) q[3];
sx q[3];
rz(-1.4963048) q[3];
sx q[3];
rz(2.1410774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0365399) q[2];
sx q[2];
rz(-0.5883216) q[2];
sx q[2];
rz(-1.499929) q[2];
rz(-2.6203652) q[3];
sx q[3];
rz(-0.14385496) q[3];
sx q[3];
rz(0.95411333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4317516) q[0];
sx q[0];
rz(-0.7702282) q[0];
sx q[0];
rz(-1.7256398) q[0];
rz(0.00090986666) q[1];
sx q[1];
rz(-1.4716499) q[1];
sx q[1];
rz(1.6843527) q[1];
rz(-1.7224689) q[2];
sx q[2];
rz(-2.5628452) q[2];
sx q[2];
rz(2.645523) q[2];
rz(-0.72003638) q[3];
sx q[3];
rz(-2.2514718) q[3];
sx q[3];
rz(2.9801647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
