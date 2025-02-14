OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0786809) q[0];
sx q[0];
rz(-0.39830387) q[0];
sx q[0];
rz(-0.31893536) q[0];
rz(2.6755264) q[1];
sx q[1];
rz(-1.4332486) q[1];
sx q[1];
rz(-0.27369174) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39616991) q[0];
sx q[0];
rz(-1.92116) q[0];
sx q[0];
rz(-0.78846447) q[0];
x q[1];
rz(-0.24342033) q[2];
sx q[2];
rz(-1.3491329) q[2];
sx q[2];
rz(1.8933715) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4727386) q[1];
sx q[1];
rz(-1.3931119) q[1];
sx q[1];
rz(-0.11858003) q[1];
rz(-pi) q[2];
rz(2.3832604) q[3];
sx q[3];
rz(-1.8852295) q[3];
sx q[3];
rz(3.0477332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29229257) q[2];
sx q[2];
rz(-0.92508525) q[2];
sx q[2];
rz(-1.383847) q[2];
rz(-0.84550953) q[3];
sx q[3];
rz(-3.129831) q[3];
sx q[3];
rz(-0.99172926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94075769) q[0];
sx q[0];
rz(-2.2469914) q[0];
sx q[0];
rz(2.1064598) q[0];
rz(1.1677405) q[1];
sx q[1];
rz(-0.20226856) q[1];
sx q[1];
rz(-0.77441961) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37062936) q[0];
sx q[0];
rz(-1.1600732) q[0];
sx q[0];
rz(2.109101) q[0];
rz(2.2317288) q[2];
sx q[2];
rz(-1.294704) q[2];
sx q[2];
rz(-0.7731437) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32774156) q[1];
sx q[1];
rz(-0.5264414) q[1];
sx q[1];
rz(1.8000283) q[1];
x q[2];
rz(-1.5205335) q[3];
sx q[3];
rz(-1.5368965) q[3];
sx q[3];
rz(-2.9933495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8059798) q[2];
sx q[2];
rz(-0.077294417) q[2];
sx q[2];
rz(0.94266164) q[2];
rz(-3.0521657) q[3];
sx q[3];
rz(-0.663203) q[3];
sx q[3];
rz(0.30557835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59219229) q[0];
sx q[0];
rz(-2.5553199) q[0];
sx q[0];
rz(2.9495682) q[0];
rz(2.5281455) q[1];
sx q[1];
rz(-2.6934721) q[1];
sx q[1];
rz(0.014785756) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0917182) q[0];
sx q[0];
rz(-0.8862859) q[0];
sx q[0];
rz(-2.6189463) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33003557) q[2];
sx q[2];
rz(-0.71612172) q[2];
sx q[2];
rz(2.696685) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0952436) q[1];
sx q[1];
rz(-2.5106437) q[1];
sx q[1];
rz(1.1693404) q[1];
x q[2];
rz(-1.5824541) q[3];
sx q[3];
rz(-1.634667) q[3];
sx q[3];
rz(-0.11430489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6141367) q[2];
sx q[2];
rz(-1.5469896) q[2];
sx q[2];
rz(-3.1318437) q[2];
rz(0.57864946) q[3];
sx q[3];
rz(-1.0520244) q[3];
sx q[3];
rz(2.3600522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0233199) q[0];
sx q[0];
rz(-0.8096205) q[0];
sx q[0];
rz(0.6672346) q[0];
rz(1.0046593) q[1];
sx q[1];
rz(-2.9861082) q[1];
sx q[1];
rz(1.7612693) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8568958) q[0];
sx q[0];
rz(-2.0486437) q[0];
sx q[0];
rz(-0.010464593) q[0];
rz(-pi) q[1];
rz(2.1491702) q[2];
sx q[2];
rz(-2.592591) q[2];
sx q[2];
rz(1.3379607) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.31691027) q[1];
sx q[1];
rz(-2.0641293) q[1];
sx q[1];
rz(1.534259) q[1];
rz(-pi) q[2];
rz(-2.9033086) q[3];
sx q[3];
rz(-2.0825948) q[3];
sx q[3];
rz(-0.82945591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37463793) q[2];
sx q[2];
rz(-1.1344323) q[2];
sx q[2];
rz(-3.0881506) q[2];
rz(-0.63353574) q[3];
sx q[3];
rz(-0.41826785) q[3];
sx q[3];
rz(3.1327278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2706914) q[0];
sx q[0];
rz(-2.5636261) q[0];
sx q[0];
rz(-2.0722678) q[0];
rz(-1.2879734) q[1];
sx q[1];
rz(-0.063871495) q[1];
sx q[1];
rz(3.0499444) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44604135) q[0];
sx q[0];
rz(-1.5484838) q[0];
sx q[0];
rz(-3.1299921) q[0];
rz(1.1149801) q[2];
sx q[2];
rz(-1.2741003) q[2];
sx q[2];
rz(-1.5199583) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.39508546) q[1];
sx q[1];
rz(-1.9347435) q[1];
sx q[1];
rz(1.1096611) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.609747) q[3];
sx q[3];
rz(-0.80128968) q[3];
sx q[3];
rz(1.0391446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.84534591) q[2];
sx q[2];
rz(-2.3643934) q[2];
sx q[2];
rz(-0.11543154) q[2];
rz(0.7655862) q[3];
sx q[3];
rz(-1.4976394) q[3];
sx q[3];
rz(0.41281858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5033443) q[0];
sx q[0];
rz(-0.64953506) q[0];
sx q[0];
rz(2.5576538) q[0];
rz(0.04843796) q[1];
sx q[1];
rz(-0.2225288) q[1];
sx q[1];
rz(-2.4979874) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7580235) q[0];
sx q[0];
rz(-1.2626881) q[0];
sx q[0];
rz(3.0917866) q[0];
rz(2.2320643) q[2];
sx q[2];
rz(-1.9478746) q[2];
sx q[2];
rz(-0.23185767) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.065066628) q[1];
sx q[1];
rz(-2.4187633) q[1];
sx q[1];
rz(-1.7711054) q[1];
x q[2];
rz(-1.7661473) q[3];
sx q[3];
rz(-0.80073154) q[3];
sx q[3];
rz(-2.8675818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.863997) q[2];
sx q[2];
rz(-2.3433351) q[2];
sx q[2];
rz(2.4994728) q[2];
rz(-0.13907214) q[3];
sx q[3];
rz(-0.13834794) q[3];
sx q[3];
rz(1.4596435) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794401) q[0];
sx q[0];
rz(-0.42069778) q[0];
sx q[0];
rz(-0.62533373) q[0];
rz(2.9026237) q[1];
sx q[1];
rz(-0.23663722) q[1];
sx q[1];
rz(-1.7854569) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1286129) q[0];
sx q[0];
rz(-1.6026359) q[0];
sx q[0];
rz(1.6252993) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2797592) q[2];
sx q[2];
rz(-1.2087529) q[2];
sx q[2];
rz(2.9933528) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.060387) q[1];
sx q[1];
rz(-1.9293509) q[1];
sx q[1];
rz(2.9286659) q[1];
rz(-pi) q[2];
rz(0.21276833) q[3];
sx q[3];
rz(-1.4179061) q[3];
sx q[3];
rz(-2.2365856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.29505342) q[2];
sx q[2];
rz(-1.876538) q[2];
sx q[2];
rz(-2.1758101) q[2];
rz(-2.8948808) q[3];
sx q[3];
rz(-1.8762981) q[3];
sx q[3];
rz(1.770796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5024289) q[0];
sx q[0];
rz(-2.736709) q[0];
sx q[0];
rz(-3.0054481) q[0];
rz(2.4038521) q[1];
sx q[1];
rz(-2.9164011) q[1];
sx q[1];
rz(-1.236261) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6039918) q[0];
sx q[0];
rz(-1.21116) q[0];
sx q[0];
rz(0.93961544) q[0];
rz(-1.8137553) q[2];
sx q[2];
rz(-2.4907673) q[2];
sx q[2];
rz(-1.6215768) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63041362) q[1];
sx q[1];
rz(-1.8466338) q[1];
sx q[1];
rz(-0.44075573) q[1];
rz(-pi) q[2];
rz(-2.2216283) q[3];
sx q[3];
rz(-1.3060089) q[3];
sx q[3];
rz(-1.277232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.69011921) q[2];
sx q[2];
rz(-1.5855007) q[2];
sx q[2];
rz(2.1915009) q[2];
rz(2.7443366) q[3];
sx q[3];
rz(-0.49644956) q[3];
sx q[3];
rz(2.6955786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5815247) q[0];
sx q[0];
rz(-0.6811322) q[0];
sx q[0];
rz(3.0233622) q[0];
rz(2.0589578) q[1];
sx q[1];
rz(-1.9397468) q[1];
sx q[1];
rz(1.7892249) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7633469) q[0];
sx q[0];
rz(-1.9634501) q[0];
sx q[0];
rz(-1.4598085) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4744583) q[2];
sx q[2];
rz(-1.0286478) q[2];
sx q[2];
rz(3.1042636) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.73001) q[1];
sx q[1];
rz(-1.6980387) q[1];
sx q[1];
rz(2.0557269) q[1];
x q[2];
rz(-2.8297573) q[3];
sx q[3];
rz(-2.0297628) q[3];
sx q[3];
rz(-0.74742693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.81789404) q[2];
sx q[2];
rz(-2.4511621) q[2];
sx q[2];
rz(0.77198088) q[2];
rz(-0.082993232) q[3];
sx q[3];
rz(-1.8738184) q[3];
sx q[3];
rz(-0.52223372) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8817187) q[0];
sx q[0];
rz(-0.75749713) q[0];
sx q[0];
rz(2.3990193) q[0];
rz(-1.2965797) q[1];
sx q[1];
rz(-0.80755889) q[1];
sx q[1];
rz(-2.6683832) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98223393) q[0];
sx q[0];
rz(-1.8372941) q[0];
sx q[0];
rz(3.052211) q[0];
rz(-1.8270566) q[2];
sx q[2];
rz(-1.8112931) q[2];
sx q[2];
rz(0.24090361) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4025583) q[1];
sx q[1];
rz(-1.1087045) q[1];
sx q[1];
rz(-0.38886221) q[1];
x q[2];
rz(2.8678611) q[3];
sx q[3];
rz(-2.377284) q[3];
sx q[3];
rz(2.4814062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.55127281) q[2];
sx q[2];
rz(-1.4164305) q[2];
sx q[2];
rz(1.9080706) q[2];
rz(0.34520712) q[3];
sx q[3];
rz(-2.749741) q[3];
sx q[3];
rz(-2.3046618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.681916) q[0];
sx q[0];
rz(-1.7727333) q[0];
sx q[0];
rz(1.75417) q[0];
rz(-3.0826898) q[1];
sx q[1];
rz(-1.7210996) q[1];
sx q[1];
rz(0.65239418) q[1];
rz(2.0116393) q[2];
sx q[2];
rz(-1.2193573) q[2];
sx q[2];
rz(-1.0680547) q[2];
rz(-0.93880063) q[3];
sx q[3];
rz(-1.1943571) q[3];
sx q[3];
rz(-0.83333701) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
