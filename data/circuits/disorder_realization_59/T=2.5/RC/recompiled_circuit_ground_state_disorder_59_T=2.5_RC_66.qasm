OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0629118) q[0];
sx q[0];
rz(-2.7432888) q[0];
sx q[0];
rz(-2.8226573) q[0];
rz(2.6755264) q[1];
sx q[1];
rz(-1.4332486) q[1];
sx q[1];
rz(2.8679009) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2994933) q[0];
sx q[0];
rz(-2.2998739) q[0];
sx q[0];
rz(1.0925348) q[0];
rz(-pi) q[1];
rz(-1.798965) q[2];
sx q[2];
rz(-1.3334477) q[2];
sx q[2];
rz(-0.37712032) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92299709) q[1];
sx q[1];
rz(-1.4540919) q[1];
sx q[1];
rz(-1.3918819) q[1];
rz(-pi) q[2];
rz(2.6998547) q[3];
sx q[3];
rz(-2.3328247) q[3];
sx q[3];
rz(-1.3490775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29229257) q[2];
sx q[2];
rz(-0.92508525) q[2];
sx q[2];
rz(1.383847) q[2];
rz(-0.84550953) q[3];
sx q[3];
rz(-0.011761646) q[3];
sx q[3];
rz(-2.1498634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.200835) q[0];
sx q[0];
rz(-2.2469914) q[0];
sx q[0];
rz(-2.1064598) q[0];
rz(-1.9738522) q[1];
sx q[1];
rz(-2.9393241) q[1];
sx q[1];
rz(-2.367173) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5308038) q[0];
sx q[0];
rz(-2.4770081) q[0];
sx q[0];
rz(-2.2749645) q[0];
x q[1];
rz(-2.7970052) q[2];
sx q[2];
rz(-2.2026052) q[2];
sx q[2];
rz(0.58877221) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8138511) q[1];
sx q[1];
rz(-2.6151513) q[1];
sx q[1];
rz(-1.8000283) q[1];
rz(-0.033942698) q[3];
sx q[3];
rz(-1.6210302) q[3];
sx q[3];
rz(1.4242581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8059798) q[2];
sx q[2];
rz(-0.077294417) q[2];
sx q[2];
rz(0.94266164) q[2];
rz(-0.089426905) q[3];
sx q[3];
rz(-2.4783897) q[3];
sx q[3];
rz(-2.8360143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59219229) q[0];
sx q[0];
rz(-2.5553199) q[0];
sx q[0];
rz(0.19202448) q[0];
rz(-2.5281455) q[1];
sx q[1];
rz(-2.6934721) q[1];
sx q[1];
rz(3.1268069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3131994) q[0];
sx q[0];
rz(-1.9678741) q[0];
sx q[0];
rz(2.3262789) q[0];
rz(-1.2959145) q[2];
sx q[2];
rz(-2.2408591) q[2];
sx q[2];
rz(2.2704147) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.80577981) q[1];
sx q[1];
rz(-1.3381914) q[1];
sx q[1];
rz(2.1628153) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5824541) q[3];
sx q[3];
rz(-1.634667) q[3];
sx q[3];
rz(3.0272878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6141367) q[2];
sx q[2];
rz(-1.5469896) q[2];
sx q[2];
rz(0.0097489348) q[2];
rz(-2.5629432) q[3];
sx q[3];
rz(-1.0520244) q[3];
sx q[3];
rz(-0.78154045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11827271) q[0];
sx q[0];
rz(-0.8096205) q[0];
sx q[0];
rz(-0.6672346) q[0];
rz(-1.0046593) q[1];
sx q[1];
rz(-0.15548448) q[1];
sx q[1];
rz(-1.3803233) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8341433) q[0];
sx q[0];
rz(-0.47795313) q[0];
sx q[0];
rz(1.5909999) q[0];
rz(2.0441846) q[2];
sx q[2];
rz(-1.2815099) q[2];
sx q[2];
rz(2.8662967) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9050153) q[1];
sx q[1];
rz(-1.5386174) q[1];
sx q[1];
rz(0.49361146) q[1];
x q[2];
rz(-1.0466688) q[3];
sx q[3];
rz(-1.7780684) q[3];
sx q[3];
rz(2.2818499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7669547) q[2];
sx q[2];
rz(-2.0071603) q[2];
sx q[2];
rz(-3.0881506) q[2];
rz(0.63353574) q[3];
sx q[3];
rz(-0.41826785) q[3];
sx q[3];
rz(-3.1327278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0165789) q[0];
sx q[0];
rz(-1.582394) q[0];
sx q[0];
rz(-1.5484823) q[0];
rz(-pi) q[1];
rz(-2.1778278) q[2];
sx q[2];
rz(-0.53812611) q[2];
sx q[2];
rz(-0.48689688) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5875467) q[1];
sx q[1];
rz(-0.57922208) q[1];
sx q[1];
rz(-0.86279561) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76988585) q[3];
sx q[3];
rz(-1.5428233) q[3];
sx q[3];
rz(-2.6370492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2962467) q[2];
sx q[2];
rz(-2.3643934) q[2];
sx q[2];
rz(-3.0261611) q[2];
rz(-0.7655862) q[3];
sx q[3];
rz(-1.4976394) q[3];
sx q[3];
rz(2.7287741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.63824832) q[0];
sx q[0];
rz(-2.4920576) q[0];
sx q[0];
rz(0.58393884) q[0];
rz(0.04843796) q[1];
sx q[1];
rz(-2.9190639) q[1];
sx q[1];
rz(2.4979874) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2023425) q[0];
sx q[0];
rz(-1.5233375) q[0];
sx q[0];
rz(-1.2623293) q[0];
rz(-pi) q[1];
rz(0.46508772) q[2];
sx q[2];
rz(-2.1784867) q[2];
sx q[2];
rz(-1.6179784) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.19932718) q[1];
sx q[1];
rz(-0.86546997) q[1];
sx q[1];
rz(2.9678515) q[1];
x q[2];
rz(-2.3619283) q[3];
sx q[3];
rz(-1.4309959) q[3];
sx q[3];
rz(-1.1598905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.863997) q[2];
sx q[2];
rz(-2.3433351) q[2];
sx q[2];
rz(0.64211988) q[2];
rz(3.0025205) q[3];
sx q[3];
rz(-0.13834794) q[3];
sx q[3];
rz(1.4596435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7794401) q[0];
sx q[0];
rz(-0.42069778) q[0];
sx q[0];
rz(-2.5162589) q[0];
rz(0.23896898) q[1];
sx q[1];
rz(-0.23663722) q[1];
sx q[1];
rz(-1.3561358) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55955333) q[0];
sx q[0];
rz(-1.516321) q[0];
sx q[0];
rz(0.031886851) q[0];
rz(-0.37650336) q[2];
sx q[2];
rz(-1.8424705) q[2];
sx q[2];
rz(1.5282549) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5290288) q[1];
sx q[1];
rz(-0.41464889) q[1];
sx q[1];
rz(-2.0842444) q[1];
rz(-1.4144355) q[3];
sx q[3];
rz(-1.7810453) q[3];
sx q[3];
rz(2.4429136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8465392) q[2];
sx q[2];
rz(-1.2650547) q[2];
sx q[2];
rz(2.1758101) q[2];
rz(0.24671181) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63916373) q[0];
sx q[0];
rz(-2.736709) q[0];
sx q[0];
rz(-3.0054481) q[0];
rz(2.4038521) q[1];
sx q[1];
rz(-2.9164011) q[1];
sx q[1];
rz(-1.236261) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.284974) q[0];
sx q[0];
rz(-2.1559733) q[0];
sx q[0];
rz(2.7057458) q[0];
rz(-pi) q[1];
rz(-2.9604016) q[2];
sx q[2];
rz(-0.94215067) q[2];
sx q[2];
rz(-1.218007) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.41676109) q[1];
sx q[1];
rz(-0.5151075) q[1];
sx q[1];
rz(-0.58578844) q[1];
rz(0.32847877) q[3];
sx q[3];
rz(-0.94624472) q[3];
sx q[3];
rz(0.49027944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.69011921) q[2];
sx q[2];
rz(-1.5560919) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5600679) q[0];
sx q[0];
rz(-0.6811322) q[0];
sx q[0];
rz(-3.0233622) q[0];
rz(2.0589578) q[1];
sx q[1];
rz(-1.9397468) q[1];
sx q[1];
rz(-1.3523678) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0467606) q[0];
sx q[0];
rz(-0.40725312) q[0];
sx q[0];
rz(-0.26131765) q[0];
x q[1];
rz(-0.1583516) q[2];
sx q[2];
rz(-2.5917946) q[2];
sx q[2];
rz(-2.9937772) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.73001) q[1];
sx q[1];
rz(-1.4435539) q[1];
sx q[1];
rz(1.0858658) q[1];
rz(-pi) q[2];
rz(-0.31183536) q[3];
sx q[3];
rz(-2.0297628) q[3];
sx q[3];
rz(0.74742693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.81789404) q[2];
sx q[2];
rz(-2.4511621) q[2];
sx q[2];
rz(0.77198088) q[2];
rz(0.082993232) q[3];
sx q[3];
rz(-1.8738184) q[3];
sx q[3];
rz(0.52223372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25987396) q[0];
sx q[0];
rz(-0.75749713) q[0];
sx q[0];
rz(0.74257332) q[0];
rz(-1.2965797) q[1];
sx q[1];
rz(-2.3340338) q[1];
sx q[1];
rz(-0.47320941) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6542235) q[0];
sx q[0];
rz(-0.28074902) q[0];
sx q[0];
rz(-1.8868179) q[0];
x q[1];
rz(2.3397105) q[2];
sx q[2];
rz(-2.7919765) q[2];
sx q[2];
rz(0.59211087) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0124346) q[1];
sx q[1];
rz(-1.9170463) q[1];
sx q[1];
rz(-1.0770256) q[1];
rz(-pi) q[2];
rz(0.74537422) q[3];
sx q[3];
rz(-1.3826118) q[3];
sx q[3];
rz(-2.4309576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.55127281) q[2];
sx q[2];
rz(-1.4164305) q[2];
sx q[2];
rz(-1.2335221) q[2];
rz(-2.7963855) q[3];
sx q[3];
rz(-0.39185169) q[3];
sx q[3];
rz(-0.83693081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45967669) q[0];
sx q[0];
rz(-1.3688594) q[0];
sx q[0];
rz(-1.3874227) q[0];
rz(-0.058902901) q[1];
sx q[1];
rz(-1.420493) q[1];
sx q[1];
rz(-2.4891985) q[1];
rz(-1.1299533) q[2];
sx q[2];
rz(-1.2193573) q[2];
sx q[2];
rz(-1.0680547) q[2];
rz(2.160499) q[3];
sx q[3];
rz(-0.72219875) q[3];
sx q[3];
rz(1.2027702) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
