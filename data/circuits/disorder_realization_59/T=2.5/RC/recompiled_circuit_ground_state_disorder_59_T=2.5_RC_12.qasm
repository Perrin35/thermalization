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
rz(2.8226573) q[0];
rz(-0.4660663) q[1];
sx q[1];
rz(4.5748413) q[1];
sx q[1];
rz(9.6984697) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84209937) q[0];
sx q[0];
rz(-0.84171879) q[0];
sx q[0];
rz(-1.0925348) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3897522) q[2];
sx q[2];
rz(-2.813857) q[2];
sx q[2];
rz(-0.40204266) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4727386) q[1];
sx q[1];
rz(-1.3931119) q[1];
sx q[1];
rz(0.11858003) q[1];
rz(-2.6998547) q[3];
sx q[3];
rz(-0.80876793) q[3];
sx q[3];
rz(1.7925151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29229257) q[2];
sx q[2];
rz(-2.2165074) q[2];
sx q[2];
rz(1.383847) q[2];
rz(2.2960831) q[3];
sx q[3];
rz(-3.129831) q[3];
sx q[3];
rz(2.1498634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94075769) q[0];
sx q[0];
rz(-2.2469914) q[0];
sx q[0];
rz(2.1064598) q[0];
rz(1.9738522) q[1];
sx q[1];
rz(-0.20226856) q[1];
sx q[1];
rz(-2.367173) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4342142) q[0];
sx q[0];
rz(-2.0601354) q[0];
sx q[0];
rz(2.6721832) q[0];
rz(-pi) q[1];
rz(-2.2317288) q[2];
sx q[2];
rz(-1.8468886) q[2];
sx q[2];
rz(-0.7731437) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0774457) q[1];
sx q[1];
rz(-1.059491) q[1];
sx q[1];
rz(-3.0102986) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1645427) q[3];
sx q[3];
rz(-3.0809743) q[3];
sx q[3];
rz(-0.82965899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3356129) q[2];
sx q[2];
rz(-0.077294417) q[2];
sx q[2];
rz(0.94266164) q[2];
rz(0.089426905) q[3];
sx q[3];
rz(-2.4783897) q[3];
sx q[3];
rz(-0.30557835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5494004) q[0];
sx q[0];
rz(-0.58627272) q[0];
sx q[0];
rz(0.19202448) q[0];
rz(-0.61344719) q[1];
sx q[1];
rz(-2.6934721) q[1];
sx q[1];
rz(0.014785756) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7887869) q[0];
sx q[0];
rz(-0.83483046) q[0];
sx q[0];
rz(2.1197181) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8115571) q[2];
sx q[2];
rz(-0.71612172) q[2];
sx q[2];
rz(-2.696685) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0952436) q[1];
sx q[1];
rz(-2.5106437) q[1];
sx q[1];
rz(1.9722523) q[1];
x q[2];
rz(1.5824541) q[3];
sx q[3];
rz(-1.634667) q[3];
sx q[3];
rz(-3.0272878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6141367) q[2];
sx q[2];
rz(-1.5469896) q[2];
sx q[2];
rz(-3.1318437) q[2];
rz(-0.57864946) q[3];
sx q[3];
rz(-2.0895683) q[3];
sx q[3];
rz(-0.78154045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0233199) q[0];
sx q[0];
rz(-0.8096205) q[0];
sx q[0];
rz(-2.4743581) q[0];
rz(1.0046593) q[1];
sx q[1];
rz(-0.15548448) q[1];
sx q[1];
rz(-1.7612693) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8568958) q[0];
sx q[0];
rz(-2.0486437) q[0];
sx q[0];
rz(0.010464593) q[0];
x q[1];
rz(-2.0441846) q[2];
sx q[2];
rz(-1.2815099) q[2];
sx q[2];
rz(-2.8662967) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9050153) q[1];
sx q[1];
rz(-1.6029753) q[1];
sx q[1];
rz(0.49361146) q[1];
rz(-0.23828407) q[3];
sx q[3];
rz(-1.0589979) q[3];
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
rz(2.5080569) q[3];
sx q[3];
rz(-0.41826785) q[3];
sx q[3];
rz(3.1327278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87090129) q[0];
sx q[0];
rz(-2.5636261) q[0];
sx q[0];
rz(-1.0693249) q[0];
rz(-1.8536192) q[1];
sx q[1];
rz(-3.0777212) q[1];
sx q[1];
rz(3.0499444) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44604135) q[0];
sx q[0];
rz(-1.5931089) q[0];
sx q[0];
rz(-3.1299921) q[0];
rz(-pi) q[1];
rz(0.32817082) q[2];
sx q[2];
rz(-2.0053021) q[2];
sx q[2];
rz(2.948394) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5875467) q[1];
sx q[1];
rz(-0.57922208) q[1];
sx q[1];
rz(-2.278797) q[1];
rz(-pi) q[2];
rz(1.5318457) q[3];
sx q[3];
rz(-2.340303) q[3];
sx q[3];
rz(-1.0391446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2962467) q[2];
sx q[2];
rz(-0.77719921) q[2];
sx q[2];
rz(-0.11543154) q[2];
rz(-0.7655862) q[3];
sx q[3];
rz(-1.4976394) q[3];
sx q[3];
rz(-0.41281858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5033443) q[0];
sx q[0];
rz(-0.64953506) q[0];
sx q[0];
rz(-0.58393884) q[0];
rz(-3.0931547) q[1];
sx q[1];
rz(-0.2225288) q[1];
sx q[1];
rz(-2.4979874) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38356919) q[0];
sx q[0];
rz(-1.2626881) q[0];
sx q[0];
rz(3.0917866) q[0];
x q[1];
rz(0.46508772) q[2];
sx q[2];
rz(-2.1784867) q[2];
sx q[2];
rz(1.5236142) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.076526) q[1];
sx q[1];
rz(-2.4187633) q[1];
sx q[1];
rz(1.3704872) q[1];
x q[2];
rz(-2.3619283) q[3];
sx q[3];
rz(-1.4309959) q[3];
sx q[3];
rz(1.9817022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27759564) q[2];
sx q[2];
rz(-0.79825753) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794401) q[0];
sx q[0];
rz(-2.7208949) q[0];
sx q[0];
rz(0.62533373) q[0];
rz(-0.23896898) q[1];
sx q[1];
rz(-0.23663722) q[1];
sx q[1];
rz(-1.7854569) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1286129) q[0];
sx q[0];
rz(-1.5389568) q[0];
sx q[0];
rz(1.6252993) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37650336) q[2];
sx q[2];
rz(-1.2991221) q[2];
sx q[2];
rz(1.6133378) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6125638) q[1];
sx q[1];
rz(-2.7269438) q[1];
sx q[1];
rz(1.0573483) q[1];
x q[2];
rz(2.9288243) q[3];
sx q[3];
rz(-1.7236865) q[3];
sx q[3];
rz(-2.2365856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29505342) q[2];
sx q[2];
rz(-1.876538) q[2];
sx q[2];
rz(-2.1758101) q[2];
rz(-2.8948808) q[3];
sx q[3];
rz(-1.2652946) q[3];
sx q[3];
rz(1.3707967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5024289) q[0];
sx q[0];
rz(-0.40488365) q[0];
sx q[0];
rz(3.0054481) q[0];
rz(-2.4038521) q[1];
sx q[1];
rz(-2.9164011) q[1];
sx q[1];
rz(-1.9053316) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.284974) q[0];
sx q[0];
rz(-2.1559733) q[0];
sx q[0];
rz(2.7057458) q[0];
rz(-pi) q[1];
rz(-2.2073123) q[2];
sx q[2];
rz(-1.4245241) q[2];
sx q[2];
rz(0.24547228) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0734399) q[1];
sx q[1];
rz(-1.9937939) q[1];
sx q[1];
rz(1.8741029) q[1];
x q[2];
rz(2.2216283) q[3];
sx q[3];
rz(-1.3060089) q[3];
sx q[3];
rz(1.277232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4514734) q[2];
sx q[2];
rz(-1.5560919) q[2];
sx q[2];
rz(-0.95009178) q[2];
rz(-2.7443366) q[3];
sx q[3];
rz(-2.6451431) q[3];
sx q[3];
rz(-0.44601405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5815247) q[0];
sx q[0];
rz(-2.4604605) q[0];
sx q[0];
rz(-0.11823046) q[0];
rz(2.0589578) q[1];
sx q[1];
rz(-1.9397468) q[1];
sx q[1];
rz(-1.3523678) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9064241) q[0];
sx q[0];
rz(-1.6733067) q[0];
sx q[0];
rz(-0.39484039) q[0];
x q[1];
rz(2.5973877) q[2];
sx q[2];
rz(-1.4883071) q[2];
sx q[2];
rz(-1.5832886) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7460359) q[1];
sx q[1];
rz(-0.50005752) q[1];
sx q[1];
rz(1.8386503) q[1];
x q[2];
rz(0.31183536) q[3];
sx q[3];
rz(-2.0297628) q[3];
sx q[3];
rz(2.3941657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3236986) q[2];
sx q[2];
rz(-2.4511621) q[2];
sx q[2];
rz(-0.77198088) q[2];
rz(3.0585994) q[3];
sx q[3];
rz(-1.8738184) q[3];
sx q[3];
rz(-0.52223372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8817187) q[0];
sx q[0];
rz(-0.75749713) q[0];
sx q[0];
rz(-2.3990193) q[0];
rz(-1.8450129) q[1];
sx q[1];
rz(-0.80755889) q[1];
sx q[1];
rz(-0.47320941) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1593587) q[0];
sx q[0];
rz(-1.8372941) q[0];
sx q[0];
rz(-0.08938163) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8270566) q[2];
sx q[2];
rz(-1.8112931) q[2];
sx q[2];
rz(-2.900689) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0124346) q[1];
sx q[1];
rz(-1.9170463) q[1];
sx q[1];
rz(1.0770256) q[1];
rz(-2.3962184) q[3];
sx q[3];
rz(-1.3826118) q[3];
sx q[3];
rz(0.71063501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5903198) q[2];
sx q[2];
rz(-1.7251622) q[2];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.681916) q[0];
sx q[0];
rz(-1.3688594) q[0];
sx q[0];
rz(-1.3874227) q[0];
rz(0.058902901) q[1];
sx q[1];
rz(-1.7210996) q[1];
sx q[1];
rz(0.65239418) q[1];
rz(-2.0116393) q[2];
sx q[2];
rz(-1.9222353) q[2];
sx q[2];
rz(2.0735379) q[2];
rz(-2.160499) q[3];
sx q[3];
rz(-2.4193939) q[3];
sx q[3];
rz(-1.9388225) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
