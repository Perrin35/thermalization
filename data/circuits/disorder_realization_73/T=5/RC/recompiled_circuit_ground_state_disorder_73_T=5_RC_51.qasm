OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0176528) q[0];
sx q[0];
rz(-2.2247563) q[0];
sx q[0];
rz(0.43489781) q[0];
rz(-1.544156) q[1];
sx q[1];
rz(-0.51900744) q[1];
sx q[1];
rz(-2.5595698) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8343024) q[0];
sx q[0];
rz(-2.8367762) q[0];
sx q[0];
rz(1.0340263) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93272722) q[2];
sx q[2];
rz(-1.1865532) q[2];
sx q[2];
rz(0.24894938) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.37108818) q[1];
sx q[1];
rz(-0.50482115) q[1];
sx q[1];
rz(1.0434898) q[1];
rz(-pi) q[2];
x q[2];
rz(1.852096) q[3];
sx q[3];
rz(-1.3996376) q[3];
sx q[3];
rz(2.8904103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5883098) q[2];
sx q[2];
rz(-1.2574235) q[2];
sx q[2];
rz(2.8498939) q[2];
rz(0.45082539) q[3];
sx q[3];
rz(-0.25142938) q[3];
sx q[3];
rz(-2.5341471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4400103) q[0];
sx q[0];
rz(-2.1511183) q[0];
sx q[0];
rz(-2.0462346) q[0];
rz(-1.1913242) q[1];
sx q[1];
rz(-2.2215863) q[1];
sx q[1];
rz(-0.90257588) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94287017) q[0];
sx q[0];
rz(-0.96317055) q[0];
sx q[0];
rz(0.81815079) q[0];
rz(-pi) q[1];
rz(-2.2392989) q[2];
sx q[2];
rz(-1.8040787) q[2];
sx q[2];
rz(1.4177314) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7831259) q[1];
sx q[1];
rz(-1.2791113) q[1];
sx q[1];
rz(-1.4452219) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9489305) q[3];
sx q[3];
rz(-1.1496667) q[3];
sx q[3];
rz(1.7647183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.1878745) q[2];
sx q[2];
rz(-0.97430054) q[2];
sx q[2];
rz(-0.11889674) q[2];
rz(-2.4925354) q[3];
sx q[3];
rz(-0.18638149) q[3];
sx q[3];
rz(-1.6868748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.651799) q[0];
sx q[0];
rz(-1.8875903) q[0];
sx q[0];
rz(2.6258262) q[0];
rz(-2.0491397) q[1];
sx q[1];
rz(-2.7383883) q[1];
sx q[1];
rz(1.8663503) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7621153) q[0];
sx q[0];
rz(-1.5526875) q[0];
sx q[0];
rz(0.041569592) q[0];
x q[1];
rz(-2.5162773) q[2];
sx q[2];
rz(-0.38039243) q[2];
sx q[2];
rz(1.8521295) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0787492) q[1];
sx q[1];
rz(-0.8991407) q[1];
sx q[1];
rz(0.25154227) q[1];
rz(-0.79408349) q[3];
sx q[3];
rz(-1.6424254) q[3];
sx q[3];
rz(0.41251999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0755997) q[2];
sx q[2];
rz(-1.3174572) q[2];
sx q[2];
rz(-1.9817748) q[2];
rz(0.48948151) q[3];
sx q[3];
rz(-1.5882086) q[3];
sx q[3];
rz(2.421853) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8664261) q[0];
sx q[0];
rz(-2.8186099) q[0];
sx q[0];
rz(-2.6599356) q[0];
rz(2.2109168) q[1];
sx q[1];
rz(-1.296867) q[1];
sx q[1];
rz(0.01893386) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4111644) q[0];
sx q[0];
rz(-1.5207964) q[0];
sx q[0];
rz(-0.020981475) q[0];
x q[1];
rz(-1.6274979) q[2];
sx q[2];
rz(-0.66326521) q[2];
sx q[2];
rz(-2.5692289) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1065797) q[1];
sx q[1];
rz(-1.0570608) q[1];
sx q[1];
rz(0.68961838) q[1];
x q[2];
rz(-0.76374028) q[3];
sx q[3];
rz(-2.1472048) q[3];
sx q[3];
rz(1.9595722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21861741) q[2];
sx q[2];
rz(-0.3564035) q[2];
sx q[2];
rz(-0.74529988) q[2];
rz(-2.7636012) q[3];
sx q[3];
rz(-1.9972921) q[3];
sx q[3];
rz(0.88855827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9208263) q[0];
sx q[0];
rz(-2.8296318) q[0];
sx q[0];
rz(1.9507324) q[0];
rz(-2.6530755) q[1];
sx q[1];
rz(-0.91594511) q[1];
sx q[1];
rz(2.3337505) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5295083) q[0];
sx q[0];
rz(-0.33984646) q[0];
sx q[0];
rz(-1.8456949) q[0];
rz(-1.135181) q[2];
sx q[2];
rz(-2.6382338) q[2];
sx q[2];
rz(-1.8653009) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3563429) q[1];
sx q[1];
rz(-1.6587573) q[1];
sx q[1];
rz(-2.6846519) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9414953) q[3];
sx q[3];
rz(-0.95733023) q[3];
sx q[3];
rz(-0.6544978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0025582) q[2];
sx q[2];
rz(-0.99271861) q[2];
sx q[2];
rz(2.9774418) q[2];
rz(0.74603355) q[3];
sx q[3];
rz(-1.7697216) q[3];
sx q[3];
rz(3.0677838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(0.98110759) q[0];
sx q[0];
rz(-1.8758513) q[0];
sx q[0];
rz(0.72738457) q[0];
rz(-2.6630867) q[1];
sx q[1];
rz(-1.6886657) q[1];
sx q[1];
rz(0.7368288) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7751303) q[0];
sx q[0];
rz(-1.7084165) q[0];
sx q[0];
rz(-0.055995106) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1428035) q[2];
sx q[2];
rz(-2.8341132) q[2];
sx q[2];
rz(-2.3124419) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2347654) q[1];
sx q[1];
rz(-1.9363469) q[1];
sx q[1];
rz(2.2254281) q[1];
x q[2];
rz(1.9584452) q[3];
sx q[3];
rz(-1.5579002) q[3];
sx q[3];
rz(-1.3059592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9937146) q[2];
sx q[2];
rz(-2.1174049) q[2];
sx q[2];
rz(2.4564339) q[2];
rz(2.1327175) q[3];
sx q[3];
rz(-0.69005552) q[3];
sx q[3];
rz(0.25079295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.93609) q[0];
sx q[0];
rz(-3.0476373) q[0];
sx q[0];
rz(-2.0482546) q[0];
rz(3.1178442) q[1];
sx q[1];
rz(-0.7370342) q[1];
sx q[1];
rz(-0.49560961) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4706866) q[0];
sx q[0];
rz(-1.5829594) q[0];
sx q[0];
rz(0.06605102) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43735403) q[2];
sx q[2];
rz(-1.3833356) q[2];
sx q[2];
rz(0.9715434) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.62140761) q[1];
sx q[1];
rz(-1.9193135) q[1];
sx q[1];
rz(1.7640616) q[1];
rz(-1.0823332) q[3];
sx q[3];
rz(-1.1853855) q[3];
sx q[3];
rz(-0.010393532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.89639837) q[2];
sx q[2];
rz(-0.79309016) q[2];
sx q[2];
rz(-0.29655656) q[2];
rz(-0.5101997) q[3];
sx q[3];
rz(-1.6161796) q[3];
sx q[3];
rz(0.27142522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9940014) q[0];
sx q[0];
rz(-2.1864102) q[0];
sx q[0];
rz(1.9872794) q[0];
rz(1.1555903) q[1];
sx q[1];
rz(-2.4842333) q[1];
sx q[1];
rz(1.4591699) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87932779) q[0];
sx q[0];
rz(-1.5478857) q[0];
sx q[0];
rz(2.4266116) q[0];
rz(-pi) q[1];
rz(2.295107) q[2];
sx q[2];
rz(-1.5716388) q[2];
sx q[2];
rz(0.049132012) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3633948) q[1];
sx q[1];
rz(-0.5372592) q[1];
sx q[1];
rz(-1.926653) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1064214) q[3];
sx q[3];
rz(-1.0058306) q[3];
sx q[3];
rz(-1.9244827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83361721) q[2];
sx q[2];
rz(-2.5994382) q[2];
sx q[2];
rz(1.3758434) q[2];
rz(-0.086325072) q[3];
sx q[3];
rz(-2.3089246) q[3];
sx q[3];
rz(1.8858006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.9432705) q[0];
sx q[0];
rz(-2.4477796) q[0];
sx q[0];
rz(-0.86945239) q[0];
rz(0.72932875) q[1];
sx q[1];
rz(-2.2980233) q[1];
sx q[1];
rz(-2.9076911) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1764404) q[0];
sx q[0];
rz(-1.6185221) q[0];
sx q[0];
rz(-2.7244669) q[0];
rz(-pi) q[1];
rz(-0.21017615) q[2];
sx q[2];
rz(-1.4769722) q[2];
sx q[2];
rz(1.2503586) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4471376) q[1];
sx q[1];
rz(-1.7181953) q[1];
sx q[1];
rz(2.6958392) q[1];
rz(-pi) q[2];
rz(-0.31932217) q[3];
sx q[3];
rz(-2.2622613) q[3];
sx q[3];
rz(-2.6806074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.087661155) q[2];
sx q[2];
rz(-1.4940741) q[2];
sx q[2];
rz(1.1762478) q[2];
rz(2.4704399) q[3];
sx q[3];
rz(-1.8920369) q[3];
sx q[3];
rz(1.6254856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.7072356) q[0];
sx q[0];
rz(-0.86152995) q[0];
sx q[0];
rz(1.0435411) q[0];
rz(3.0442944) q[1];
sx q[1];
rz(-1.0568591) q[1];
sx q[1];
rz(-1.1204488) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0138884) q[0];
sx q[0];
rz(-0.59895016) q[0];
sx q[0];
rz(1.4636135) q[0];
x q[1];
rz(-2.3557243) q[2];
sx q[2];
rz(-0.95986569) q[2];
sx q[2];
rz(-0.20139209) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7880558) q[1];
sx q[1];
rz(-0.18419838) q[1];
sx q[1];
rz(-2.8723529) q[1];
rz(-pi) q[2];
rz(-2.9531526) q[3];
sx q[3];
rz(-1.095929) q[3];
sx q[3];
rz(-1.7295854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4572767) q[2];
sx q[2];
rz(-2.5280759) q[2];
sx q[2];
rz(-1.786001) q[2];
rz(-1.5654303) q[3];
sx q[3];
rz(-1.0836982) q[3];
sx q[3];
rz(2.1255597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23715699) q[0];
sx q[0];
rz(-0.56023993) q[0];
sx q[0];
rz(2.3518363) q[0];
rz(-0.65658983) q[1];
sx q[1];
rz(-2.0273392) q[1];
sx q[1];
rz(-0.56710342) q[1];
rz(-0.97548998) q[2];
sx q[2];
rz(-2.3885623) q[2];
sx q[2];
rz(1.9293712) q[2];
rz(1.7651032) q[3];
sx q[3];
rz(-1.8560709) q[3];
sx q[3];
rz(-0.85154497) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
