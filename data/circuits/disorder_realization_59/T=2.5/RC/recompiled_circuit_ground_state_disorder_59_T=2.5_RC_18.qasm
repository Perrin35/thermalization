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
rz(2.7432888) q[0];
sx q[0];
rz(9.7437133) q[0];
rz(2.6755264) q[1];
sx q[1];
rz(-1.4332486) q[1];
sx q[1];
rz(-0.27369174) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6382054) q[0];
sx q[0];
rz(-2.2944258) q[0];
sx q[0];
rz(-2.6658325) q[0];
rz(-pi) q[1];
rz(0.24342033) q[2];
sx q[2];
rz(-1.7924597) q[2];
sx q[2];
rz(1.8933715) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66885405) q[1];
sx q[1];
rz(-1.7484807) q[1];
sx q[1];
rz(3.0230126) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9919649) q[3];
sx q[3];
rz(-2.2836489) q[3];
sx q[3];
rz(1.9496535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29229257) q[2];
sx q[2];
rz(-0.92508525) q[2];
sx q[2];
rz(-1.383847) q[2];
rz(0.84550953) q[3];
sx q[3];
rz(-0.011761646) q[3];
sx q[3];
rz(-0.99172926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94075769) q[0];
sx q[0];
rz(-0.89460129) q[0];
sx q[0];
rz(2.1064598) q[0];
rz(1.9738522) q[1];
sx q[1];
rz(-2.9393241) q[1];
sx q[1];
rz(2.367173) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5308038) q[0];
sx q[0];
rz(-2.4770081) q[0];
sx q[0];
rz(-2.2749645) q[0];
rz(-0.34458745) q[2];
sx q[2];
rz(-0.93898749) q[2];
sx q[2];
rz(0.58877221) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.442123) q[1];
sx q[1];
rz(-1.6852196) q[1];
sx q[1];
rz(-2.0857986) q[1];
rz(-pi) q[2];
rz(1.6210591) q[3];
sx q[3];
rz(-1.6046962) q[3];
sx q[3];
rz(2.9933495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3356129) q[2];
sx q[2];
rz(-0.077294417) q[2];
sx q[2];
rz(0.94266164) q[2];
rz(-0.089426905) q[3];
sx q[3];
rz(-2.4783897) q[3];
sx q[3];
rz(0.30557835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5494004) q[0];
sx q[0];
rz(-2.5553199) q[0];
sx q[0];
rz(-0.19202448) q[0];
rz(2.5281455) q[1];
sx q[1];
rz(-2.6934721) q[1];
sx q[1];
rz(-3.1268069) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0498745) q[0];
sx q[0];
rz(-0.8862859) q[0];
sx q[0];
rz(0.52264638) q[0];
rz(2.8115571) q[2];
sx q[2];
rz(-0.71612172) q[2];
sx q[2];
rz(0.44490769) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.61122429) q[1];
sx q[1];
rz(-2.1448128) q[1];
sx q[1];
rz(0.27807971) q[1];
rz(-pi) q[2];
rz(1.5591386) q[3];
sx q[3];
rz(-1.5069256) q[3];
sx q[3];
rz(0.11430489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6141367) q[2];
sx q[2];
rz(-1.5946031) q[2];
sx q[2];
rz(0.0097489348) q[2];
rz(-0.57864946) q[3];
sx q[3];
rz(-2.0895683) q[3];
sx q[3];
rz(-0.78154045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0233199) q[0];
sx q[0];
rz(-0.8096205) q[0];
sx q[0];
rz(-2.4743581) q[0];
rz(-1.0046593) q[1];
sx q[1];
rz(-2.9861082) q[1];
sx q[1];
rz(1.3803233) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2846968) q[0];
sx q[0];
rz(-1.0929489) q[0];
sx q[0];
rz(-0.010464593) q[0];
rz(-pi) q[1];
rz(2.818872) q[2];
sx q[2];
rz(-2.0230133) q[2];
sx q[2];
rz(1.9911901) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.31691027) q[1];
sx q[1];
rz(-2.0641293) q[1];
sx q[1];
rz(-1.534259) q[1];
x q[2];
rz(-2.9033086) q[3];
sx q[3];
rz(-1.0589979) q[3];
sx q[3];
rz(0.82945591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7669547) q[2];
sx q[2];
rz(-1.1344323) q[2];
sx q[2];
rz(0.053442027) q[2];
rz(-2.5080569) q[3];
sx q[3];
rz(-2.7233248) q[3];
sx q[3];
rz(3.1327278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87090129) q[0];
sx q[0];
rz(-0.57796657) q[0];
sx q[0];
rz(-1.0693249) q[0];
rz(-1.8536192) q[1];
sx q[1];
rz(-0.063871495) q[1];
sx q[1];
rz(0.091648253) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033460334) q[0];
sx q[0];
rz(-0.025147557) q[0];
sx q[0];
rz(1.0914241) q[0];
rz(-0.32817082) q[2];
sx q[2];
rz(-2.0053021) q[2];
sx q[2];
rz(-2.948394) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7908259) q[1];
sx q[1];
rz(-1.1419527) q[1];
sx q[1];
rz(0.40216202) q[1];
x q[2];
rz(2.3717068) q[3];
sx q[3];
rz(-1.5987694) q[3];
sx q[3];
rz(-0.5045435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2962467) q[2];
sx q[2];
rz(-0.77719921) q[2];
sx q[2];
rz(0.11543154) q[2];
rz(-2.3760065) q[3];
sx q[3];
rz(-1.6439532) q[3];
sx q[3];
rz(-0.41281858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63824832) q[0];
sx q[0];
rz(-0.64953506) q[0];
sx q[0];
rz(0.58393884) q[0];
rz(-0.04843796) q[1];
sx q[1];
rz(-2.9190639) q[1];
sx q[1];
rz(0.64360523) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2023425) q[0];
sx q[0];
rz(-1.5233375) q[0];
sx q[0];
rz(-1.8792633) q[0];
rz(-0.90952833) q[2];
sx q[2];
rz(-1.9478746) q[2];
sx q[2];
rz(2.909735) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4847624) q[1];
sx q[1];
rz(-1.4387913) q[1];
sx q[1];
rz(2.2836192) q[1];
rz(-pi) q[2];
rz(0.19754628) q[3];
sx q[3];
rz(-2.3521083) q[3];
sx q[3];
rz(-0.55093605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.27759564) q[2];
sx q[2];
rz(-0.79825753) q[2];
sx q[2];
rz(-2.4994728) q[2];
rz(0.13907214) q[3];
sx q[3];
rz(-3.0032447) q[3];
sx q[3];
rz(-1.6819491) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794401) q[0];
sx q[0];
rz(-0.42069778) q[0];
sx q[0];
rz(2.5162589) q[0];
rz(2.9026237) q[1];
sx q[1];
rz(-2.9049554) q[1];
sx q[1];
rz(1.7854569) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029615868) q[0];
sx q[0];
rz(-3.0784791) q[0];
sx q[0];
rz(-2.099865) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37650336) q[2];
sx q[2];
rz(-1.8424705) q[2];
sx q[2];
rz(-1.5282549) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.060387) q[1];
sx q[1];
rz(-1.2122417) q[1];
sx q[1];
rz(0.21292673) q[1];
rz(-pi) q[2];
rz(-1.4144355) q[3];
sx q[3];
rz(-1.3605474) q[3];
sx q[3];
rz(-2.4429136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29505342) q[2];
sx q[2];
rz(-1.876538) q[2];
sx q[2];
rz(0.96578252) q[2];
rz(-0.24671181) q[3];
sx q[3];
rz(-1.2652946) q[3];
sx q[3];
rz(1.770796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5024289) q[0];
sx q[0];
rz(-2.736709) q[0];
sx q[0];
rz(3.0054481) q[0];
rz(0.73774058) q[1];
sx q[1];
rz(-2.9164011) q[1];
sx q[1];
rz(1.236261) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6039918) q[0];
sx q[0];
rz(-1.9304327) q[0];
sx q[0];
rz(2.2019772) q[0];
rz(-pi) q[1];
rz(0.18119104) q[2];
sx q[2];
rz(-0.94215067) q[2];
sx q[2];
rz(-1.218007) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.41676109) q[1];
sx q[1];
rz(-2.6264852) q[1];
sx q[1];
rz(-0.58578844) q[1];
x q[2];
rz(-2.2216283) q[3];
sx q[3];
rz(-1.3060089) q[3];
sx q[3];
rz(1.8643606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69011921) q[2];
sx q[2];
rz(-1.5855007) q[2];
sx q[2];
rz(-0.95009178) q[2];
rz(-0.39725605) q[3];
sx q[3];
rz(-0.49644956) q[3];
sx q[3];
rz(-0.44601405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5600679) q[0];
sx q[0];
rz(-2.4604605) q[0];
sx q[0];
rz(0.11823046) q[0];
rz(1.0826348) q[1];
sx q[1];
rz(-1.2018459) q[1];
sx q[1];
rz(-1.3523678) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0948321) q[0];
sx q[0];
rz(-2.7343395) q[0];
sx q[0];
rz(-2.880275) q[0];
rz(-2.5973877) q[2];
sx q[2];
rz(-1.4883071) q[2];
sx q[2];
rz(-1.5583041) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.73001) q[1];
sx q[1];
rz(-1.4435539) q[1];
sx q[1];
rz(2.0557269) q[1];
rz(-pi) q[2];
rz(-2.1264137) q[3];
sx q[3];
rz(-0.54856442) q[3];
sx q[3];
rz(-3.0231904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.81789404) q[2];
sx q[2];
rz(-0.69043058) q[2];
sx q[2];
rz(-0.77198088) q[2];
rz(-3.0585994) q[3];
sx q[3];
rz(-1.2677742) q[3];
sx q[3];
rz(2.6193589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25987396) q[0];
sx q[0];
rz(-0.75749713) q[0];
sx q[0];
rz(-0.74257332) q[0];
rz(-1.2965797) q[1];
sx q[1];
rz(-0.80755889) q[1];
sx q[1];
rz(-2.6683832) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61215996) q[0];
sx q[0];
rz(-1.4845779) q[0];
sx q[0];
rz(-1.8383121) q[0];
rz(1.3145361) q[2];
sx q[2];
rz(-1.8112931) q[2];
sx q[2];
rz(-2.900689) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1291581) q[1];
sx q[1];
rz(-1.9170463) q[1];
sx q[1];
rz(1.0770256) q[1];
rz(-pi) q[2];
rz(2.8678611) q[3];
sx q[3];
rz(-2.377284) q[3];
sx q[3];
rz(2.4814062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5903198) q[2];
sx q[2];
rz(-1.7251622) q[2];
sx q[2];
rz(1.9080706) q[2];
rz(0.34520712) q[3];
sx q[3];
rz(-0.39185169) q[3];
sx q[3];
rz(-0.83693081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45967669) q[0];
sx q[0];
rz(-1.3688594) q[0];
sx q[0];
rz(-1.3874227) q[0];
rz(-3.0826898) q[1];
sx q[1];
rz(-1.7210996) q[1];
sx q[1];
rz(0.65239418) q[1];
rz(-2.7564213) q[2];
sx q[2];
rz(-1.9829911) q[2];
sx q[2];
rz(0.34172716) q[2];
rz(0.4555489) q[3];
sx q[3];
rz(-2.1524317) q[3];
sx q[3];
rz(0.47453415) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
