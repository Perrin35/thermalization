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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2994933) q[0];
sx q[0];
rz(-2.2998739) q[0];
sx q[0];
rz(1.0925348) q[0];
x q[1];
rz(-0.24342033) q[2];
sx q[2];
rz(-1.7924597) q[2];
sx q[2];
rz(-1.8933715) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2185956) q[1];
sx q[1];
rz(-1.6875008) q[1];
sx q[1];
rz(1.7497108) q[1];
rz(1.9919649) q[3];
sx q[3];
rz(-2.2836489) q[3];
sx q[3];
rz(-1.9496535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29229257) q[2];
sx q[2];
rz(-0.92508525) q[2];
sx q[2];
rz(1.383847) q[2];
rz(2.2960831) q[3];
sx q[3];
rz(-0.011761646) q[3];
sx q[3];
rz(-2.1498634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.200835) q[0];
sx q[0];
rz(-0.89460129) q[0];
sx q[0];
rz(2.1064598) q[0];
rz(1.1677405) q[1];
sx q[1];
rz(-0.20226856) q[1];
sx q[1];
rz(-0.77441961) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7709633) q[0];
sx q[0];
rz(-1.1600732) q[0];
sx q[0];
rz(-1.0324917) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2317288) q[2];
sx q[2];
rz(-1.294704) q[2];
sx q[2];
rz(-0.7731437) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6994696) q[1];
sx q[1];
rz(-1.456373) q[1];
sx q[1];
rz(-2.0857986) q[1];
rz(-pi) q[2];
rz(-2.1645427) q[3];
sx q[3];
rz(-3.0809743) q[3];
sx q[3];
rz(-0.82965899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8059798) q[2];
sx q[2];
rz(-0.077294417) q[2];
sx q[2];
rz(2.198931) q[2];
rz(-3.0521657) q[3];
sx q[3];
rz(-0.663203) q[3];
sx q[3];
rz(-2.8360143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5494004) q[0];
sx q[0];
rz(-0.58627272) q[0];
sx q[0];
rz(-2.9495682) q[0];
rz(-2.5281455) q[1];
sx q[1];
rz(-0.44812056) q[1];
sx q[1];
rz(-3.1268069) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0917182) q[0];
sx q[0];
rz(-2.2553068) q[0];
sx q[0];
rz(-2.6189463) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8456782) q[2];
sx q[2];
rz(-2.2408591) q[2];
sx q[2];
rz(2.2704147) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3358128) q[1];
sx q[1];
rz(-1.3381914) q[1];
sx q[1];
rz(2.1628153) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5824541) q[3];
sx q[3];
rz(-1.5069256) q[3];
sx q[3];
rz(3.0272878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6141367) q[2];
sx q[2];
rz(-1.5946031) q[2];
sx q[2];
rz(-0.0097489348) q[2];
rz(-2.5629432) q[3];
sx q[3];
rz(-1.0520244) q[3];
sx q[3];
rz(2.3600522) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-2.9861082) q[1];
sx q[1];
rz(1.3803233) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8603056) q[0];
sx q[0];
rz(-1.5800887) q[0];
sx q[0];
rz(-1.0929266) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99242248) q[2];
sx q[2];
rz(-2.592591) q[2];
sx q[2];
rz(1.8036319) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9050153) q[1];
sx q[1];
rz(-1.6029753) q[1];
sx q[1];
rz(-0.49361146) q[1];
rz(-pi) q[2];
rz(-1.0466688) q[3];
sx q[3];
rz(-1.3635242) q[3];
sx q[3];
rz(0.85974271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.37463793) q[2];
sx q[2];
rz(-1.1344323) q[2];
sx q[2];
rz(3.0881506) q[2];
rz(-0.63353574) q[3];
sx q[3];
rz(-0.41826785) q[3];
sx q[3];
rz(3.1327278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87090129) q[0];
sx q[0];
rz(-0.57796657) q[0];
sx q[0];
rz(2.0722678) q[0];
rz(-1.2879734) q[1];
sx q[1];
rz(-0.063871495) q[1];
sx q[1];
rz(-0.091648253) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1250138) q[0];
sx q[0];
rz(-1.582394) q[0];
sx q[0];
rz(-1.5931104) q[0];
x q[1];
rz(1.1149801) q[2];
sx q[2];
rz(-1.2741003) q[2];
sx q[2];
rz(-1.5199583) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7465072) q[1];
sx q[1];
rz(-1.2068492) q[1];
sx q[1];
rz(-2.0319315) q[1];
rz(-pi) q[2];
x q[2];
rz(1.609747) q[3];
sx q[3];
rz(-0.80128968) q[3];
sx q[3];
rz(-1.0391446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2962467) q[2];
sx q[2];
rz(-0.77719921) q[2];
sx q[2];
rz(-0.11543154) q[2];
rz(-0.7655862) q[3];
sx q[3];
rz(-1.6439532) q[3];
sx q[3];
rz(0.41281858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5033443) q[0];
sx q[0];
rz(-2.4920576) q[0];
sx q[0];
rz(-2.5576538) q[0];
rz(3.0931547) q[1];
sx q[1];
rz(-2.9190639) q[1];
sx q[1];
rz(0.64360523) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9209401) q[0];
sx q[0];
rz(-0.31198129) q[0];
sx q[0];
rz(1.4156154) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90952833) q[2];
sx q[2];
rz(-1.9478746) q[2];
sx q[2];
rz(-2.909735) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4847624) q[1];
sx q[1];
rz(-1.7028013) q[1];
sx q[1];
rz(0.85797347) q[1];
x q[2];
rz(0.19754628) q[3];
sx q[3];
rz(-0.7894844) q[3];
sx q[3];
rz(0.55093605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.863997) q[2];
sx q[2];
rz(-2.3433351) q[2];
sx q[2];
rz(0.64211988) q[2];
rz(-0.13907214) q[3];
sx q[3];
rz(-3.0032447) q[3];
sx q[3];
rz(1.6819491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7794401) q[0];
sx q[0];
rz(-2.7208949) q[0];
sx q[0];
rz(2.5162589) q[0];
rz(-0.23896898) q[1];
sx q[1];
rz(-2.9049554) q[1];
sx q[1];
rz(-1.3561358) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5820393) q[0];
sx q[0];
rz(-1.516321) q[0];
sx q[0];
rz(3.1097058) q[0];
x q[1];
rz(0.37650336) q[2];
sx q[2];
rz(-1.2991221) q[2];
sx q[2];
rz(1.5282549) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7069088) q[1];
sx q[1];
rz(-1.3715991) q[1];
sx q[1];
rz(-1.2046709) q[1];
rz(-pi) q[2];
rz(-0.21276833) q[3];
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
rz(-1.2650547) q[2];
sx q[2];
rz(2.1758101) q[2];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63916373) q[0];
sx q[0];
rz(-0.40488365) q[0];
sx q[0];
rz(-0.13614458) q[0];
rz(2.4038521) q[1];
sx q[1];
rz(-2.9164011) q[1];
sx q[1];
rz(1.9053316) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5571283) q[0];
sx q[0];
rz(-2.4274913) q[0];
sx q[0];
rz(-2.1380928) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8137553) q[2];
sx q[2];
rz(-2.4907673) q[2];
sx q[2];
rz(1.5200159) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7248316) q[1];
sx q[1];
rz(-2.6264852) q[1];
sx q[1];
rz(-2.5558042) q[1];
x q[2];
rz(-1.9916204) q[3];
sx q[3];
rz(-2.4462786) q[3];
sx q[3];
rz(-3.1041404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.69011921) q[2];
sx q[2];
rz(-1.5560919) q[2];
sx q[2];
rz(-2.1915009) q[2];
rz(0.39725605) q[3];
sx q[3];
rz(-2.6451431) q[3];
sx q[3];
rz(2.6955786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5815247) q[0];
sx q[0];
rz(-2.4604605) q[0];
sx q[0];
rz(-0.11823046) q[0];
rz(-1.0826348) q[1];
sx q[1];
rz(-1.2018459) q[1];
sx q[1];
rz(1.3523678) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0467606) q[0];
sx q[0];
rz(-2.7343395) q[0];
sx q[0];
rz(0.26131765) q[0];
rz(-pi) q[1];
rz(2.5973877) q[2];
sx q[2];
rz(-1.6532856) q[2];
sx q[2];
rz(-1.5583041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4115826) q[1];
sx q[1];
rz(-1.6980387) q[1];
sx q[1];
rz(-1.0858658) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0919052) q[3];
sx q[3];
rz(-1.8494431) q[3];
sx q[3];
rz(2.1763733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3236986) q[2];
sx q[2];
rz(-0.69043058) q[2];
sx q[2];
rz(0.77198088) q[2];
rz(-3.0585994) q[3];
sx q[3];
rz(-1.2677742) q[3];
sx q[3];
rz(2.6193589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8817187) q[0];
sx q[0];
rz(-2.3840955) q[0];
sx q[0];
rz(0.74257332) q[0];
rz(-1.2965797) q[1];
sx q[1];
rz(-2.3340338) q[1];
sx q[1];
rz(-0.47320941) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4873692) q[0];
sx q[0];
rz(-0.28074902) q[0];
sx q[0];
rz(-1.8868179) q[0];
x q[1];
rz(1.8270566) q[2];
sx q[2];
rz(-1.8112931) q[2];
sx q[2];
rz(2.900689) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1373945) q[1];
sx q[1];
rz(-0.59474218) q[1];
sx q[1];
rz(-2.221446) q[1];
rz(-pi) q[2];
rz(-0.74537422) q[3];
sx q[3];
rz(-1.7589809) q[3];
sx q[3];
rz(-2.4309576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5903198) q[2];
sx q[2];
rz(-1.4164305) q[2];
sx q[2];
rz(-1.2335221) q[2];
rz(0.34520712) q[3];
sx q[3];
rz(-2.749741) q[3];
sx q[3];
rz(-2.3046618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.681916) q[0];
sx q[0];
rz(-1.7727333) q[0];
sx q[0];
rz(1.75417) q[0];
rz(0.058902901) q[1];
sx q[1];
rz(-1.7210996) q[1];
sx q[1];
rz(0.65239418) q[1];
rz(0.38517135) q[2];
sx q[2];
rz(-1.9829911) q[2];
sx q[2];
rz(0.34172716) q[2];
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
