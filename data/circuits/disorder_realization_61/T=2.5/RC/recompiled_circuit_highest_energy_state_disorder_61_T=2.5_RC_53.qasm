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
rz(-2.6391368) q[0];
sx q[0];
rz(-1.1075736) q[0];
sx q[0];
rz(-1.3753608) q[0];
rz(2.8227168) q[1];
sx q[1];
rz(-1.5006637) q[1];
sx q[1];
rz(-0.73135102) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36424822) q[0];
sx q[0];
rz(-0.70667446) q[0];
sx q[0];
rz(1.0663435) q[0];
x q[1];
rz(-0.66197239) q[2];
sx q[2];
rz(-1.5675002) q[2];
sx q[2];
rz(-2.4087083) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7544985) q[1];
sx q[1];
rz(-1.6162786) q[1];
sx q[1];
rz(-2.3552853) q[1];
x q[2];
rz(-2.6177789) q[3];
sx q[3];
rz(-2.9546545) q[3];
sx q[3];
rz(-2.4863249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2607164) q[2];
sx q[2];
rz(-1.9274351) q[2];
sx q[2];
rz(-2.1991849) q[2];
rz(-2.2830394) q[3];
sx q[3];
rz(-2.5311354) q[3];
sx q[3];
rz(2.4220991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061676625) q[0];
sx q[0];
rz(-2.583857) q[0];
sx q[0];
rz(-0.058636531) q[0];
rz(0.05642852) q[1];
sx q[1];
rz(-1.7516878) q[1];
sx q[1];
rz(2.0243534) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4898713) q[0];
sx q[0];
rz(-0.84017032) q[0];
sx q[0];
rz(1.0307167) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12877474) q[2];
sx q[2];
rz(-1.2733545) q[2];
sx q[2];
rz(2.8853438) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8630115) q[1];
sx q[1];
rz(-1.3323235) q[1];
sx q[1];
rz(3.0568074) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9100203) q[3];
sx q[3];
rz(-0.48677126) q[3];
sx q[3];
rz(2.6269905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3650018) q[2];
sx q[2];
rz(-0.33010179) q[2];
sx q[2];
rz(-2.7351725) q[2];
rz(0.31600076) q[3];
sx q[3];
rz(-1.4603115) q[3];
sx q[3];
rz(2.3587904) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86421788) q[0];
sx q[0];
rz(-0.45806956) q[0];
sx q[0];
rz(1.8055441) q[0];
rz(-0.27851963) q[1];
sx q[1];
rz(-1.3986992) q[1];
sx q[1];
rz(-2.1727402) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1135546) q[0];
sx q[0];
rz(-0.73532205) q[0];
sx q[0];
rz(-2.5707629) q[0];
rz(-pi) q[1];
rz(-0.55726643) q[2];
sx q[2];
rz(-0.87199482) q[2];
sx q[2];
rz(-3.0521309) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.75894657) q[1];
sx q[1];
rz(-0.29932705) q[1];
sx q[1];
rz(-1.4274197) q[1];
x q[2];
rz(-2.8519451) q[3];
sx q[3];
rz(-2.395219) q[3];
sx q[3];
rz(-2.1044855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1442673) q[2];
sx q[2];
rz(-0.48339016) q[2];
sx q[2];
rz(2.9664795) q[2];
rz(-1.8363606) q[3];
sx q[3];
rz(-0.94279083) q[3];
sx q[3];
rz(1.85873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1433732) q[0];
sx q[0];
rz(-2.0046736) q[0];
sx q[0];
rz(1.0293707) q[0];
rz(0.81946212) q[1];
sx q[1];
rz(-1.1623323) q[1];
sx q[1];
rz(-2.3447461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7356257) q[0];
sx q[0];
rz(-1.7628364) q[0];
sx q[0];
rz(2.6940932) q[0];
x q[1];
rz(1.5605075) q[2];
sx q[2];
rz(-2.6136189) q[2];
sx q[2];
rz(-0.99425232) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.5552434) q[1];
sx q[1];
rz(-1.8420911) q[1];
sx q[1];
rz(0.21315141) q[1];
rz(2.4791777) q[3];
sx q[3];
rz(-1.9180505) q[3];
sx q[3];
rz(-1.9062689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0113819) q[2];
sx q[2];
rz(-3.104976) q[2];
sx q[2];
rz(1.2288564) q[2];
rz(2.7025488) q[3];
sx q[3];
rz(-1.5758347) q[3];
sx q[3];
rz(-1.9108093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48290408) q[0];
sx q[0];
rz(-1.6084325) q[0];
sx q[0];
rz(-0.80832344) q[0];
rz(1.3574379) q[1];
sx q[1];
rz(-2.3517377) q[1];
sx q[1];
rz(-2.435991) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1418512) q[0];
sx q[0];
rz(-1.6054285) q[0];
sx q[0];
rz(-2.8409728) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86807495) q[2];
sx q[2];
rz(-1.3931615) q[2];
sx q[2];
rz(-2.9736819) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8452705) q[1];
sx q[1];
rz(-0.25320881) q[1];
sx q[1];
rz(1.5814387) q[1];
rz(-pi) q[2];
rz(-2.5217149) q[3];
sx q[3];
rz(-2.2408463) q[3];
sx q[3];
rz(2.1170962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6610873) q[2];
sx q[2];
rz(-1.2837807) q[2];
sx q[2];
rz(-2.535848) q[2];
rz(0.12039603) q[3];
sx q[3];
rz(-1.2984637) q[3];
sx q[3];
rz(-2.3275183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8814988) q[0];
sx q[0];
rz(-0.81387481) q[0];
sx q[0];
rz(3.1392198) q[0];
rz(2.782605) q[1];
sx q[1];
rz(-1.2009883) q[1];
sx q[1];
rz(2.733309) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4160894) q[0];
sx q[0];
rz(-1.6973087) q[0];
sx q[0];
rz(-3.1055272) q[0];
x q[1];
rz(-0.26770182) q[2];
sx q[2];
rz(-2.1609339) q[2];
sx q[2];
rz(2.6538284) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5376512) q[1];
sx q[1];
rz(-1.3408325) q[1];
sx q[1];
rz(1.5020788) q[1];
rz(-pi) q[2];
rz(-0.8114641) q[3];
sx q[3];
rz(-2.0272019) q[3];
sx q[3];
rz(-1.2450258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6458873) q[2];
sx q[2];
rz(-2.3453823) q[2];
sx q[2];
rz(3.0118946) q[2];
rz(0.4194704) q[3];
sx q[3];
rz(-1.9286112) q[3];
sx q[3];
rz(2.3278842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7523338) q[0];
sx q[0];
rz(-2.3663754) q[0];
sx q[0];
rz(-2.4580521) q[0];
rz(-0.93841249) q[1];
sx q[1];
rz(-1.7529731) q[1];
sx q[1];
rz(1.9155115) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6594567) q[0];
sx q[0];
rz(-0.99471751) q[0];
sx q[0];
rz(-2.0123737) q[0];
rz(-pi) q[1];
rz(1.3834882) q[2];
sx q[2];
rz(-2.1626327) q[2];
sx q[2];
rz(-2.1746152) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.214512) q[1];
sx q[1];
rz(-2.0190372) q[1];
sx q[1];
rz(2.8238676) q[1];
x q[2];
rz(-0.55415537) q[3];
sx q[3];
rz(-0.99923493) q[3];
sx q[3];
rz(-3.1031648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40743264) q[2];
sx q[2];
rz(-1.7075044) q[2];
sx q[2];
rz(0.70719353) q[2];
rz(-1.4736942) q[3];
sx q[3];
rz(-0.67631045) q[3];
sx q[3];
rz(1.7128806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68607512) q[0];
sx q[0];
rz(-2.7954743) q[0];
sx q[0];
rz(-0.92380512) q[0];
rz(-1.0651945) q[1];
sx q[1];
rz(-1.9520452) q[1];
sx q[1];
rz(-2.0069897) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3649639) q[0];
sx q[0];
rz(-2.7641649) q[0];
sx q[0];
rz(-1.6799742) q[0];
rz(-pi) q[1];
rz(-2.2840225) q[2];
sx q[2];
rz(-2.3152707) q[2];
sx q[2];
rz(-0.37304953) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75443134) q[1];
sx q[1];
rz(-0.5335156) q[1];
sx q[1];
rz(1.4656161) q[1];
rz(-pi) q[2];
rz(0.60875684) q[3];
sx q[3];
rz(-1.8894102) q[3];
sx q[3];
rz(0.99005885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0040466641) q[2];
sx q[2];
rz(-1.4038439) q[2];
sx q[2];
rz(-2.498632) q[2];
rz(0.17063394) q[3];
sx q[3];
rz(-2.4891487) q[3];
sx q[3];
rz(-2.046106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74649015) q[0];
sx q[0];
rz(-3.009142) q[0];
sx q[0];
rz(-2.3353031) q[0];
rz(-0.0099446615) q[1];
sx q[1];
rz(-0.6178304) q[1];
sx q[1];
rz(2.8795805) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1073955) q[0];
sx q[0];
rz(-2.5881564) q[0];
sx q[0];
rz(1.8897927) q[0];
rz(-pi) q[1];
rz(-1.7315032) q[2];
sx q[2];
rz(-0.93148684) q[2];
sx q[2];
rz(-2.8455545) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1478473) q[1];
sx q[1];
rz(-2.6585852) q[1];
sx q[1];
rz(1.1149393) q[1];
rz(-2.8132961) q[3];
sx q[3];
rz(-2.0496164) q[3];
sx q[3];
rz(1.6678068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.749873) q[2];
sx q[2];
rz(-2.5842857) q[2];
sx q[2];
rz(-0.14287512) q[2];
rz(2.2376132) q[3];
sx q[3];
rz(-1.6429792) q[3];
sx q[3];
rz(-2.2601295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3751462) q[0];
sx q[0];
rz(-1.2297933) q[0];
sx q[0];
rz(-2.4850856) q[0];
rz(0.87908602) q[1];
sx q[1];
rz(-1.1708941) q[1];
sx q[1];
rz(2.5118929) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.499959) q[0];
sx q[0];
rz(-0.90454209) q[0];
sx q[0];
rz(2.0883661) q[0];
rz(-pi) q[1];
rz(-0.26816396) q[2];
sx q[2];
rz(-2.0879474) q[2];
sx q[2];
rz(-1.6059396) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.15805298) q[1];
sx q[1];
rz(-1/(6*pi)) q[1];
sx q[1];
rz(-2.1237462) q[1];
x q[2];
rz(-1.0694169) q[3];
sx q[3];
rz(-0.58749858) q[3];
sx q[3];
rz(-2.0794433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62275824) q[2];
sx q[2];
rz(-2.3164985) q[2];
sx q[2];
rz(-0.68010124) q[2];
rz(-0.71637362) q[3];
sx q[3];
rz(-0.10321897) q[3];
sx q[3];
rz(-1.1187925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5730561) q[0];
sx q[0];
rz(-1.9087044) q[0];
sx q[0];
rz(-2.9317324) q[0];
rz(-2.9987891) q[1];
sx q[1];
rz(-1.0719943) q[1];
sx q[1];
rz(0.73678585) q[1];
rz(-1.8129195) q[2];
sx q[2];
rz(-1.6264781) q[2];
sx q[2];
rz(-1.3774058) q[2];
rz(1.5132001) q[3];
sx q[3];
rz(-0.61020281) q[3];
sx q[3];
rz(-1.4635066) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
