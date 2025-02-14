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
rz(-1.708344) q[1];
sx q[1];
rz(0.27369174) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84209937) q[0];
sx q[0];
rz(-0.84171879) q[0];
sx q[0];
rz(2.0490579) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8981723) q[2];
sx q[2];
rz(-1.7924597) q[2];
sx q[2];
rz(1.8933715) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.92299709) q[1];
sx q[1];
rz(-1.4540919) q[1];
sx q[1];
rz(-1.7497108) q[1];
rz(-pi) q[2];
rz(-1.9919649) q[3];
sx q[3];
rz(-0.85794373) q[3];
sx q[3];
rz(1.1919392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8493001) q[2];
sx q[2];
rz(-0.92508525) q[2];
sx q[2];
rz(-1.383847) q[2];
rz(2.2960831) q[3];
sx q[3];
rz(-0.011761646) q[3];
sx q[3];
rz(0.99172926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94075769) q[0];
sx q[0];
rz(-2.2469914) q[0];
sx q[0];
rz(1.0351329) q[0];
rz(-1.9738522) q[1];
sx q[1];
rz(-2.9393241) q[1];
sx q[1];
rz(-2.367173) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7073785) q[0];
sx q[0];
rz(-2.0601354) q[0];
sx q[0];
rz(0.46940946) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34458745) q[2];
sx q[2];
rz(-0.93898749) q[2];
sx q[2];
rz(-2.5528204) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.442123) q[1];
sx q[1];
rz(-1.456373) q[1];
sx q[1];
rz(2.0857986) q[1];
x q[2];
rz(1.6210591) q[3];
sx q[3];
rz(-1.6046962) q[3];
sx q[3];
rz(-0.14824319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8059798) q[2];
sx q[2];
rz(-3.0642982) q[2];
sx q[2];
rz(-2.198931) q[2];
rz(3.0521657) q[3];
sx q[3];
rz(-0.663203) q[3];
sx q[3];
rz(2.8360143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59219229) q[0];
sx q[0];
rz(-2.5553199) q[0];
sx q[0];
rz(-2.9495682) q[0];
rz(2.5281455) q[1];
sx q[1];
rz(-2.6934721) q[1];
sx q[1];
rz(-3.1268069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0498745) q[0];
sx q[0];
rz(-2.2553068) q[0];
sx q[0];
rz(0.52264638) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8115571) q[2];
sx q[2];
rz(-0.71612172) q[2];
sx q[2];
rz(0.44490769) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61122429) q[1];
sx q[1];
rz(-2.1448128) q[1];
sx q[1];
rz(-0.27807971) q[1];
rz(-pi) q[2];
rz(0.1802894) q[3];
sx q[3];
rz(-0.064924463) q[3];
sx q[3];
rz(0.29496671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.52745596) q[2];
sx q[2];
rz(-1.5946031) q[2];
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
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0233199) q[0];
sx q[0];
rz(-0.8096205) q[0];
sx q[0];
rz(-0.6672346) q[0];
rz(-1.0046593) q[1];
sx q[1];
rz(-0.15548448) q[1];
sx q[1];
rz(-1.3803233) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2846968) q[0];
sx q[0];
rz(-1.0929489) q[0];
sx q[0];
rz(3.1311281) q[0];
x q[1];
rz(2.0441846) q[2];
sx q[2];
rz(-1.8600827) q[2];
sx q[2];
rz(0.27529596) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7476469) q[1];
sx q[1];
rz(-0.4945728) q[1];
sx q[1];
rz(-3.073758) q[1];
x q[2];
rz(-1.0466688) q[3];
sx q[3];
rz(-1.7780684) q[3];
sx q[3];
rz(2.2818499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.37463793) q[2];
sx q[2];
rz(-1.1344323) q[2];
sx q[2];
rz(0.053442027) q[2];
rz(0.63353574) q[3];
sx q[3];
rz(-2.7233248) q[3];
sx q[3];
rz(3.1327278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87090129) q[0];
sx q[0];
rz(-0.57796657) q[0];
sx q[0];
rz(2.0722678) q[0];
rz(-1.8536192) q[1];
sx q[1];
rz(-3.0777212) q[1];
sx q[1];
rz(3.0499444) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6955513) q[0];
sx q[0];
rz(-1.5484838) q[0];
sx q[0];
rz(-0.011600539) q[0];
x q[1];
rz(-0.32817082) q[2];
sx q[2];
rz(-1.1362906) q[2];
sx q[2];
rz(-0.19319867) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7465072) q[1];
sx q[1];
rz(-1.9347435) q[1];
sx q[1];
rz(-2.0319315) q[1];
x q[2];
rz(-0.76988585) q[3];
sx q[3];
rz(-1.5428233) q[3];
sx q[3];
rz(-2.6370492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.84534591) q[2];
sx q[2];
rz(-0.77719921) q[2];
sx q[2];
rz(-3.0261611) q[2];
rz(2.3760065) q[3];
sx q[3];
rz(-1.4976394) q[3];
sx q[3];
rz(-0.41281858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5033443) q[0];
sx q[0];
rz(-0.64953506) q[0];
sx q[0];
rz(-2.5576538) q[0];
rz(0.04843796) q[1];
sx q[1];
rz(-0.2225288) q[1];
sx q[1];
rz(-2.4979874) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22065255) q[0];
sx q[0];
rz(-0.31198129) q[0];
sx q[0];
rz(-1.7259773) q[0];
rz(-0.46508772) q[2];
sx q[2];
rz(-2.1784867) q[2];
sx q[2];
rz(-1.5236142) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9422655) q[1];
sx q[1];
rz(-0.86546997) q[1];
sx q[1];
rz(-0.17374111) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77966431) q[3];
sx q[3];
rz(-1.4309959) q[3];
sx q[3];
rz(1.1598905) q[3];
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
rz(-0.64211988) q[2];
rz(0.13907214) q[3];
sx q[3];
rz(-3.0032447) q[3];
sx q[3];
rz(1.4596435) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794401) q[0];
sx q[0];
rz(-2.7208949) q[0];
sx q[0];
rz(-2.5162589) q[0];
rz(-0.23896898) q[1];
sx q[1];
rz(-2.9049554) q[1];
sx q[1];
rz(1.7854569) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55955333) q[0];
sx q[0];
rz(-1.6252717) q[0];
sx q[0];
rz(-0.031886851) q[0];
rz(-pi) q[1];
rz(-1.2797592) q[2];
sx q[2];
rz(-1.2087529) q[2];
sx q[2];
rz(0.14823981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.43468388) q[1];
sx q[1];
rz(-1.3715991) q[1];
sx q[1];
rz(-1.2046709) q[1];
rz(1.7271572) q[3];
sx q[3];
rz(-1.3605474) q[3];
sx q[3];
rz(-2.4429136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8465392) q[2];
sx q[2];
rz(-1.2650547) q[2];
sx q[2];
rz(-0.96578252) q[2];
rz(2.8948808) q[3];
sx q[3];
rz(-1.8762981) q[3];
sx q[3];
rz(-1.770796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5024289) q[0];
sx q[0];
rz(-2.736709) q[0];
sx q[0];
rz(0.13614458) q[0];
rz(-2.4038521) q[1];
sx q[1];
rz(-0.22519153) q[1];
sx q[1];
rz(1.9053316) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5571283) q[0];
sx q[0];
rz(-0.71410134) q[0];
sx q[0];
rz(-2.1380928) q[0];
rz(-1.8137553) q[2];
sx q[2];
rz(-2.4907673) q[2];
sx q[2];
rz(1.5200159) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7248316) q[1];
sx q[1];
rz(-0.5151075) q[1];
sx q[1];
rz(0.58578844) q[1];
rz(-0.32847877) q[3];
sx q[3];
rz(-2.1953479) q[3];
sx q[3];
rz(0.49027944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69011921) q[2];
sx q[2];
rz(-1.5560919) q[2];
sx q[2];
rz(-2.1915009) q[2];
rz(-0.39725605) q[3];
sx q[3];
rz(-2.6451431) q[3];
sx q[3];
rz(-2.6955786) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5600679) q[0];
sx q[0];
rz(-0.6811322) q[0];
sx q[0];
rz(-0.11823046) q[0];
rz(-1.0826348) q[1];
sx q[1];
rz(-1.2018459) q[1];
sx q[1];
rz(-1.7892249) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7633469) q[0];
sx q[0];
rz(-1.9634501) q[0];
sx q[0];
rz(-1.6817842) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1583516) q[2];
sx q[2];
rz(-2.5917946) q[2];
sx q[2];
rz(2.9937772) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.39555672) q[1];
sx q[1];
rz(-2.6415351) q[1];
sx q[1];
rz(-1.3029424) q[1];
rz(-pi) q[2];
rz(-0.31183536) q[3];
sx q[3];
rz(-1.1118299) q[3];
sx q[3];
rz(-0.74742693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.81789404) q[2];
sx q[2];
rz(-0.69043058) q[2];
sx q[2];
rz(2.3696118) q[2];
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
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-2.8817187) q[0];
sx q[0];
rz(-2.3840955) q[0];
sx q[0];
rz(0.74257332) q[0];
rz(1.2965797) q[1];
sx q[1];
rz(-2.3340338) q[1];
sx q[1];
rz(0.47320941) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1593587) q[0];
sx q[0];
rz(-1.3042985) q[0];
sx q[0];
rz(3.052211) q[0];
rz(-pi) q[1];
rz(1.8270566) q[2];
sx q[2];
rz(-1.3302996) q[2];
sx q[2];
rz(0.24090361) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1373945) q[1];
sx q[1];
rz(-2.5468505) q[1];
sx q[1];
rz(-2.221446) q[1];
x q[2];
rz(1.3172174) q[3];
sx q[3];
rz(-2.2999956) q[3];
sx q[3];
rz(1.031145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.55127281) q[2];
sx q[2];
rz(-1.4164305) q[2];
sx q[2];
rz(1.2335221) q[2];
rz(-0.34520712) q[3];
sx q[3];
rz(-0.39185169) q[3];
sx q[3];
rz(0.83693081) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.681916) q[0];
sx q[0];
rz(-1.7727333) q[0];
sx q[0];
rz(1.75417) q[0];
rz(-0.058902901) q[1];
sx q[1];
rz(-1.420493) q[1];
sx q[1];
rz(-2.4891985) q[1];
rz(-0.38517135) q[2];
sx q[2];
rz(-1.1586015) q[2];
sx q[2];
rz(-2.7998655) q[2];
rz(-2.6860438) q[3];
sx q[3];
rz(-2.1524317) q[3];
sx q[3];
rz(0.47453415) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
