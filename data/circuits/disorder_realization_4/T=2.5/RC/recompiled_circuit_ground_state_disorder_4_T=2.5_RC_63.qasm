OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57920116) q[0];
sx q[0];
rz(-0.76051036) q[0];
sx q[0];
rz(1.0150681) q[0];
rz(-1.1633582) q[1];
sx q[1];
rz(-2.7203163) q[1];
sx q[1];
rz(-1.2383229) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9559798) q[0];
sx q[0];
rz(-2.2816814) q[0];
sx q[0];
rz(0.86998633) q[0];
x q[1];
rz(-0.86482817) q[2];
sx q[2];
rz(-2.6061686) q[2];
sx q[2];
rz(-3.0147417) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.94436344) q[1];
sx q[1];
rz(-0.60603332) q[1];
sx q[1];
rz(-2.2418749) q[1];
rz(2.5776776) q[3];
sx q[3];
rz(-1.7132572) q[3];
sx q[3];
rz(3.0092653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6282661) q[2];
sx q[2];
rz(-1.3506177) q[2];
sx q[2];
rz(-0.70293054) q[2];
rz(-2.2859196) q[3];
sx q[3];
rz(-2.296505) q[3];
sx q[3];
rz(1.2969016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4502451) q[0];
sx q[0];
rz(-2.0424728) q[0];
sx q[0];
rz(-0.74478373) q[0];
rz(-1.5738457) q[1];
sx q[1];
rz(-2.4796922) q[1];
sx q[1];
rz(0.94211284) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5013393) q[0];
sx q[0];
rz(-0.82385175) q[0];
sx q[0];
rz(1.8477316) q[0];
rz(-1.6555384) q[2];
sx q[2];
rz(-0.7543482) q[2];
sx q[2];
rz(-0.95065439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.86650634) q[1];
sx q[1];
rz(-0.77436354) q[1];
sx q[1];
rz(-1.3759717) q[1];
rz(-pi) q[2];
rz(-0.09047507) q[3];
sx q[3];
rz(-2.5535085) q[3];
sx q[3];
rz(-0.17234853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.88910237) q[2];
sx q[2];
rz(-2.1847051) q[2];
sx q[2];
rz(1.0955742) q[2];
rz(-1.6628294) q[3];
sx q[3];
rz(-2.3507599) q[3];
sx q[3];
rz(1.7553294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11109322) q[0];
sx q[0];
rz(-1.4249304) q[0];
sx q[0];
rz(-1.4053364) q[0];
rz(1.3966712) q[1];
sx q[1];
rz(-1.3537355) q[1];
sx q[1];
rz(-0.90000802) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25728658) q[0];
sx q[0];
rz(-2.2955756) q[0];
sx q[0];
rz(-1.0777362) q[0];
rz(-pi) q[1];
rz(-0.75863691) q[2];
sx q[2];
rz(-0.73497811) q[2];
sx q[2];
rz(-1.0204362) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.450728) q[1];
sx q[1];
rz(-0.83800661) q[1];
sx q[1];
rz(1.8045252) q[1];
rz(-2.1846287) q[3];
sx q[3];
rz(-1.004289) q[3];
sx q[3];
rz(1.5885799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6806543) q[2];
sx q[2];
rz(-0.21549455) q[2];
sx q[2];
rz(-0.94949618) q[2];
rz(0.62120581) q[3];
sx q[3];
rz(-0.92531365) q[3];
sx q[3];
rz(1.1533823) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7200274) q[0];
sx q[0];
rz(-1.8864487) q[0];
sx q[0];
rz(-0.048728745) q[0];
rz(0.85982927) q[1];
sx q[1];
rz(-0.33165926) q[1];
sx q[1];
rz(-2.8299832) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6831419) q[0];
sx q[0];
rz(-1.8425757) q[0];
sx q[0];
rz(-1.0807476) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34398244) q[2];
sx q[2];
rz(-1.0687573) q[2];
sx q[2];
rz(0.12884049) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3606122) q[1];
sx q[1];
rz(-2.2192839) q[1];
sx q[1];
rz(2.8070252) q[1];
x q[2];
rz(2.1457003) q[3];
sx q[3];
rz(-0.51255915) q[3];
sx q[3];
rz(2.9957774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9638046) q[2];
sx q[2];
rz(-2.9298941) q[2];
sx q[2];
rz(1.4908028) q[2];
rz(-0.35495159) q[3];
sx q[3];
rz(-1.7345411) q[3];
sx q[3];
rz(-1.8420334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0896924) q[0];
sx q[0];
rz(-2.7841452) q[0];
sx q[0];
rz(-1.8748913) q[0];
rz(-3.025324) q[1];
sx q[1];
rz(-0.99028844) q[1];
sx q[1];
rz(-1.6273392) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26402347) q[0];
sx q[0];
rz(-2.2248587) q[0];
sx q[0];
rz(-0.39790776) q[0];
x q[1];
rz(-0.090574663) q[2];
sx q[2];
rz(-1.3514869) q[2];
sx q[2];
rz(-2.5136869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7092412) q[1];
sx q[1];
rz(-0.54644692) q[1];
sx q[1];
rz(-1.6691471) q[1];
rz(-pi) q[2];
rz(2.299304) q[3];
sx q[3];
rz(-1.5986048) q[3];
sx q[3];
rz(0.14735809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9194455) q[2];
sx q[2];
rz(-0.47480348) q[2];
sx q[2];
rz(-1.9661281) q[2];
rz(-2.2373824) q[3];
sx q[3];
rz(-1.892482) q[3];
sx q[3];
rz(-0.97833943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0670052) q[0];
sx q[0];
rz(-0.33127221) q[0];
sx q[0];
rz(-2.7401155) q[0];
rz(2.0145156) q[1];
sx q[1];
rz(-1.3104442) q[1];
sx q[1];
rz(-2.3289767) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3881587) q[0];
sx q[0];
rz(-1.5666612) q[0];
sx q[0];
rz(-2.2049516) q[0];
x q[1];
rz(0.27490669) q[2];
sx q[2];
rz(-1.3454352) q[2];
sx q[2];
rz(1.7079086) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.26906313) q[1];
sx q[1];
rz(-1.3663962) q[1];
sx q[1];
rz(-0.60592954) q[1];
rz(-pi) q[2];
rz(-0.74515588) q[3];
sx q[3];
rz(-2.0260915) q[3];
sx q[3];
rz(-0.32768341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.20069417) q[2];
sx q[2];
rz(-2.5914067) q[2];
sx q[2];
rz(-1.5563439) q[2];
rz(2.9511792) q[3];
sx q[3];
rz(-2.3868581) q[3];
sx q[3];
rz(-1.2079027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.3170526) q[0];
sx q[0];
rz(-0.38924488) q[0];
sx q[0];
rz(-3.0188766) q[0];
rz(1.9484733) q[1];
sx q[1];
rz(-0.90843186) q[1];
sx q[1];
rz(2.3741123) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7372663) q[0];
sx q[0];
rz(-0.64868673) q[0];
sx q[0];
rz(-1.655633) q[0];
rz(-0.20920472) q[2];
sx q[2];
rz(-2.7520617) q[2];
sx q[2];
rz(-0.65702932) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1134935) q[1];
sx q[1];
rz(-0.65526795) q[1];
sx q[1];
rz(-1.6428309) q[1];
x q[2];
rz(3.1311099) q[3];
sx q[3];
rz(-0.40626486) q[3];
sx q[3];
rz(2.3113113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3976589) q[2];
sx q[2];
rz(-3.1199516) q[2];
sx q[2];
rz(1.5896612) q[2];
rz(-1.4895561) q[3];
sx q[3];
rz(-1.7347615) q[3];
sx q[3];
rz(2.0261197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81812304) q[0];
sx q[0];
rz(-2.7796845) q[0];
sx q[0];
rz(-0.37539151) q[0];
rz(-0.88343945) q[1];
sx q[1];
rz(-1.6049623) q[1];
sx q[1];
rz(2.4129131) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.323384) q[0];
sx q[0];
rz(-2.2012156) q[0];
sx q[0];
rz(-0.66108836) q[0];
x q[1];
rz(0.23587464) q[2];
sx q[2];
rz(-1.4067003) q[2];
sx q[2];
rz(-0.99170384) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2051712) q[1];
sx q[1];
rz(-1.2146307) q[1];
sx q[1];
rz(2.2486399) q[1];
rz(-pi) q[2];
rz(-2.3865984) q[3];
sx q[3];
rz(-1.0858616) q[3];
sx q[3];
rz(-0.86673966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6931307) q[2];
sx q[2];
rz(-1.5631661) q[2];
sx q[2];
rz(0.7473839) q[2];
rz(1.4061617) q[3];
sx q[3];
rz(-1.9236671) q[3];
sx q[3];
rz(1.5555443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2324227) q[0];
sx q[0];
rz(-2.2012043) q[0];
sx q[0];
rz(-0.7255834) q[0];
rz(1.3369417) q[1];
sx q[1];
rz(-1.888211) q[1];
sx q[1];
rz(0.57428378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1385571) q[0];
sx q[0];
rz(-1.7302365) q[0];
sx q[0];
rz(-0.11388679) q[0];
x q[1];
rz(-1.0011739) q[2];
sx q[2];
rz(-1.2204079) q[2];
sx q[2];
rz(-0.97942615) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6114823) q[1];
sx q[1];
rz(-2.2951295) q[1];
sx q[1];
rz(0.15165374) q[1];
x q[2];
rz(0.78403715) q[3];
sx q[3];
rz(-1.136354) q[3];
sx q[3];
rz(-1.2617574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1672704) q[2];
sx q[2];
rz(-1.3784626) q[2];
sx q[2];
rz(0.66217011) q[2];
rz(2.3769489) q[3];
sx q[3];
rz(-1.60631) q[3];
sx q[3];
rz(1.8173328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8793256) q[0];
sx q[0];
rz(-0.58795324) q[0];
sx q[0];
rz(1.3978488) q[0];
rz(-2.0715879) q[1];
sx q[1];
rz(-1.1281697) q[1];
sx q[1];
rz(-1.8096583) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4797213) q[0];
sx q[0];
rz(-2.9878231) q[0];
sx q[0];
rz(-1.3698306) q[0];
rz(1.1564163) q[2];
sx q[2];
rz(-1.3480901) q[2];
sx q[2];
rz(-0.10939314) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9402658) q[1];
sx q[1];
rz(-1.801071) q[1];
sx q[1];
rz(0.23576945) q[1];
rz(-pi) q[2];
rz(0.58140786) q[3];
sx q[3];
rz(-2.1801342) q[3];
sx q[3];
rz(1.4319624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8269044) q[2];
sx q[2];
rz(-1.1004227) q[2];
sx q[2];
rz(2.9617214) q[2];
rz(-1.7449069) q[3];
sx q[3];
rz(-1.4254009) q[3];
sx q[3];
rz(0.21656187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77072813) q[0];
sx q[0];
rz(-0.48038078) q[0];
sx q[0];
rz(0.90850716) q[0];
rz(-2.6188359) q[1];
sx q[1];
rz(-2.0261384) q[1];
sx q[1];
rz(-1.1631858) q[1];
rz(2.5066067) q[2];
sx q[2];
rz(-1.9188835) q[2];
sx q[2];
rz(-1.9139482) q[2];
rz(1.2046075) q[3];
sx q[3];
rz(-1.1259176) q[3];
sx q[3];
rz(2.068145) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
