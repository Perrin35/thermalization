OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0319808) q[0];
sx q[0];
rz(4.1469753) q[0];
sx q[0];
rz(8.4561705) q[0];
rz(-2.3321416) q[1];
sx q[1];
rz(4.6757698) q[1];
sx q[1];
rz(9.4546131) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.580509) q[0];
sx q[0];
rz(-0.59701118) q[0];
sx q[0];
rz(2.719814) q[0];
rz(-pi) q[1];
rz(2.3230629) q[2];
sx q[2];
rz(-1.6224183) q[2];
sx q[2];
rz(0.35721401) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3181495) q[1];
sx q[1];
rz(-1.9572721) q[1];
sx q[1];
rz(2.3837368) q[1];
rz(-pi) q[2];
rz(-2.5621595) q[3];
sx q[3];
rz(-2.979485) q[3];
sx q[3];
rz(2.1539089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5030824) q[2];
sx q[2];
rz(-0.57380399) q[2];
sx q[2];
rz(-0.32162515) q[2];
rz(-2.5173748) q[3];
sx q[3];
rz(-1.3892684) q[3];
sx q[3];
rz(1.6577087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6931848) q[0];
sx q[0];
rz(-0.25633651) q[0];
sx q[0];
rz(-0.93628991) q[0];
rz(1.869465) q[1];
sx q[1];
rz(-1.5446168) q[1];
sx q[1];
rz(1.7535271) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6427338) q[0];
sx q[0];
rz(-1.6279477) q[0];
sx q[0];
rz(1.3512035) q[0];
rz(-pi) q[1];
rz(0.77735591) q[2];
sx q[2];
rz(-1.1749975) q[2];
sx q[2];
rz(1.1071902) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9656877) q[1];
sx q[1];
rz(-2.2908194) q[1];
sx q[1];
rz(-2.4026924) q[1];
rz(-pi) q[2];
rz(0.21556385) q[3];
sx q[3];
rz(-1.8290534) q[3];
sx q[3];
rz(0.35999957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0442514) q[2];
sx q[2];
rz(-2.1772431) q[2];
sx q[2];
rz(-0.67330366) q[2];
rz(-1.3168969) q[3];
sx q[3];
rz(-1.8484867) q[3];
sx q[3];
rz(2.3015658) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8116542) q[0];
sx q[0];
rz(-1.5978156) q[0];
sx q[0];
rz(3.0918308) q[0];
rz(2.9352169) q[1];
sx q[1];
rz(-1.754909) q[1];
sx q[1];
rz(-1.6516364) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39078907) q[0];
sx q[0];
rz(-2.514578) q[0];
sx q[0];
rz(0.20558447) q[0];
rz(-pi) q[1];
rz(-1.1300722) q[2];
sx q[2];
rz(-2.080938) q[2];
sx q[2];
rz(-0.59512344) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.083410376) q[1];
sx q[1];
rz(-1.8334926) q[1];
sx q[1];
rz(-1.865126) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95155119) q[3];
sx q[3];
rz(-2.5971322) q[3];
sx q[3];
rz(-1.3109372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9516307) q[2];
sx q[2];
rz(-2.5962574) q[2];
sx q[2];
rz(3.087431) q[2];
rz(-1.9870029) q[3];
sx q[3];
rz(-1.2879939) q[3];
sx q[3];
rz(-2.0258928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7390249) q[0];
sx q[0];
rz(-1.4250647) q[0];
sx q[0];
rz(-2.2836852) q[0];
rz(2.4671381) q[1];
sx q[1];
rz(-0.94991389) q[1];
sx q[1];
rz(1.3309006) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20786634) q[0];
sx q[0];
rz(-0.27369341) q[0];
sx q[0];
rz(-2.6143608) q[0];
rz(-2.0728378) q[2];
sx q[2];
rz(-0.34345657) q[2];
sx q[2];
rz(-2.1844027) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7965296) q[1];
sx q[1];
rz(-1.3994675) q[1];
sx q[1];
rz(-1.0061128) q[1];
rz(1.1515297) q[3];
sx q[3];
rz(-1.5738838) q[3];
sx q[3];
rz(-3.1379238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7669507) q[2];
sx q[2];
rz(-2.9810814) q[2];
sx q[2];
rz(3.0373419) q[2];
rz(-2.5637964) q[3];
sx q[3];
rz(-1.0775074) q[3];
sx q[3];
rz(-0.92208636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0084956) q[0];
sx q[0];
rz(-1.47559) q[0];
sx q[0];
rz(0.9444899) q[0];
rz(-1.3635483) q[1];
sx q[1];
rz(-1.3282158) q[1];
sx q[1];
rz(1.69467) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5529842) q[0];
sx q[0];
rz(-0.7352587) q[0];
sx q[0];
rz(1.2862434) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1051272) q[2];
sx q[2];
rz(-1.8210095) q[2];
sx q[2];
rz(-1.6584876) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8268298) q[1];
sx q[1];
rz(-2.3896425) q[1];
sx q[1];
rz(0.58677267) q[1];
rz(-0.99926104) q[3];
sx q[3];
rz(-2.0700698) q[3];
sx q[3];
rz(-1.5443791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.038534433) q[2];
sx q[2];
rz(-2.6699622) q[2];
sx q[2];
rz(2.4209957) q[2];
rz(-0.12901148) q[3];
sx q[3];
rz(-1.7804264) q[3];
sx q[3];
rz(-1.2989929) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98220372) q[0];
sx q[0];
rz(-2.144618) q[0];
sx q[0];
rz(-2.5250315) q[0];
rz(-2.189134) q[1];
sx q[1];
rz(-1.1943694) q[1];
sx q[1];
rz(1.2430826) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6641419) q[0];
sx q[0];
rz(-1.7074025) q[0];
sx q[0];
rz(1.8325498) q[0];
rz(-0.41105777) q[2];
sx q[2];
rz(-0.90663821) q[2];
sx q[2];
rz(2.1020087) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3042708) q[1];
sx q[1];
rz(-2.4189286) q[1];
sx q[1];
rz(3.0501306) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9493306) q[3];
sx q[3];
rz(-1.6311181) q[3];
sx q[3];
rz(2.5669135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9055966) q[2];
sx q[2];
rz(-2.8541028) q[2];
sx q[2];
rz(3.0671885) q[2];
rz(-1.0807886) q[3];
sx q[3];
rz(-1.4932884) q[3];
sx q[3];
rz(-2.1006445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3022795) q[0];
sx q[0];
rz(-1.1458719) q[0];
sx q[0];
rz(0.062051274) q[0];
rz(-1.7265559) q[1];
sx q[1];
rz(-2.3686385) q[1];
sx q[1];
rz(0.089990377) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9190529) q[0];
sx q[0];
rz(-2.3947859) q[0];
sx q[0];
rz(1.1928305) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7715459) q[2];
sx q[2];
rz(-2.4436969) q[2];
sx q[2];
rz(3.0332886) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18971009) q[1];
sx q[1];
rz(-1.9196577) q[1];
sx q[1];
rz(-2.2079218) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3004933) q[3];
sx q[3];
rz(-0.78783082) q[3];
sx q[3];
rz(-2.6177518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.0061079582) q[2];
sx q[2];
rz(-1.2697479) q[2];
sx q[2];
rz(-0.6832068) q[2];
rz(-1.3225383) q[3];
sx q[3];
rz(-2.1893978) q[3];
sx q[3];
rz(-1.3873842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(1.4676062) q[0];
sx q[0];
rz(-1.645393) q[0];
sx q[0];
rz(2.6737387) q[0];
rz(-2.4521008) q[1];
sx q[1];
rz(-0.71527022) q[1];
sx q[1];
rz(2.668344) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1263862) q[0];
sx q[0];
rz(-1.6412982) q[0];
sx q[0];
rz(1.3944993) q[0];
rz(-pi) q[1];
rz(2.8125561) q[2];
sx q[2];
rz(-0.46571443) q[2];
sx q[2];
rz(-2.8570115) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2040904) q[1];
sx q[1];
rz(-0.023462208) q[1];
sx q[1];
rz(0.84659441) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7353294) q[3];
sx q[3];
rz(-2.7711754) q[3];
sx q[3];
rz(0.35740023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1271992) q[2];
sx q[2];
rz(-1.3213804) q[2];
sx q[2];
rz(-0.42352208) q[2];
rz(-1.8262156) q[3];
sx q[3];
rz(-0.49153057) q[3];
sx q[3];
rz(1.2929085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92726436) q[0];
sx q[0];
rz(-2.1926227) q[0];
sx q[0];
rz(1.8073136) q[0];
rz(-1.8233874) q[1];
sx q[1];
rz(-2.2315836) q[1];
sx q[1];
rz(-0.013890161) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62187884) q[0];
sx q[0];
rz(-1.5962999) q[0];
sx q[0];
rz(-1.9344781) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35160323) q[2];
sx q[2];
rz(-1.7187864) q[2];
sx q[2];
rz(-0.44662145) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.23194961) q[1];
sx q[1];
rz(-1.0000537) q[1];
sx q[1];
rz(0.29517031) q[1];
rz(0.84133103) q[3];
sx q[3];
rz(-1.4784665) q[3];
sx q[3];
rz(-2.009237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5615329) q[2];
sx q[2];
rz(-1.1416898) q[2];
sx q[2];
rz(-2.0884936) q[2];
rz(2.2140391) q[3];
sx q[3];
rz(-1.901123) q[3];
sx q[3];
rz(-2.659306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7841566) q[0];
sx q[0];
rz(-0.3967537) q[0];
sx q[0];
rz(0.60047737) q[0];
rz(-0.26957574) q[1];
sx q[1];
rz(-0.69200626) q[1];
sx q[1];
rz(-1.3788266) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10256448) q[0];
sx q[0];
rz(-1.9223273) q[0];
sx q[0];
rz(0.28798703) q[0];
rz(-0.68280934) q[2];
sx q[2];
rz(-0.46418846) q[2];
sx q[2];
rz(0.95303553) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.87390806) q[1];
sx q[1];
rz(-0.46894858) q[1];
sx q[1];
rz(-2.7280366) q[1];
rz(-0.58568256) q[3];
sx q[3];
rz(-2.437311) q[3];
sx q[3];
rz(1.3552988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.025734162) q[2];
sx q[2];
rz(-0.33418843) q[2];
sx q[2];
rz(0.96246976) q[2];
rz(2.1879503) q[3];
sx q[3];
rz(-1.0735268) q[3];
sx q[3];
rz(1.7474705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7643323) q[0];
sx q[0];
rz(-0.88484103) q[0];
sx q[0];
rz(-1.9405889) q[0];
rz(-1.7878905) q[1];
sx q[1];
rz(-1.9722912) q[1];
sx q[1];
rz(-0.73696662) q[1];
rz(-2.3804661) q[2];
sx q[2];
rz(-0.82999753) q[2];
sx q[2];
rz(0.33090811) q[2];
rz(-0.50696452) q[3];
sx q[3];
rz(-1.4486113) q[3];
sx q[3];
rz(2.2511528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
