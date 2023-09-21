OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8945583) q[0];
sx q[0];
rz(-1.2556827) q[0];
sx q[0];
rz(-0.32796252) q[0];
rz(-0.2344996) q[1];
sx q[1];
rz(3.3426715) q[1];
sx q[1];
rz(9.3333416) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8592005) q[0];
sx q[0];
rz(-2.1452367) q[0];
sx q[0];
rz(-2.0496856) q[0];
rz(-pi) q[1];
rz(2.6684127) q[2];
sx q[2];
rz(-0.28943974) q[2];
sx q[2];
rz(-2.6602886) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5603197) q[1];
sx q[1];
rz(-0.25174192) q[1];
sx q[1];
rz(1.9074349) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8880635) q[3];
sx q[3];
rz(-2.5228365) q[3];
sx q[3];
rz(1.3878824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.66449195) q[2];
sx q[2];
rz(-0.97350073) q[2];
sx q[2];
rz(-2.0155902) q[2];
rz(-2.8664355) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1317516) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(-0.69357187) q[0];
rz(1.0961078) q[1];
sx q[1];
rz(-0.98384905) q[1];
sx q[1];
rz(0.19031659) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7305671) q[0];
sx q[0];
rz(-1.9733631) q[0];
sx q[0];
rz(2.9257724) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6066577) q[2];
sx q[2];
rz(-1.9133647) q[2];
sx q[2];
rz(1.5869706) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.39365921) q[1];
sx q[1];
rz(-0.84888443) q[1];
sx q[1];
rz(3.1118803) q[1];
x q[2];
rz(2.6211561) q[3];
sx q[3];
rz(-0.79870975) q[3];
sx q[3];
rz(-0.48898104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.50743121) q[2];
sx q[2];
rz(-0.71175152) q[2];
sx q[2];
rz(-1.921839) q[2];
rz(-2.9988585) q[3];
sx q[3];
rz(-1.0120564) q[3];
sx q[3];
rz(0.90863168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.398657) q[0];
sx q[0];
rz(-2.0516899) q[0];
sx q[0];
rz(-2.3213342) q[0];
rz(2.8495158) q[1];
sx q[1];
rz(-2.0675817) q[1];
sx q[1];
rz(1.2480199) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4688063) q[0];
sx q[0];
rz(-0.60657036) q[0];
sx q[0];
rz(0.16288217) q[0];
x q[1];
rz(-0.29675656) q[2];
sx q[2];
rz(-2.9251758) q[2];
sx q[2];
rz(1.4518472) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6229912) q[1];
sx q[1];
rz(-0.90497436) q[1];
sx q[1];
rz(1.4856505) q[1];
rz(-pi) q[2];
rz(-2.2218024) q[3];
sx q[3];
rz(-1.8120822) q[3];
sx q[3];
rz(-1.9581219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8212006) q[2];
sx q[2];
rz(-2.0911066) q[2];
sx q[2];
rz(1.2134264) q[2];
rz(2.976867) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40959013) q[0];
sx q[0];
rz(-1.859917) q[0];
sx q[0];
rz(0.18606342) q[0];
rz(0.20446725) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(-1.8444555) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0766749) q[0];
sx q[0];
rz(-0.13741048) q[0];
sx q[0];
rz(-1.8924367) q[0];
rz(-pi) q[1];
rz(2.0879891) q[2];
sx q[2];
rz(-2.811811) q[2];
sx q[2];
rz(-1.5151378) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.101673) q[1];
sx q[1];
rz(-0.97517255) q[1];
sx q[1];
rz(2.4458829) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6468871) q[3];
sx q[3];
rz(-2.9615059) q[3];
sx q[3];
rz(2.0538581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90594784) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(-2.2223991) q[2];
rz(-0.32133189) q[3];
sx q[3];
rz(-2.0733757) q[3];
sx q[3];
rz(-1.2478158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677251) q[0];
sx q[0];
rz(-2.6350832) q[0];
sx q[0];
rz(-2.2633973) q[0];
rz(1.8163266) q[1];
sx q[1];
rz(-1.4122496) q[1];
sx q[1];
rz(-1.426288) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25709897) q[0];
sx q[0];
rz(-1.4887267) q[0];
sx q[0];
rz(-0.9135855) q[0];
x q[1];
rz(2.9039188) q[2];
sx q[2];
rz(-1.875068) q[2];
sx q[2];
rz(-1.8397457) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6225699) q[1];
sx q[1];
rz(-1.4624603) q[1];
sx q[1];
rz(-2.3035754) q[1];
rz(-pi) q[2];
rz(0.77107314) q[3];
sx q[3];
rz(-1.3372476) q[3];
sx q[3];
rz(1.1119103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2720126) q[2];
sx q[2];
rz(-11*pi/13) q[2];
sx q[2];
rz(-0.041794725) q[2];
rz(-3.0801008) q[3];
sx q[3];
rz(-2.0525335) q[3];
sx q[3];
rz(0.4367691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7866661) q[0];
sx q[0];
rz(-3.0556398) q[0];
sx q[0];
rz(-2.1110995) q[0];
rz(-0.73973918) q[1];
sx q[1];
rz(-1.5286427) q[1];
sx q[1];
rz(0.57156634) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33486903) q[0];
sx q[0];
rz(-2.3248701) q[0];
sx q[0];
rz(2.4094894) q[0];
rz(-pi) q[1];
rz(-1.7320485) q[2];
sx q[2];
rz(-0.67220062) q[2];
sx q[2];
rz(1.4816928) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.188365) q[1];
sx q[1];
rz(-1.1168915) q[1];
sx q[1];
rz(-0.26785775) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3011175) q[3];
sx q[3];
rz(-2.2395036) q[3];
sx q[3];
rz(-0.83276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.17343865) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(2.5773933) q[2];
rz(-0.12600222) q[3];
sx q[3];
rz(-1.4583476) q[3];
sx q[3];
rz(2.8469767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0512222) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(2.4293161) q[0];
rz(0.5258711) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(2.4760822) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77008477) q[0];
sx q[0];
rz(-1.5697644) q[0];
sx q[0];
rz(-1.3285711) q[0];
rz(-pi) q[1];
x q[1];
rz(2.268157) q[2];
sx q[2];
rz(-1.4566112) q[2];
sx q[2];
rz(0.21542491) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4588786) q[1];
sx q[1];
rz(-2.6280118) q[1];
sx q[1];
rz(2.0290124) q[1];
rz(1.0379167) q[3];
sx q[3];
rz(-0.27927342) q[3];
sx q[3];
rz(2.6168407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9843288) q[2];
sx q[2];
rz(-0.66005808) q[2];
sx q[2];
rz(1.7187913) q[2];
rz(3.026399) q[3];
sx q[3];
rz(-0.51518232) q[3];
sx q[3];
rz(2.9490024) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0916864) q[0];
sx q[0];
rz(-2.005907) q[0];
sx q[0];
rz(-0.19083047) q[0];
rz(-2.514839) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(-0.33871067) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9324547) q[0];
sx q[0];
rz(-1.7212241) q[0];
sx q[0];
rz(-1.223279) q[0];
rz(1.3391101) q[2];
sx q[2];
rz(-2.0121687) q[2];
sx q[2];
rz(-2.1053134) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.45704493) q[1];
sx q[1];
rz(-1.0984048) q[1];
sx q[1];
rz(3.0169675) q[1];
rz(-pi) q[2];
rz(-2.8645105) q[3];
sx q[3];
rz(-1.2556561) q[3];
sx q[3];
rz(1.806123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64615858) q[2];
sx q[2];
rz(-2.50714) q[2];
sx q[2];
rz(2.4411566) q[2];
rz(-2.2436079) q[3];
sx q[3];
rz(-1.8549517) q[3];
sx q[3];
rz(0.17351304) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53196466) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(-0.65834808) q[0];
rz(-2.530653) q[1];
sx q[1];
rz(-1.2555723) q[1];
sx q[1];
rz(-0.13959612) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9717279) q[0];
sx q[0];
rz(-1.5908049) q[0];
sx q[0];
rz(1.6318897) q[0];
x q[1];
rz(3.082824) q[2];
sx q[2];
rz(-0.81548703) q[2];
sx q[2];
rz(-1.8975443) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4715251) q[1];
sx q[1];
rz(-1.0499665) q[1];
sx q[1];
rz(-1.4937558) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0890907) q[3];
sx q[3];
rz(-2.2146261) q[3];
sx q[3];
rz(-0.98464636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5247941) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(2.810478) q[2];
rz(-0.75774276) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(3.0537135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963294) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(2.9560126) q[0];
rz(1.0962076) q[1];
sx q[1];
rz(-0.21462333) q[1];
sx q[1];
rz(1.4846444) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.164924) q[0];
sx q[0];
rz(-0.52545588) q[0];
sx q[0];
rz(2.1303961) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8469641) q[2];
sx q[2];
rz(-1.3535) q[2];
sx q[2];
rz(-1.5664139) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85877307) q[1];
sx q[1];
rz(-0.38837896) q[1];
sx q[1];
rz(-1.8358843) q[1];
rz(-pi) q[2];
rz(-0.86587571) q[3];
sx q[3];
rz(-1.5670334) q[3];
sx q[3];
rz(0.83434425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2075656) q[2];
sx q[2];
rz(-2.2521314) q[2];
sx q[2];
rz(-0.55220848) q[2];
rz(-0.77783716) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(-2.623693) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26375297) q[0];
sx q[0];
rz(-1.629566) q[0];
sx q[0];
rz(-1.4019479) q[0];
rz(2.9539625) q[1];
sx q[1];
rz(-1.3627121) q[1];
sx q[1];
rz(-0.77308853) q[1];
rz(-1.7170231) q[2];
sx q[2];
rz(-1.29371) q[2];
sx q[2];
rz(0.84809662) q[2];
rz(-2.256176) q[3];
sx q[3];
rz(-2.1366742) q[3];
sx q[3];
rz(-2.4921806) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];