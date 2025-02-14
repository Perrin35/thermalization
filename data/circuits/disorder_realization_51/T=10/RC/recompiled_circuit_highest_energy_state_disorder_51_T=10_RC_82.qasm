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
rz(0.32938862) q[0];
sx q[0];
rz(-2.5100799) q[0];
sx q[0];
rz(0.00087498571) q[0];
rz(2.5007091) q[1];
sx q[1];
rz(-2.1407318) q[1];
sx q[1];
rz(-0.34520087) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5986818) q[0];
sx q[0];
rz(-0.094176725) q[0];
sx q[0];
rz(-1.3046632) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7852374) q[2];
sx q[2];
rz(-0.40268597) q[2];
sx q[2];
rz(0.25549437) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7789014) q[1];
sx q[1];
rz(-2.3786015) q[1];
sx q[1];
rz(-2.1982026) q[1];
x q[2];
rz(-3.1253417) q[3];
sx q[3];
rz(-1.4964087) q[3];
sx q[3];
rz(0.95916884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.22902809) q[2];
sx q[2];
rz(-1.8077069) q[2];
sx q[2];
rz(-2.933617) q[2];
rz(2.8424756) q[3];
sx q[3];
rz(-0.57204539) q[3];
sx q[3];
rz(-2.0007029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3308554) q[0];
sx q[0];
rz(-2.9006697) q[0];
sx q[0];
rz(2.148707) q[0];
rz(-1.8244686) q[1];
sx q[1];
rz(-2.8254852) q[1];
sx q[1];
rz(0.77450007) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74906384) q[0];
sx q[0];
rz(-1.3925902) q[0];
sx q[0];
rz(0.011179608) q[0];
x q[1];
rz(0.80747898) q[2];
sx q[2];
rz(-1.9614855) q[2];
sx q[2];
rz(1.1451858) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5238873) q[1];
sx q[1];
rz(-1.7951709) q[1];
sx q[1];
rz(-2.7588506) q[1];
rz(-2.0404313) q[3];
sx q[3];
rz(-1.7723933) q[3];
sx q[3];
rz(-0.93702173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.896686) q[2];
sx q[2];
rz(-1.8550355) q[2];
sx q[2];
rz(-1.0150821) q[2];
rz(-0.24909881) q[3];
sx q[3];
rz(-2.2738012) q[3];
sx q[3];
rz(0.36756137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86431137) q[0];
sx q[0];
rz(-0.69121498) q[0];
sx q[0];
rz(2.9840898) q[0];
rz(-2.1306254) q[1];
sx q[1];
rz(-0.33282655) q[1];
sx q[1];
rz(-0.13883042) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.39108) q[0];
sx q[0];
rz(-1.8939928) q[0];
sx q[0];
rz(1.1884407) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9857668) q[2];
sx q[2];
rz(-0.89902821) q[2];
sx q[2];
rz(0.49723724) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.638354) q[1];
sx q[1];
rz(-1.279083) q[1];
sx q[1];
rz(2.9364763) q[1];
x q[2];
rz(-1.3957455) q[3];
sx q[3];
rz(-1.2915466) q[3];
sx q[3];
rz(1.4917013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.46131721) q[2];
sx q[2];
rz(-0.91247827) q[2];
sx q[2];
rz(-2.7834564) q[2];
rz(-0.67874587) q[3];
sx q[3];
rz(-2.1923784) q[3];
sx q[3];
rz(-2.1024735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.045227483) q[0];
sx q[0];
rz(-0.97310936) q[0];
sx q[0];
rz(-2.8367693) q[0];
rz(-2.7335956) q[1];
sx q[1];
rz(-1.4381189) q[1];
sx q[1];
rz(0.92794424) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5261032) q[0];
sx q[0];
rz(-2.9080221) q[0];
sx q[0];
rz(-2.156267) q[0];
x q[1];
rz(2.1963652) q[2];
sx q[2];
rz(-2.0778227) q[2];
sx q[2];
rz(-0.27238174) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2603307) q[1];
sx q[1];
rz(-1.7219543) q[1];
sx q[1];
rz(2.4511454) q[1];
rz(-2.0870861) q[3];
sx q[3];
rz(-2.1355503) q[3];
sx q[3];
rz(1.3206583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9503595) q[2];
sx q[2];
rz(-0.57098907) q[2];
sx q[2];
rz(2.9534269) q[2];
rz(-1.865271) q[3];
sx q[3];
rz(-1.3151582) q[3];
sx q[3];
rz(1.4307384) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6998049) q[0];
sx q[0];
rz(-2.7975174) q[0];
sx q[0];
rz(0.019388327) q[0];
rz(-1.060932) q[1];
sx q[1];
rz(-1.7273936) q[1];
sx q[1];
rz(1.2394989) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083154924) q[0];
sx q[0];
rz(-0.64955901) q[0];
sx q[0];
rz(1.0970647) q[0];
rz(-pi) q[1];
rz(-2.4024348) q[2];
sx q[2];
rz(-2.0817167) q[2];
sx q[2];
rz(0.3581274) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8722788) q[1];
sx q[1];
rz(-1.8511103) q[1];
sx q[1];
rz(2.6107236) q[1];
x q[2];
rz(2.7321759) q[3];
sx q[3];
rz(-2.0338661) q[3];
sx q[3];
rz(-0.45209979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1218607) q[2];
sx q[2];
rz(-1.3132361) q[2];
sx q[2];
rz(-1.8149553) q[2];
rz(0.2462247) q[3];
sx q[3];
rz(-1.1028057) q[3];
sx q[3];
rz(-2.2897913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5354079) q[0];
sx q[0];
rz(-0.98537213) q[0];
sx q[0];
rz(2.4024409) q[0];
rz(-0.34573653) q[1];
sx q[1];
rz(-1.3290936) q[1];
sx q[1];
rz(-2.2703222) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6046048) q[0];
sx q[0];
rz(-2.3571627) q[0];
sx q[0];
rz(0.098172763) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75052336) q[2];
sx q[2];
rz(-2.0156246) q[2];
sx q[2];
rz(1.9409279) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.024684357) q[1];
sx q[1];
rz(-1.4977807) q[1];
sx q[1];
rz(1.1817929) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5302251) q[3];
sx q[3];
rz(-1.1725559) q[3];
sx q[3];
rz(-2.7643725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.31200108) q[2];
sx q[2];
rz(-2.5025554) q[2];
sx q[2];
rz(-2.269022) q[2];
rz(-1.5729337) q[3];
sx q[3];
rz(-2.5376153) q[3];
sx q[3];
rz(-2.2085371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0503814) q[0];
sx q[0];
rz(-0.25930431) q[0];
sx q[0];
rz(0.76164371) q[0];
rz(2.7161982) q[1];
sx q[1];
rz(-2.0014706) q[1];
sx q[1];
rz(-2.2147307) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8445209) q[0];
sx q[0];
rz(-1.335134) q[0];
sx q[0];
rz(-1.336444) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2144187) q[2];
sx q[2];
rz(-2.1774051) q[2];
sx q[2];
rz(-1.7546897) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35173479) q[1];
sx q[1];
rz(-2.7119293) q[1];
sx q[1];
rz(-1.0142782) q[1];
rz(-pi) q[2];
rz(-1.6377444) q[3];
sx q[3];
rz(-2.5173325) q[3];
sx q[3];
rz(0.74893307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0285792) q[2];
sx q[2];
rz(-2.4455652) q[2];
sx q[2];
rz(2.2704303) q[2];
rz(-0.29843676) q[3];
sx q[3];
rz(-1.8961743) q[3];
sx q[3];
rz(-1.8106073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4153862) q[0];
sx q[0];
rz(-2.6415249) q[0];
sx q[0];
rz(-0.4183847) q[0];
rz(-0.2977953) q[1];
sx q[1];
rz(-1.9458385) q[1];
sx q[1];
rz(0.13430886) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9878085) q[0];
sx q[0];
rz(-2.1001108) q[0];
sx q[0];
rz(-0.80940009) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1573345) q[2];
sx q[2];
rz(-0.43817876) q[2];
sx q[2];
rz(2.7587121) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5843552) q[1];
sx q[1];
rz(-1.3882625) q[1];
sx q[1];
rz(1.4339158) q[1];
rz(-pi) q[2];
rz(2.8577096) q[3];
sx q[3];
rz(-1.1177269) q[3];
sx q[3];
rz(2.4456152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.40992752) q[2];
sx q[2];
rz(-1.2890559) q[2];
sx q[2];
rz(0.64201391) q[2];
rz(-1.0632473) q[3];
sx q[3];
rz(-0.65170538) q[3];
sx q[3];
rz(2.1003907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6006271) q[0];
sx q[0];
rz(-0.0890812) q[0];
sx q[0];
rz(-0.11216057) q[0];
rz(-1.8070096) q[1];
sx q[1];
rz(-2.5490675) q[1];
sx q[1];
rz(1.0293915) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4698668) q[0];
sx q[0];
rz(-0.70505667) q[0];
sx q[0];
rz(-1.5211578) q[0];
rz(-pi) q[1];
rz(-1.5310982) q[2];
sx q[2];
rz(-1.8744812) q[2];
sx q[2];
rz(-3.1162709) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7299906) q[1];
sx q[1];
rz(-1.8008324) q[1];
sx q[1];
rz(0.045129808) q[1];
rz(-pi) q[2];
rz(-1.6624032) q[3];
sx q[3];
rz(-1.1044958) q[3];
sx q[3];
rz(2.8331851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7182497) q[2];
sx q[2];
rz(-3.0666879) q[2];
sx q[2];
rz(3.1035799) q[2];
rz(2.5293317) q[3];
sx q[3];
rz(-0.93658787) q[3];
sx q[3];
rz(-0.14122252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.8925979) q[0];
sx q[0];
rz(-1.4709512) q[0];
sx q[0];
rz(-0.41879642) q[0];
rz(1.2576125) q[1];
sx q[1];
rz(-1.6323171) q[1];
sx q[1];
rz(-0.14258252) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9197971) q[0];
sx q[0];
rz(-2.9807973) q[0];
sx q[0];
rz(0.17531403) q[0];
rz(-1.293574) q[2];
sx q[2];
rz(-2.1810594) q[2];
sx q[2];
rz(2.0137613) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6302418) q[1];
sx q[1];
rz(-0.6410743) q[1];
sx q[1];
rz(1.1180693) q[1];
rz(-1.110465) q[3];
sx q[3];
rz(-1.2602196) q[3];
sx q[3];
rz(-0.47720695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7964145) q[2];
sx q[2];
rz(-3.0609481) q[2];
sx q[2];
rz(-0.26819116) q[2];
rz(1.4043407) q[3];
sx q[3];
rz(-1.0920478) q[3];
sx q[3];
rz(1.0442737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1254697) q[0];
sx q[0];
rz(-0.37624993) q[0];
sx q[0];
rz(-2.7952623) q[0];
rz(0.12915962) q[1];
sx q[1];
rz(-1.2551413) q[1];
sx q[1];
rz(-1.7358949) q[1];
rz(2.6665192) q[2];
sx q[2];
rz(-1.8865841) q[2];
sx q[2];
rz(0.24812698) q[2];
rz(-1.3553452) q[3];
sx q[3];
rz(-2.8402495) q[3];
sx q[3];
rz(0.80339669) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
