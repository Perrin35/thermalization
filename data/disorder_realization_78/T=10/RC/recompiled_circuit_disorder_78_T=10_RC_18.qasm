OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.55943638) q[0];
sx q[0];
rz(-2.5506033) q[0];
sx q[0];
rz(-0.58340573) q[0];
rz(2.9572339) q[1];
sx q[1];
rz(-0.98362041) q[1];
sx q[1];
rz(2.2489927) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8164506) q[0];
sx q[0];
rz(-2.4976375) q[0];
sx q[0];
rz(-2.08026) q[0];
x q[1];
rz(-2.5298654) q[2];
sx q[2];
rz(-0.76843843) q[2];
sx q[2];
rz(-0.4380463) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0639227) q[1];
sx q[1];
rz(-2.3706672) q[1];
sx q[1];
rz(2.0248807) q[1];
rz(-pi) q[2];
rz(-0.75099545) q[3];
sx q[3];
rz(-1.6023811) q[3];
sx q[3];
rz(0.75814523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9154174) q[2];
sx q[2];
rz(-1.6025275) q[2];
sx q[2];
rz(-1.1606476) q[2];
rz(-0.21696572) q[3];
sx q[3];
rz(-2.6187077) q[3];
sx q[3];
rz(1.0552361) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0579257) q[0];
sx q[0];
rz(-0.96494976) q[0];
sx q[0];
rz(0.57587409) q[0];
rz(1.8946164) q[1];
sx q[1];
rz(-1.8449102) q[1];
sx q[1];
rz(-1.974568) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28536797) q[0];
sx q[0];
rz(-1.8282187) q[0];
sx q[0];
rz(3.0811937) q[0];
x q[1];
rz(-2.1115233) q[2];
sx q[2];
rz(-1.3582555) q[2];
sx q[2];
rz(-2.9589047) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6845219) q[1];
sx q[1];
rz(-1.6324537) q[1];
sx q[1];
rz(0.80656273) q[1];
rz(-2.4715273) q[3];
sx q[3];
rz(-1.1747922) q[3];
sx q[3];
rz(-2.0585287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1077659) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(-0.21437422) q[2];
rz(-0.073444627) q[3];
sx q[3];
rz(-2.6918604) q[3];
sx q[3];
rz(0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4784933) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(-1.9146772) q[0];
rz(-2.7413209) q[1];
sx q[1];
rz(-1.2534671) q[1];
sx q[1];
rz(1.0148369) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0815711) q[0];
sx q[0];
rz(-2.2131753) q[0];
sx q[0];
rz(-1.6409671) q[0];
rz(-1.9760518) q[2];
sx q[2];
rz(-1.793867) q[2];
sx q[2];
rz(-2.9802393) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.866232) q[1];
sx q[1];
rz(-0.97517698) q[1];
sx q[1];
rz(-3.0443405) q[1];
rz(-0.26168163) q[3];
sx q[3];
rz(-2.2398584) q[3];
sx q[3];
rz(-2.3467968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0456475) q[2];
sx q[2];
rz(-2.084705) q[2];
sx q[2];
rz(-2.2423559) q[2];
rz(-2.4441161) q[3];
sx q[3];
rz(-1.2858425) q[3];
sx q[3];
rz(0.60788679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65524453) q[0];
sx q[0];
rz(-1.0826033) q[0];
sx q[0];
rz(0.43193257) q[0];
rz(2.5090384) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(-2.5057709) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.302127) q[0];
sx q[0];
rz(-1.7502022) q[0];
sx q[0];
rz(2.0749712) q[0];
rz(-pi) q[1];
rz(-2.2976774) q[2];
sx q[2];
rz(-1.6944431) q[2];
sx q[2];
rz(0.67019776) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.35866657) q[1];
sx q[1];
rz(-1.3861321) q[1];
sx q[1];
rz(-0.77628805) q[1];
rz(-pi) q[2];
rz(-0.28663978) q[3];
sx q[3];
rz(-1.7754103) q[3];
sx q[3];
rz(2.3972547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2137961) q[2];
sx q[2];
rz(-2.9979604) q[2];
sx q[2];
rz(-2.4528743) q[2];
rz(-2.8074746) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(0.14373246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9976945) q[0];
sx q[0];
rz(-0.69902885) q[0];
sx q[0];
rz(2.0671663) q[0];
rz(0.74514666) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(0.27854663) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6644088) q[0];
sx q[0];
rz(-1.9946949) q[0];
sx q[0];
rz(-2.3197078) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61168806) q[2];
sx q[2];
rz(-1.4809161) q[2];
sx q[2];
rz(2.9181366) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.46036938) q[1];
sx q[1];
rz(-1.9378098) q[1];
sx q[1];
rz(0.91194921) q[1];
x q[2];
rz(2.3260818) q[3];
sx q[3];
rz(-1.1482571) q[3];
sx q[3];
rz(0.58327196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4429861) q[2];
sx q[2];
rz(-1.3977945) q[2];
sx q[2];
rz(2.8473575) q[2];
rz(0.081929835) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(0.036227139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1086403) q[0];
sx q[0];
rz(-0.8529129) q[0];
sx q[0];
rz(0.0090573514) q[0];
rz(-2.5065705) q[1];
sx q[1];
rz(-2.451684) q[1];
sx q[1];
rz(3.0335398) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8979643) q[0];
sx q[0];
rz(-1.122323) q[0];
sx q[0];
rz(1.2178221) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4322386) q[2];
sx q[2];
rz(-1.1922622) q[2];
sx q[2];
rz(-2.0356503) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3012645) q[1];
sx q[1];
rz(-0.64126188) q[1];
sx q[1];
rz(-0.68323369) q[1];
x q[2];
rz(0.26515682) q[3];
sx q[3];
rz(-1.0401871) q[3];
sx q[3];
rz(-1.6924072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99284995) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(1.8704869) q[2];
rz(-3.0631915) q[3];
sx q[3];
rz(-1.4784808) q[3];
sx q[3];
rz(2.0097849) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8909797) q[0];
sx q[0];
rz(-1.5761292) q[0];
sx q[0];
rz(-2.4839731) q[0];
rz(-1.3972067) q[1];
sx q[1];
rz(-1.1746635) q[1];
sx q[1];
rz(-2.2479642) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.87003) q[0];
sx q[0];
rz(-1.8989519) q[0];
sx q[0];
rz(-2.4782655) q[0];
x q[1];
rz(1.0412752) q[2];
sx q[2];
rz(-2.6011701) q[2];
sx q[2];
rz(-1.4735917) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.212008) q[1];
sx q[1];
rz(-1.3699023) q[1];
sx q[1];
rz(0.68540539) q[1];
rz(3.0393533) q[3];
sx q[3];
rz(-2.3792301) q[3];
sx q[3];
rz(-1.8765212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.39067337) q[2];
sx q[2];
rz(-2.1942287) q[2];
sx q[2];
rz(2.612109) q[2];
rz(2.6654065) q[3];
sx q[3];
rz(-1.809285) q[3];
sx q[3];
rz(0.34255323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7628409) q[0];
sx q[0];
rz(-1.6289926) q[0];
sx q[0];
rz(-1.09028) q[0];
rz(3.0293363) q[1];
sx q[1];
rz(-1.1039762) q[1];
sx q[1];
rz(-1.9876678) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.740828) q[0];
sx q[0];
rz(-0.65069288) q[0];
sx q[0];
rz(-0.056218938) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9032193) q[2];
sx q[2];
rz(-0.56406883) q[2];
sx q[2];
rz(0.43018815) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.47626074) q[1];
sx q[1];
rz(-1.5090824) q[1];
sx q[1];
rz(-2.3803821) q[1];
rz(-pi) q[2];
rz(2.7525547) q[3];
sx q[3];
rz(-2.2044047) q[3];
sx q[3];
rz(0.47215677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.551679) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(0.38044688) q[2];
rz(2.0137265) q[3];
sx q[3];
rz(-1.7815855) q[3];
sx q[3];
rz(-1.5650704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34148759) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(-2.5019116) q[0];
rz(-1.2387964) q[1];
sx q[1];
rz(-1.4150554) q[1];
sx q[1];
rz(1.170084) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026924883) q[0];
sx q[0];
rz(-1.9244734) q[0];
sx q[0];
rz(-0.35004079) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3652472) q[2];
sx q[2];
rz(-1.3753969) q[2];
sx q[2];
rz(2.5459144) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1994836) q[1];
sx q[1];
rz(-2.0523199) q[1];
sx q[1];
rz(1.3660618) q[1];
rz(0.89948489) q[3];
sx q[3];
rz(-1.6791108) q[3];
sx q[3];
rz(1.2182106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3141979) q[2];
sx q[2];
rz(-0.65076995) q[2];
sx q[2];
rz(2.7588552) q[2];
rz(2.2132204) q[3];
sx q[3];
rz(-1.9675156) q[3];
sx q[3];
rz(-2.476957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2255573) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(0.25892192) q[0];
rz(-0.71031538) q[1];
sx q[1];
rz(-2.0878891) q[1];
sx q[1];
rz(2.6616667) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9057248) q[0];
sx q[0];
rz(-0.94434443) q[0];
sx q[0];
rz(-2.7263374) q[0];
rz(-2.9051022) q[2];
sx q[2];
rz(-0.91954008) q[2];
sx q[2];
rz(1.0054393) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.48206115) q[1];
sx q[1];
rz(-1.9226942) q[1];
sx q[1];
rz(0.54606502) q[1];
rz(-pi) q[2];
rz(-2.6142526) q[3];
sx q[3];
rz(-2.1747327) q[3];
sx q[3];
rz(-2.7690943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.84247983) q[2];
sx q[2];
rz(-0.92876902) q[2];
sx q[2];
rz(1.3170362) q[2];
rz(-1.8995829) q[3];
sx q[3];
rz(-2.1879991) q[3];
sx q[3];
rz(0.73808134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.15923545) q[0];
sx q[0];
rz(-0.032082162) q[0];
sx q[0];
rz(-1.4627009) q[0];
rz(2.1622529) q[1];
sx q[1];
rz(-2.0420488) q[1];
sx q[1];
rz(2.2534823) q[1];
rz(0.48158823) q[2];
sx q[2];
rz(-1.4143741) q[2];
sx q[2];
rz(1.0798567) q[2];
rz(0.76673037) q[3];
sx q[3];
rz(-2.7646716) q[3];
sx q[3];
rz(-0.9203831) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
