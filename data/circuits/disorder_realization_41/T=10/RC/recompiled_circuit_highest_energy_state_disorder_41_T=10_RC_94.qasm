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
rz(0.32349411) q[0];
sx q[0];
rz(-0.20322023) q[0];
sx q[0];
rz(-2.8414371) q[0];
rz(2.7872941) q[1];
sx q[1];
rz(-1.2343255) q[1];
sx q[1];
rz(-0.77114463) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010492652) q[0];
sx q[0];
rz(-0.26828921) q[0];
sx q[0];
rz(0.77666755) q[0];
x q[1];
rz(-1.096719) q[2];
sx q[2];
rz(-1.0869622) q[2];
sx q[2];
rz(-0.7686309) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.59471666) q[1];
sx q[1];
rz(-1.1961814) q[1];
sx q[1];
rz(2.0105816) q[1];
rz(-pi) q[2];
rz(1.7158137) q[3];
sx q[3];
rz(-0.43802625) q[3];
sx q[3];
rz(1.6159542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.99042088) q[2];
sx q[2];
rz(-2.6800818) q[2];
sx q[2];
rz(-1.0837466) q[2];
rz(0.57461965) q[3];
sx q[3];
rz(-1.5950404) q[3];
sx q[3];
rz(2.3641724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8271178) q[0];
sx q[0];
rz(-0.62905335) q[0];
sx q[0];
rz(0.3984867) q[0];
rz(1.9972948) q[1];
sx q[1];
rz(-0.60595787) q[1];
sx q[1];
rz(0.95432895) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7643549) q[0];
sx q[0];
rz(-1.5713552) q[0];
sx q[0];
rz(0.0045545983) q[0];
rz(-0.09403566) q[2];
sx q[2];
rz(-2.2263029) q[2];
sx q[2];
rz(1.2588225) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7415175) q[1];
sx q[1];
rz(-2.1360101) q[1];
sx q[1];
rz(-2.8533622) q[1];
rz(-2.0404194) q[3];
sx q[3];
rz(-2.0631115) q[3];
sx q[3];
rz(2.5143169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2362471) q[2];
sx q[2];
rz(-0.50062847) q[2];
sx q[2];
rz(-0.11032571) q[2];
rz(-0.8253544) q[3];
sx q[3];
rz(-2.3767411) q[3];
sx q[3];
rz(-0.74025214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1684619) q[0];
sx q[0];
rz(-0.2178807) q[0];
sx q[0];
rz(-2.2595898) q[0];
rz(1.0285671) q[1];
sx q[1];
rz(-2.5879637) q[1];
sx q[1];
rz(-2.2291768) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0093632657) q[0];
sx q[0];
rz(-2.6840586) q[0];
sx q[0];
rz(1.7347764) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4503612) q[2];
sx q[2];
rz(-0.85589534) q[2];
sx q[2];
rz(-0.038106136) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6710224) q[1];
sx q[1];
rz(-0.67364866) q[1];
sx q[1];
rz(-2.1191052) q[1];
x q[2];
rz(-1.0718143) q[3];
sx q[3];
rz(-2.6423965) q[3];
sx q[3];
rz(-2.5170391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2622111) q[2];
sx q[2];
rz(-2.8935581) q[2];
sx q[2];
rz(2.3233419) q[2];
rz(-2.7815871) q[3];
sx q[3];
rz(-1.2986978) q[3];
sx q[3];
rz(-0.037809614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70571947) q[0];
sx q[0];
rz(-2.191045) q[0];
sx q[0];
rz(1.1608423) q[0];
rz(-2.1707161) q[1];
sx q[1];
rz(-2.2258046) q[1];
sx q[1];
rz(0.11347778) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1176096) q[0];
sx q[0];
rz(-1.9736088) q[0];
sx q[0];
rz(-0.65448032) q[0];
rz(-pi) q[1];
rz(-3.1283896) q[2];
sx q[2];
rz(-1.535386) q[2];
sx q[2];
rz(-1.7865045) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7474688) q[1];
sx q[1];
rz(-1.040084) q[1];
sx q[1];
rz(-2.3783724) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3498126) q[3];
sx q[3];
rz(-0.80884514) q[3];
sx q[3];
rz(0.57989449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0299782) q[2];
sx q[2];
rz(-1.8054211) q[2];
sx q[2];
rz(0.90799904) q[2];
rz(2.8152483) q[3];
sx q[3];
rz(-0.55485266) q[3];
sx q[3];
rz(-2.7197796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57080764) q[0];
sx q[0];
rz(-2.3796005) q[0];
sx q[0];
rz(-2.1245891) q[0];
rz(2.6550338) q[1];
sx q[1];
rz(-2.1162972) q[1];
sx q[1];
rz(-0.09495458) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9236889) q[0];
sx q[0];
rz(-0.14553864) q[0];
sx q[0];
rz(-0.83011629) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3516897) q[2];
sx q[2];
rz(-0.75828248) q[2];
sx q[2];
rz(0.74196434) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.44014318) q[1];
sx q[1];
rz(-2.6924163) q[1];
sx q[1];
rz(-0.67235948) q[1];
x q[2];
rz(-2.8311555) q[3];
sx q[3];
rz(-1.2708029) q[3];
sx q[3];
rz(-2.7315549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.93593705) q[2];
sx q[2];
rz(-2.50694) q[2];
sx q[2];
rz(-0.46979365) q[2];
rz(0.69822407) q[3];
sx q[3];
rz(-1.2263115) q[3];
sx q[3];
rz(-2.8214112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6093589) q[0];
sx q[0];
rz(-0.73901743) q[0];
sx q[0];
rz(-0.21482378) q[0];
rz(0.23669446) q[1];
sx q[1];
rz(-2.5645945) q[1];
sx q[1];
rz(2.7223041) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12690645) q[0];
sx q[0];
rz(-1.4923411) q[0];
sx q[0];
rz(-0.04562745) q[0];
rz(-pi) q[1];
rz(-0.31089155) q[2];
sx q[2];
rz(-1.0513554) q[2];
sx q[2];
rz(2.3497557) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57945573) q[1];
sx q[1];
rz(-0.39135763) q[1];
sx q[1];
rz(2.1298903) q[1];
x q[2];
rz(-0.22844875) q[3];
sx q[3];
rz(-0.72586021) q[3];
sx q[3];
rz(-1.6231692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.145891) q[2];
sx q[2];
rz(-0.13549165) q[2];
sx q[2];
rz(0.091391407) q[2];
rz(-2.7898096) q[3];
sx q[3];
rz(-0.7472977) q[3];
sx q[3];
rz(-2.6763797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2202989) q[0];
sx q[0];
rz(-1.1722246) q[0];
sx q[0];
rz(1.9660796) q[0];
rz(-2.562404) q[1];
sx q[1];
rz(-0.80668956) q[1];
sx q[1];
rz(-0.18956345) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2611147) q[0];
sx q[0];
rz(-1.3700598) q[0];
sx q[0];
rz(1.6693382) q[0];
rz(-pi) q[1];
rz(2.2722424) q[2];
sx q[2];
rz(-1.5044208) q[2];
sx q[2];
rz(-0.36176031) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.8622562) q[1];
sx q[1];
rz(-1.5784043) q[1];
sx q[1];
rz(-0.81806446) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.060154288) q[3];
sx q[3];
rz(-1.0511479) q[3];
sx q[3];
rz(-1.026767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8976979) q[2];
sx q[2];
rz(-1.8420668) q[2];
sx q[2];
rz(-0.50590903) q[2];
rz(-0.77887744) q[3];
sx q[3];
rz(-2.7454822) q[3];
sx q[3];
rz(-1.2946543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8706354) q[0];
sx q[0];
rz(-2.0896235) q[0];
sx q[0];
rz(1.1554931) q[0];
rz(-0.0011860154) q[1];
sx q[1];
rz(-2.3540034) q[1];
sx q[1];
rz(-2.7364065) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7195994) q[0];
sx q[0];
rz(-1.1047537) q[0];
sx q[0];
rz(0.58902044) q[0];
rz(1.1616906) q[2];
sx q[2];
rz(-1.4295345) q[2];
sx q[2];
rz(-2.2337332) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5444774) q[1];
sx q[1];
rz(-1.6469203) q[1];
sx q[1];
rz(-2.3352469) q[1];
x q[2];
rz(-0.12801068) q[3];
sx q[3];
rz(-1.6388005) q[3];
sx q[3];
rz(-1.8314513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.81795168) q[2];
sx q[2];
rz(-1.1672856) q[2];
sx q[2];
rz(2.6200068) q[2];
rz(-0.096605435) q[3];
sx q[3];
rz(-1.01869) q[3];
sx q[3];
rz(-0.49083138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7535962) q[0];
sx q[0];
rz(-0.8803519) q[0];
sx q[0];
rz(0.11465797) q[0];
rz(1.313734) q[1];
sx q[1];
rz(-1.6856245) q[1];
sx q[1];
rz(0.1350666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5324482) q[0];
sx q[0];
rz(-1.3125064) q[0];
sx q[0];
rz(1.3468955) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7497063) q[2];
sx q[2];
rz(-1.8183437) q[2];
sx q[2];
rz(-0.72617578) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.8986021) q[1];
sx q[1];
rz(-2.1861316) q[1];
sx q[1];
rz(-2.9427285) q[1];
rz(-pi) q[2];
rz(1.4536132) q[3];
sx q[3];
rz(-0.56196594) q[3];
sx q[3];
rz(-2.9148852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8257398) q[2];
sx q[2];
rz(-1.9261381) q[2];
sx q[2];
rz(2.1960171) q[2];
rz(2.8816176) q[3];
sx q[3];
rz(-1.0228478) q[3];
sx q[3];
rz(-2.0135349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0876227) q[0];
sx q[0];
rz(-1.0902417) q[0];
sx q[0];
rz(-1.4029652) q[0];
rz(2.2258017) q[1];
sx q[1];
rz(-2.5251838) q[1];
sx q[1];
rz(2.2260407) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4175891) q[0];
sx q[0];
rz(-1.6678828) q[0];
sx q[0];
rz(-1.7018945) q[0];
rz(-pi) q[1];
rz(-0.19277566) q[2];
sx q[2];
rz(-1.0571684) q[2];
sx q[2];
rz(0.35266587) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8232223) q[1];
sx q[1];
rz(-0.61152667) q[1];
sx q[1];
rz(2.3182858) q[1];
rz(-pi) q[2];
rz(-2.9718982) q[3];
sx q[3];
rz(-2.2732919) q[3];
sx q[3];
rz(0.22512057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.12330684) q[2];
sx q[2];
rz(-2.4985963) q[2];
sx q[2];
rz(2.7389738) q[2];
rz(1.1766524) q[3];
sx q[3];
rz(-2.8549356) q[3];
sx q[3];
rz(2.3736242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9993512) q[0];
sx q[0];
rz(-0.9724697) q[0];
sx q[0];
rz(-1.019626) q[0];
rz(0.6877407) q[1];
sx q[1];
rz(-0.80324695) q[1];
sx q[1];
rz(-1.3245503) q[1];
rz(-2.4551212) q[2];
sx q[2];
rz(-2.1778637) q[2];
sx q[2];
rz(-1.7354497) q[2];
rz(-0.63078337) q[3];
sx q[3];
rz(-1.7787686) q[3];
sx q[3];
rz(-2.6566915) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
