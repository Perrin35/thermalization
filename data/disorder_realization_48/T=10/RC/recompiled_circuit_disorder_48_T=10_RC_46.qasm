OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7735908) q[0];
sx q[0];
rz(3.9324023) q[0];
sx q[0];
rz(12.232236) q[0];
rz(2.6842527) q[1];
sx q[1];
rz(-2.1973124) q[1];
sx q[1];
rz(-1.2184719) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1760362) q[0];
sx q[0];
rz(-0.82057014) q[0];
sx q[0];
rz(1.6198938) q[0];
rz(-pi) q[1];
rz(2.676744) q[2];
sx q[2];
rz(-1.6072304) q[2];
sx q[2];
rz(0.51793098) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.54108809) q[1];
sx q[1];
rz(-2.0588015) q[1];
sx q[1];
rz(1.0469251) q[1];
rz(2.9524607) q[3];
sx q[3];
rz(-1.3526219) q[3];
sx q[3];
rz(1.8060341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5228287) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(-0.5775601) q[2];
rz(1.1497568) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(0.66453385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1537271) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(-0.38744774) q[0];
rz(2.2024343) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(-1.4025677) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1419066) q[0];
sx q[0];
rz(-1.5748595) q[0];
sx q[0];
rz(3.1121029) q[0];
x q[1];
rz(-0.55949975) q[2];
sx q[2];
rz(-1.3845978) q[2];
sx q[2];
rz(-2.6452243) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.95784159) q[1];
sx q[1];
rz(-1.1693923) q[1];
sx q[1];
rz(0.5387696) q[1];
rz(-1.6658695) q[3];
sx q[3];
rz(-2.2007696) q[3];
sx q[3];
rz(2.7271987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7188321) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(0.31769162) q[2];
rz(-2.9348532) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(2.3247705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7725672) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(1.7279708) q[0];
rz(2.6638022) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(0.40107045) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7558407) q[0];
sx q[0];
rz(-2.7999561) q[0];
sx q[0];
rz(1.2582448) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7602073) q[2];
sx q[2];
rz(-0.93511287) q[2];
sx q[2];
rz(0.96178255) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.89951) q[1];
sx q[1];
rz(-2.3009355) q[1];
sx q[1];
rz(-1.5622557) q[1];
x q[2];
rz(0.19604711) q[3];
sx q[3];
rz(-1.184433) q[3];
sx q[3];
rz(-2.8783609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.59427375) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(-2.5857914) q[2];
rz(-2.1650971) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(0.78021375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34898409) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(1.4439616) q[0];
rz(-1.6216888) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(2.8881883) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4356416) q[0];
sx q[0];
rz(-0.55103978) q[0];
sx q[0];
rz(-1.3735229) q[0];
x q[1];
rz(-1.0850111) q[2];
sx q[2];
rz(-1.5554785) q[2];
sx q[2];
rz(-0.77335301) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.464444) q[1];
sx q[1];
rz(-0.943999) q[1];
sx q[1];
rz(-1.4518751) q[1];
rz(-0.6494135) q[3];
sx q[3];
rz(-1.6522437) q[3];
sx q[3];
rz(2.7922975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.110934) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(0.78732642) q[2];
rz(-0.91286719) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(1.0095989) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71516365) q[0];
sx q[0];
rz(-2.5029095) q[0];
sx q[0];
rz(-0.062967904) q[0];
rz(-0.12403034) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(-2.6834992) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3375181) q[0];
sx q[0];
rz(-0.44263698) q[0];
sx q[0];
rz(-0.0054365693) q[0];
x q[1];
rz(2.6890254) q[2];
sx q[2];
rz(-1.0107702) q[2];
sx q[2];
rz(-2.3134311) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.749436) q[1];
sx q[1];
rz(-1.4538308) q[1];
sx q[1];
rz(2.3163296) q[1];
x q[2];
rz(3.0466988) q[3];
sx q[3];
rz(-2.0695544) q[3];
sx q[3];
rz(-0.14548485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9690341) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(-2.9120973) q[2];
rz(0.0028006639) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(1.3389448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(2.2221185) q[0];
rz(-3.0793076) q[1];
sx q[1];
rz(-2.1376164) q[1];
sx q[1];
rz(-1.8744291) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29612449) q[0];
sx q[0];
rz(-2.2757029) q[0];
sx q[0];
rz(-2.092093) q[0];
rz(-pi) q[1];
rz(1.7153347) q[2];
sx q[2];
rz(-0.83206165) q[2];
sx q[2];
rz(-1.2046255) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3082723) q[1];
sx q[1];
rz(-0.38494021) q[1];
sx q[1];
rz(-1.0233378) q[1];
x q[2];
rz(-2.0503067) q[3];
sx q[3];
rz(-1.7877842) q[3];
sx q[3];
rz(-3.019141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1798114) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(0.58376694) q[2];
rz(2.4328655) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(-3.1183929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0631183) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(2.069058) q[0];
rz(0.5468927) q[1];
sx q[1];
rz(-1.8976338) q[1];
sx q[1];
rz(2.0297208) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81387732) q[0];
sx q[0];
rz(-2.0534671) q[0];
sx q[0];
rz(-3.0359603) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9260336) q[2];
sx q[2];
rz(-1.0676427) q[2];
sx q[2];
rz(-0.1728729) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8060311) q[1];
sx q[1];
rz(-2.8415488) q[1];
sx q[1];
rz(-2.8801444) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92915793) q[3];
sx q[3];
rz(-0.60141364) q[3];
sx q[3];
rz(-1.2777559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.83773461) q[2];
sx q[2];
rz(-1.6656817) q[2];
sx q[2];
rz(1.973935) q[2];
rz(-1.5363103) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(-1.7355828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7325571) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(-2.4107966) q[0];
rz(-0.90019512) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(-0.75497595) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98159957) q[0];
sx q[0];
rz(-0.85708517) q[0];
sx q[0];
rz(-0.1937565) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72889363) q[2];
sx q[2];
rz(-1.1155827) q[2];
sx q[2];
rz(0.58065562) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6983812) q[1];
sx q[1];
rz(-1.4581919) q[1];
sx q[1];
rz(0.57971445) q[1];
rz(-0.98324361) q[3];
sx q[3];
rz(-1.1508815) q[3];
sx q[3];
rz(2.0411232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76453152) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(-2.5047452) q[2];
rz(2.8751255) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72717845) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(2.0027347) q[0];
rz(-2.3873734) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(3.1220904) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0855904) q[0];
sx q[0];
rz(-1.7770355) q[0];
sx q[0];
rz(-3.0598559) q[0];
rz(0.15140622) q[2];
sx q[2];
rz(-1.3726241) q[2];
sx q[2];
rz(-2.1422447) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2724185) q[1];
sx q[1];
rz(-2.4915016) q[1];
sx q[1];
rz(1.8723349) q[1];
rz(-1.8832302) q[3];
sx q[3];
rz(-1.7282681) q[3];
sx q[3];
rz(-1.0851932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1372244) q[2];
sx q[2];
rz(-1.7251816) q[2];
sx q[2];
rz(2.2926889) q[2];
rz(0.38765872) q[3];
sx q[3];
rz(-1.1281697) q[3];
sx q[3];
rz(1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7983109) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(2.6570901) q[0];
rz(-1.7548521) q[1];
sx q[1];
rz(-1.4258899) q[1];
sx q[1];
rz(-1.1482931) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1185547) q[0];
sx q[0];
rz(-1.430129) q[0];
sx q[0];
rz(-1.5492357) q[0];
rz(1.3482773) q[2];
sx q[2];
rz(-0.95138022) q[2];
sx q[2];
rz(-2.7811108) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3155568) q[1];
sx q[1];
rz(-2.450374) q[1];
sx q[1];
rz(-1.8408937) q[1];
x q[2];
rz(-0.14264543) q[3];
sx q[3];
rz(-0.35257617) q[3];
sx q[3];
rz(0.25776097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2293573) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(-0.36515507) q[2];
rz(0.12864104) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0317595) q[0];
sx q[0];
rz(-0.84072996) q[0];
sx q[0];
rz(1.6054556) q[0];
rz(0.96314349) q[1];
sx q[1];
rz(-1.8704725) q[1];
sx q[1];
rz(2.0830547) q[1];
rz(2.4776496) q[2];
sx q[2];
rz(-1.2524111) q[2];
sx q[2];
rz(-1.8128957) q[2];
rz(0.036988463) q[3];
sx q[3];
rz(-1.0109517) q[3];
sx q[3];
rz(3.0299822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];