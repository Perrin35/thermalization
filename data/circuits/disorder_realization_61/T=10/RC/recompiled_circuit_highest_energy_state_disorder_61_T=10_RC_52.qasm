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
rz(-1.6725809) q[0];
sx q[0];
rz(-2.2218158) q[0];
sx q[0];
rz(2.4866009) q[0];
rz(-1.4870149) q[1];
sx q[1];
rz(-1.4227285) q[1];
sx q[1];
rz(1.3230327) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5736377) q[0];
sx q[0];
rz(-2.1186189) q[0];
sx q[0];
rz(1.9325038) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62686129) q[2];
sx q[2];
rz(-0.38598362) q[2];
sx q[2];
rz(2.4311993) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3016794) q[1];
sx q[1];
rz(-0.29407802) q[1];
sx q[1];
rz(-1.2204351) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5526974) q[3];
sx q[3];
rz(-0.82020226) q[3];
sx q[3];
rz(0.39862788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.46372867) q[2];
sx q[2];
rz(-2.4552796) q[2];
sx q[2];
rz(-2.4832671) q[2];
rz(-2.5074734) q[3];
sx q[3];
rz(-1.442797) q[3];
sx q[3];
rz(1.3707976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089791678) q[0];
sx q[0];
rz(-2.8103516) q[0];
sx q[0];
rz(3.064503) q[0];
rz(2.5568621) q[1];
sx q[1];
rz(-1.3405776) q[1];
sx q[1];
rz(-1.9416521) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9532419) q[0];
sx q[0];
rz(-2.6380499) q[0];
sx q[0];
rz(1.153323) q[0];
x q[1];
rz(-1.7351772) q[2];
sx q[2];
rz(-1.5998249) q[2];
sx q[2];
rz(2.8433702) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.89948326) q[1];
sx q[1];
rz(-0.98190183) q[1];
sx q[1];
rz(-0.45581006) q[1];
rz(-pi) q[2];
rz(2.582483) q[3];
sx q[3];
rz(-1.5190795) q[3];
sx q[3];
rz(2.8962108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8889019) q[2];
sx q[2];
rz(-0.98068792) q[2];
sx q[2];
rz(2.4661031) q[2];
rz(0.0096970079) q[3];
sx q[3];
rz(-2.8986425) q[3];
sx q[3];
rz(-0.01827904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49552396) q[0];
sx q[0];
rz(-1.4718453) q[0];
sx q[0];
rz(-0.75463265) q[0];
rz(-3.1309639) q[1];
sx q[1];
rz(-1.3898712) q[1];
sx q[1];
rz(-1.0108112) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1026238) q[0];
sx q[0];
rz(-1.083923) q[0];
sx q[0];
rz(0.34619934) q[0];
x q[1];
rz(-1.6057683) q[2];
sx q[2];
rz(-0.681923) q[2];
sx q[2];
rz(-2.4833895) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.52151187) q[1];
sx q[1];
rz(-0.82640582) q[1];
sx q[1];
rz(-1.6537731) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1819886) q[3];
sx q[3];
rz(-0.90590796) q[3];
sx q[3];
rz(1.5184107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.90091554) q[2];
sx q[2];
rz(-0.42833504) q[2];
sx q[2];
rz(-2.6348616) q[2];
rz(-1.2725376) q[3];
sx q[3];
rz(-1.3576018) q[3];
sx q[3];
rz(0.27819628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5946567) q[0];
sx q[0];
rz(-1.2242067) q[0];
sx q[0];
rz(1.1154255) q[0];
rz(-1.8755272) q[1];
sx q[1];
rz(-1.5645809) q[1];
sx q[1];
rz(0.87475264) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0717569) q[0];
sx q[0];
rz(-1.6594145) q[0];
sx q[0];
rz(-1.8231234) q[0];
x q[1];
rz(0.018847887) q[2];
sx q[2];
rz(-0.91721469) q[2];
sx q[2];
rz(0.0073405618) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0656888) q[1];
sx q[1];
rz(-1.836628) q[1];
sx q[1];
rz(-0.42215438) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.277209) q[3];
sx q[3];
rz(-1.4339851) q[3];
sx q[3];
rz(-2.9446326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0227585) q[2];
sx q[2];
rz(-0.47670445) q[2];
sx q[2];
rz(-1.9258707) q[2];
rz(1.7440965) q[3];
sx q[3];
rz(-1.5236866) q[3];
sx q[3];
rz(1.7106748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27595156) q[0];
sx q[0];
rz(-0.10198274) q[0];
sx q[0];
rz(1.4707461) q[0];
rz(1.3784846) q[1];
sx q[1];
rz(-2.3536317) q[1];
sx q[1];
rz(1.4102304) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8277692) q[0];
sx q[0];
rz(-0.54049379) q[0];
sx q[0];
rz(-1.9825516) q[0];
x q[1];
rz(2.8431434) q[2];
sx q[2];
rz(-0.69576925) q[2];
sx q[2];
rz(-1.6844105) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8157364) q[1];
sx q[1];
rz(-1.5724025) q[1];
sx q[1];
rz(-2.8291507) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35378176) q[3];
sx q[3];
rz(-1.1615586) q[3];
sx q[3];
rz(1.080846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0146279) q[2];
sx q[2];
rz(-2.2829599) q[2];
sx q[2];
rz(1.3237759) q[2];
rz(-1.7823559) q[3];
sx q[3];
rz(-1.2991354) q[3];
sx q[3];
rz(0.38715473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1112261) q[0];
sx q[0];
rz(-1.4596326) q[0];
sx q[0];
rz(0.02221814) q[0];
rz(-2.4413595) q[1];
sx q[1];
rz(-1.8902238) q[1];
sx q[1];
rz(-0.95692316) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2882521) q[0];
sx q[0];
rz(-1.7159675) q[0];
sx q[0];
rz(-2.4397544) q[0];
rz(-pi) q[1];
rz(-1.3889203) q[2];
sx q[2];
rz(-1.3564912) q[2];
sx q[2];
rz(-0.3265115) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0007504) q[1];
sx q[1];
rz(-0.4943119) q[1];
sx q[1];
rz(0.96338455) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1311573) q[3];
sx q[3];
rz(-1.2613858) q[3];
sx q[3];
rz(0.65195665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9243098) q[2];
sx q[2];
rz(-2.0192912) q[2];
sx q[2];
rz(1.699532) q[2];
rz(1.1066655) q[3];
sx q[3];
rz(-2.536074) q[3];
sx q[3];
rz(1.9895915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.7149413) q[0];
sx q[0];
rz(-1.277667) q[0];
sx q[0];
rz(0.46698025) q[0];
rz(-1.3106208) q[1];
sx q[1];
rz(-0.95342031) q[1];
sx q[1];
rz(2.7511645) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6735839) q[0];
sx q[0];
rz(-1.4303741) q[0];
sx q[0];
rz(-1.6037081) q[0];
x q[1];
rz(-1.8653581) q[2];
sx q[2];
rz(-1.4678363) q[2];
sx q[2];
rz(-2.0077133) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5891287) q[1];
sx q[1];
rz(-0.93163449) q[1];
sx q[1];
rz(0.57560779) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3620699) q[3];
sx q[3];
rz(-0.62297076) q[3];
sx q[3];
rz(1.2869715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9499669) q[2];
sx q[2];
rz(-1.9248631) q[2];
sx q[2];
rz(1.0170271) q[2];
rz(-0.93938604) q[3];
sx q[3];
rz(-2.2685969) q[3];
sx q[3];
rz(0.92923195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.7153213) q[0];
sx q[0];
rz(-0.66522288) q[0];
sx q[0];
rz(-1.2128879) q[0];
rz(-0.60316482) q[1];
sx q[1];
rz(-2.0096571) q[1];
sx q[1];
rz(2.718198) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0807665) q[0];
sx q[0];
rz(-0.87196022) q[0];
sx q[0];
rz(0.10159512) q[0];
rz(-pi) q[1];
rz(-3.064761) q[2];
sx q[2];
rz(-1.4203826) q[2];
sx q[2];
rz(0.39932775) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4881542) q[1];
sx q[1];
rz(-0.27380015) q[1];
sx q[1];
rz(2.4523712) q[1];
rz(-0.16984197) q[3];
sx q[3];
rz(-1.5780996) q[3];
sx q[3];
rz(2.9057692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49070552) q[2];
sx q[2];
rz(-0.31279534) q[2];
sx q[2];
rz(2.1631964) q[2];
rz(-2.4466416) q[3];
sx q[3];
rz(-1.2000822) q[3];
sx q[3];
rz(-1.2945226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0438743) q[0];
sx q[0];
rz(-0.23625034) q[0];
sx q[0];
rz(0.043721113) q[0];
rz(1.1737191) q[1];
sx q[1];
rz(-0.81467384) q[1];
sx q[1];
rz(-2.3497605) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4872015) q[0];
sx q[0];
rz(-2.374361) q[0];
sx q[0];
rz(2.5228398) q[0];
x q[1];
rz(1.2709684) q[2];
sx q[2];
rz(-1.3770896) q[2];
sx q[2];
rz(1.9676339) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2481195) q[1];
sx q[1];
rz(-1.9933102) q[1];
sx q[1];
rz(-0.36568195) q[1];
rz(-2.6641162) q[3];
sx q[3];
rz(-1.9100185) q[3];
sx q[3];
rz(2.4216258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.23058471) q[2];
sx q[2];
rz(-2.8318996) q[2];
sx q[2];
rz(-1.9045551) q[2];
rz(-2.6593995) q[3];
sx q[3];
rz(-0.92415205) q[3];
sx q[3];
rz(2.6528416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.024254) q[0];
sx q[0];
rz(-0.92702213) q[0];
sx q[0];
rz(-1.3379958) q[0];
rz(0.85211873) q[1];
sx q[1];
rz(-1.9133277) q[1];
sx q[1];
rz(-1.9911912) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.710708) q[0];
sx q[0];
rz(-1.0477433) q[0];
sx q[0];
rz(-0.56351785) q[0];
rz(1.8878292) q[2];
sx q[2];
rz(-1.6765127) q[2];
sx q[2];
rz(-2.6144165) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2472966) q[1];
sx q[1];
rz(-1.9622231) q[1];
sx q[1];
rz(2.9904234) q[1];
rz(-pi) q[2];
rz(2.1086333) q[3];
sx q[3];
rz(-2.4696484) q[3];
sx q[3];
rz(0.6368466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.368025) q[2];
sx q[2];
rz(-1.0057534) q[2];
sx q[2];
rz(-2.8311938) q[2];
rz(-2.5144905) q[3];
sx q[3];
rz(-0.98395625) q[3];
sx q[3];
rz(-1.1712801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2720168) q[0];
sx q[0];
rz(-1.5740812) q[0];
sx q[0];
rz(1.5480702) q[0];
rz(1.9238453) q[1];
sx q[1];
rz(-1.4532614) q[1];
sx q[1];
rz(-0.90167602) q[1];
rz(-2.0704253) q[2];
sx q[2];
rz(-1.7555321) q[2];
sx q[2];
rz(2.1499014) q[2];
rz(1.6321833) q[3];
sx q[3];
rz(-1.8922378) q[3];
sx q[3];
rz(0.031382244) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
