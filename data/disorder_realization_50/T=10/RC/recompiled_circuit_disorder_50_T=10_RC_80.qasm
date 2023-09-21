OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.87941909) q[0];
sx q[0];
rz(-1.3949431) q[0];
sx q[0];
rz(-3.1403132) q[0];
rz(1.4446422) q[1];
sx q[1];
rz(-1.1029707) q[1];
sx q[1];
rz(-0.77494088) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2431508) q[0];
sx q[0];
rz(-1.4282707) q[0];
sx q[0];
rz(-1.6210763) q[0];
rz(-pi) q[1];
rz(-2.2968729) q[2];
sx q[2];
rz(-1.2436927) q[2];
sx q[2];
rz(2.6127882) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85841132) q[1];
sx q[1];
rz(-1.5309257) q[1];
sx q[1];
rz(-1.7130501) q[1];
rz(-pi) q[2];
rz(1.7232056) q[3];
sx q[3];
rz(-1.2890352) q[3];
sx q[3];
rz(1.8834653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0960192) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(1.9457031) q[2];
rz(1.9879509) q[3];
sx q[3];
rz(-0.78912815) q[3];
sx q[3];
rz(-1.3886064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7213223) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(1.8700245) q[0];
rz(2.0416073) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(1.3756479) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7398867) q[0];
sx q[0];
rz(-1.4835843) q[0];
sx q[0];
rz(0.4102104) q[0];
rz(1.3324845) q[2];
sx q[2];
rz(-1.5268227) q[2];
sx q[2];
rz(2.6251052) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.980629) q[1];
sx q[1];
rz(-1.7501608) q[1];
sx q[1];
rz(2.7705926) q[1];
rz(2.7820285) q[3];
sx q[3];
rz(-1.0457977) q[3];
sx q[3];
rz(-1.1121225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21330825) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(0.48970547) q[2];
rz(-1.1335763) q[3];
sx q[3];
rz(-0.19530345) q[3];
sx q[3];
rz(-2.074923) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5415444) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(2.9009853) q[0];
rz(2.799017) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(1.2352357) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9440207) q[0];
sx q[0];
rz(-2.4320142) q[0];
sx q[0];
rz(-2.5463085) q[0];
rz(1.3182993) q[2];
sx q[2];
rz(-1.2090948) q[2];
sx q[2];
rz(-0.30644882) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.55360824) q[1];
sx q[1];
rz(-1.8318892) q[1];
sx q[1];
rz(1.8334465) q[1];
x q[2];
rz(-2.5123185) q[3];
sx q[3];
rz(-0.73262239) q[3];
sx q[3];
rz(-2.7936558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.76413313) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(-1.4366478) q[2];
rz(1.4012339) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(-2.5333372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075994611) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(-0.80379379) q[0];
rz(0.94961387) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(-3.0217357) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0474931) q[0];
sx q[0];
rz(-2.1764628) q[0];
sx q[0];
rz(3.0981307) q[0];
rz(-pi) q[1];
rz(-1.4255964) q[2];
sx q[2];
rz(-1.2135266) q[2];
sx q[2];
rz(1.4841339) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2377492) q[1];
sx q[1];
rz(-1.3386054) q[1];
sx q[1];
rz(3.1234427) q[1];
x q[2];
rz(2.1060137) q[3];
sx q[3];
rz(-1.03973) q[3];
sx q[3];
rz(-2.5553184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5084761) q[2];
sx q[2];
rz(-1.2458331) q[2];
sx q[2];
rz(-0.072337739) q[2];
rz(-2.7667601) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(1.1981296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(-0.48686349) q[0];
rz(-2.411719) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(1.9015076) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70532521) q[0];
sx q[0];
rz(-1.8099394) q[0];
sx q[0];
rz(1.5572085) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24057062) q[2];
sx q[2];
rz(-1.8261989) q[2];
sx q[2];
rz(1.0742998) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8023194) q[1];
sx q[1];
rz(-0.60802751) q[1];
sx q[1];
rz(3.0250711) q[1];
rz(-pi) q[2];
rz(2.1754856) q[3];
sx q[3];
rz(-1.9786414) q[3];
sx q[3];
rz(-0.84433489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1029677) q[2];
sx q[2];
rz(-0.90927783) q[2];
sx q[2];
rz(2.4385578) q[2];
rz(-1.4098343) q[3];
sx q[3];
rz(-1.792428) q[3];
sx q[3];
rz(-0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1773961) q[0];
sx q[0];
rz(-1.2921768) q[0];
sx q[0];
rz(0.99037209) q[0];
rz(-3.0888427) q[1];
sx q[1];
rz(-0.92664781) q[1];
sx q[1];
rz(1.4809158) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4303362) q[0];
sx q[0];
rz(-0.38861409) q[0];
sx q[0];
rz(-1.3769763) q[0];
rz(-2.0159857) q[2];
sx q[2];
rz(-1.5700969) q[2];
sx q[2];
rz(2.6631151) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33038352) q[1];
sx q[1];
rz(-0.57302176) q[1];
sx q[1];
rz(-1.1793544) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.048599676) q[3];
sx q[3];
rz(-0.89266333) q[3];
sx q[3];
rz(0.62852678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1331553) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(2.5793502) q[2];
rz(1.0605313) q[3];
sx q[3];
rz(-2.3980467) q[3];
sx q[3];
rz(2.8806768) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9437207) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(-1.6725756) q[0];
rz(2.127227) q[1];
sx q[1];
rz(-2.1140153) q[1];
sx q[1];
rz(1.8168824) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0990471) q[0];
sx q[0];
rz(-1.5481661) q[0];
sx q[0];
rz(0.12734539) q[0];
x q[1];
rz(2.6323694) q[2];
sx q[2];
rz(-1.0668796) q[2];
sx q[2];
rz(0.06171209) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.087254698) q[1];
sx q[1];
rz(-1.195765) q[1];
sx q[1];
rz(-2.7021728) q[1];
rz(-pi) q[2];
rz(-2.2277101) q[3];
sx q[3];
rz(-1.2607288) q[3];
sx q[3];
rz(3.0403746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5614732) q[2];
sx q[2];
rz(-1.3886398) q[2];
sx q[2];
rz(0.44357792) q[2];
rz(-2.1929072) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(-2.366812) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7396486) q[0];
sx q[0];
rz(-1.0111324) q[0];
sx q[0];
rz(2.6810714) q[0];
rz(3.0415688) q[1];
sx q[1];
rz(-0.99383751) q[1];
sx q[1];
rz(1.8519648) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8136918) q[0];
sx q[0];
rz(-2.0223589) q[0];
sx q[0];
rz(-0.62664647) q[0];
rz(2.4256267) q[2];
sx q[2];
rz(-2.1858474) q[2];
sx q[2];
rz(-1.8383319) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.928927) q[1];
sx q[1];
rz(-1.7576828) q[1];
sx q[1];
rz(3.0287663) q[1];
x q[2];
rz(1.7361705) q[3];
sx q[3];
rz(-0.91455063) q[3];
sx q[3];
rz(-1.0429494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.148968) q[2];
sx q[2];
rz(-0.31105369) q[2];
sx q[2];
rz(1.0827433) q[2];
rz(-0.096171245) q[3];
sx q[3];
rz(-1.440719) q[3];
sx q[3];
rz(1.910803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37725317) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(0.33426958) q[0];
rz(1.9175247) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(2.8589378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44336549) q[0];
sx q[0];
rz(-2.3400314) q[0];
sx q[0];
rz(0.65936868) q[0];
x q[1];
rz(2.8106868) q[2];
sx q[2];
rz(-1.3661642) q[2];
sx q[2];
rz(-0.76588878) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4947195) q[1];
sx q[1];
rz(-2.6962526) q[1];
sx q[1];
rz(-2.7237707) q[1];
rz(-pi) q[2];
rz(-2.2563124) q[3];
sx q[3];
rz(-1.4797398) q[3];
sx q[3];
rz(2.1212999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9294372) q[2];
sx q[2];
rz(-1.1634469) q[2];
sx q[2];
rz(2.4003417) q[2];
rz(-0.50179982) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(-0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.042645) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(-0.50869554) q[0];
rz(0.11518654) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(0.68181109) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1334575) q[0];
sx q[0];
rz(-1.7383766) q[0];
sx q[0];
rz(-2.125678) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82294686) q[2];
sx q[2];
rz(-1.032885) q[2];
sx q[2];
rz(0.049494628) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0360003) q[1];
sx q[1];
rz(-0.75165527) q[1];
sx q[1];
rz(2.4113301) q[1];
x q[2];
rz(-0.52313157) q[3];
sx q[3];
rz(-2.0683859) q[3];
sx q[3];
rz(0.79062068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.29601413) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(-0.34109035) q[2];
rz(2.0579445) q[3];
sx q[3];
rz(-2.3458979) q[3];
sx q[3];
rz(-0.17102374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1682128) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(-2.534261) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(0.66424673) q[2];
sx q[2];
rz(-1.1462117) q[2];
sx q[2];
rz(0.29427634) q[2];
rz(-1.3689465) q[3];
sx q[3];
rz(-1.9484083) q[3];
sx q[3];
rz(1.2231135) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];