OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.26602715) q[0];
sx q[0];
rz(-0.53524435) q[0];
sx q[0];
rz(0.75403655) q[0];
rz(-5.4929805) q[1];
sx q[1];
rz(5.0561855) q[1];
sx q[1];
rz(8.2639134) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55863419) q[0];
sx q[0];
rz(-1.6003803) q[0];
sx q[0];
rz(2.6628859) q[0];
rz(-pi) q[1];
rz(0.16205807) q[2];
sx q[2];
rz(-1.098212) q[2];
sx q[2];
rz(-0.87640793) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4143715) q[1];
sx q[1];
rz(-1.9804269) q[1];
sx q[1];
rz(-1.0135256) q[1];
x q[2];
rz(0.17920223) q[3];
sx q[3];
rz(-1.7342907) q[3];
sx q[3];
rz(-1.4737827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51241088) q[2];
sx q[2];
rz(-1.8449731) q[2];
sx q[2];
rz(-0.66649246) q[2];
rz(2.6317821) q[3];
sx q[3];
rz(-1.1504268) q[3];
sx q[3];
rz(1.2734909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.78012413) q[0];
sx q[0];
rz(-1.402401) q[0];
sx q[0];
rz(-1.1126888) q[0];
rz(2.9878222) q[1];
sx q[1];
rz(-2.230481) q[1];
sx q[1];
rz(-1.8033093) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3801549) q[0];
sx q[0];
rz(-2.2614711) q[0];
sx q[0];
rz(-1.3209016) q[0];
x q[1];
rz(2.0427225) q[2];
sx q[2];
rz(-2.3510691) q[2];
sx q[2];
rz(0.84004842) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3160623) q[1];
sx q[1];
rz(-1.9426553) q[1];
sx q[1];
rz(-1.4447681) q[1];
rz(0.23985858) q[3];
sx q[3];
rz(-2.5407102) q[3];
sx q[3];
rz(-0.89752737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0587557) q[2];
sx q[2];
rz(-2.357491) q[2];
sx q[2];
rz(-1.178297) q[2];
rz(-0.96238771) q[3];
sx q[3];
rz(-2.048384) q[3];
sx q[3];
rz(0.66550955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1798379) q[0];
sx q[0];
rz(-2.6340155) q[0];
sx q[0];
rz(-1.5270365) q[0];
rz(-0.64287341) q[1];
sx q[1];
rz(-1.0765272) q[1];
sx q[1];
rz(-0.33338526) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6282435) q[0];
sx q[0];
rz(-1.3643364) q[0];
sx q[0];
rz(2.2702361) q[0];
x q[1];
rz(-1.9579499) q[2];
sx q[2];
rz(-1.6946812) q[2];
sx q[2];
rz(-1.3882335) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.15409878) q[1];
sx q[1];
rz(-1.0990267) q[1];
sx q[1];
rz(3.0497453) q[1];
rz(-pi) q[2];
rz(-1.3473347) q[3];
sx q[3];
rz(-0.94933214) q[3];
sx q[3];
rz(-0.37924757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.12144111) q[2];
sx q[2];
rz(-1.8276428) q[2];
sx q[2];
rz(2.3392759) q[2];
rz(2.9004167) q[3];
sx q[3];
rz(-0.69883385) q[3];
sx q[3];
rz(-0.10087092) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0543095) q[0];
sx q[0];
rz(-1.5232975) q[0];
sx q[0];
rz(2.5774082) q[0];
rz(2.5634649) q[1];
sx q[1];
rz(-1.4923948) q[1];
sx q[1];
rz(-0.50813466) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3343737) q[0];
sx q[0];
rz(-2.5426572) q[0];
sx q[0];
rz(0.81143023) q[0];
x q[1];
rz(2.5035985) q[2];
sx q[2];
rz(-1.5882512) q[2];
sx q[2];
rz(-2.2538315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8453464) q[1];
sx q[1];
rz(-2.2987662) q[1];
sx q[1];
rz(-1.216757) q[1];
rz(-pi) q[2];
rz(0.99677892) q[3];
sx q[3];
rz(-0.32696163) q[3];
sx q[3];
rz(0.60861482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6886787) q[2];
sx q[2];
rz(-2.8082509) q[2];
sx q[2];
rz(2.508146) q[2];
rz(-0.59988919) q[3];
sx q[3];
rz(-1.9918631) q[3];
sx q[3];
rz(-1.5002804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3863581) q[0];
sx q[0];
rz(-2.8383377) q[0];
sx q[0];
rz(-2.9454943) q[0];
rz(-1.261699) q[1];
sx q[1];
rz(-0.82048565) q[1];
sx q[1];
rz(-2.0702147) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5539861) q[0];
sx q[0];
rz(-1.1452132) q[0];
sx q[0];
rz(-2.2527184) q[0];
rz(1.1159665) q[2];
sx q[2];
rz(-2.358846) q[2];
sx q[2];
rz(2.6775529) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.57429806) q[1];
sx q[1];
rz(-1.8540566) q[1];
sx q[1];
rz(1.1351372) q[1];
rz(-pi) q[2];
rz(-2.3134872) q[3];
sx q[3];
rz(-2.5525186) q[3];
sx q[3];
rz(2.7554054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.034996899) q[2];
sx q[2];
rz(-1.804616) q[2];
sx q[2];
rz(0.98199797) q[2];
rz(-2.9563831) q[3];
sx q[3];
rz(-2.2976112) q[3];
sx q[3];
rz(-1.265032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5380602) q[0];
sx q[0];
rz(-2.2221727) q[0];
sx q[0];
rz(-2.61125) q[0];
rz(1.2999339) q[1];
sx q[1];
rz(-1.329774) q[1];
sx q[1];
rz(2.9249654) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3743065) q[0];
sx q[0];
rz(-1.8348798) q[0];
sx q[0];
rz(0.81861511) q[0];
x q[1];
rz(1.8703307) q[2];
sx q[2];
rz(-0.96200633) q[2];
sx q[2];
rz(-1.8744206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4649012) q[1];
sx q[1];
rz(-1.562582) q[1];
sx q[1];
rz(-1.2354922) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2413007) q[3];
sx q[3];
rz(-1.4196463) q[3];
sx q[3];
rz(2.8063262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4914322) q[2];
sx q[2];
rz(-1.6314793) q[2];
sx q[2];
rz(1.023863) q[2];
rz(-2.2653545) q[3];
sx q[3];
rz(-0.70228464) q[3];
sx q[3];
rz(-1.8700301) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4890471) q[0];
sx q[0];
rz(-1.9788195) q[0];
sx q[0];
rz(2.6126557) q[0];
rz(-1.5286998) q[1];
sx q[1];
rz(-1.9493608) q[1];
sx q[1];
rz(1.0891917) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3750317) q[0];
sx q[0];
rz(-1.8683109) q[0];
sx q[0];
rz(3.1356698) q[0];
x q[1];
rz(-2.2839374) q[2];
sx q[2];
rz(-2.0115888) q[2];
sx q[2];
rz(-2.8105274) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.293922) q[1];
sx q[1];
rz(-1.5772181) q[1];
sx q[1];
rz(0.033172219) q[1];
x q[2];
rz(0.090737061) q[3];
sx q[3];
rz(-1.4416579) q[3];
sx q[3];
rz(1.4108301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.12525325) q[2];
sx q[2];
rz(-0.80417997) q[2];
sx q[2];
rz(1.9160697) q[2];
rz(1.4922173) q[3];
sx q[3];
rz(-1.7047313) q[3];
sx q[3];
rz(0.15587458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4847223) q[0];
sx q[0];
rz(-3.0796034) q[0];
sx q[0];
rz(-2.2739676) q[0];
rz(3.0743657) q[1];
sx q[1];
rz(-1.0311238) q[1];
sx q[1];
rz(-0.19518383) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0474437) q[0];
sx q[0];
rz(-0.35956811) q[0];
sx q[0];
rz(-0.6028428) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5848265) q[2];
sx q[2];
rz(-1.1696891) q[2];
sx q[2];
rz(-0.48497981) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.47123779) q[1];
sx q[1];
rz(-1.1413304) q[1];
sx q[1];
rz(-3.0958789) q[1];
rz(-pi) q[2];
rz(-2.1956452) q[3];
sx q[3];
rz(-2.5759856) q[3];
sx q[3];
rz(0.92029508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.247867) q[2];
sx q[2];
rz(-0.21394955) q[2];
sx q[2];
rz(-0.90551886) q[2];
rz(1.9838105) q[3];
sx q[3];
rz(-1.1685305) q[3];
sx q[3];
rz(-2.982443) q[3];
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
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4147707) q[0];
sx q[0];
rz(-2.045571) q[0];
sx q[0];
rz(-1.0409521) q[0];
rz(-3.0629311) q[1];
sx q[1];
rz(-0.18053308) q[1];
sx q[1];
rz(0.35531607) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8126292) q[0];
sx q[0];
rz(-1.2041429) q[0];
sx q[0];
rz(-3.0366412) q[0];
x q[1];
rz(1.1460733) q[2];
sx q[2];
rz(-0.5364843) q[2];
sx q[2];
rz(-1.2338961) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3231343) q[1];
sx q[1];
rz(-1.6272021) q[1];
sx q[1];
rz(1.8947381) q[1];
x q[2];
rz(2.5411685) q[3];
sx q[3];
rz(-1.7250337) q[3];
sx q[3];
rz(-2.141181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4650402) q[2];
sx q[2];
rz(-0.6074473) q[2];
sx q[2];
rz(1.4036277) q[2];
rz(-3.1353531) q[3];
sx q[3];
rz(-1.798636) q[3];
sx q[3];
rz(-1.5374373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86826098) q[0];
sx q[0];
rz(-1.1678168) q[0];
sx q[0];
rz(1.2982752) q[0];
rz(2.6955993) q[1];
sx q[1];
rz(-2.0263717) q[1];
sx q[1];
rz(2.840852) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8200127) q[0];
sx q[0];
rz(-1.0386779) q[0];
sx q[0];
rz(-0.26474712) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5955574) q[2];
sx q[2];
rz(-1.9249501) q[2];
sx q[2];
rz(-1.491577) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4423351) q[1];
sx q[1];
rz(-1.8619969) q[1];
sx q[1];
rz(-1.553781) q[1];
rz(-pi) q[2];
rz(-0.91046393) q[3];
sx q[3];
rz(-1.2372036) q[3];
sx q[3];
rz(-0.098284483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5131502) q[2];
sx q[2];
rz(-2.0437888) q[2];
sx q[2];
rz(-2.002031) q[2];
rz(1.3509753) q[3];
sx q[3];
rz(-0.59777483) q[3];
sx q[3];
rz(-1.1736419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.448485) q[0];
sx q[0];
rz(-1.2379452) q[0];
sx q[0];
rz(-2.2647279) q[0];
rz(1.7383472) q[1];
sx q[1];
rz(-1.2259903) q[1];
sx q[1];
rz(-1.2798053) q[1];
rz(3.1115361) q[2];
sx q[2];
rz(-3.0855784) q[2];
sx q[2];
rz(-2.2655178) q[2];
rz(-0.26294796) q[3];
sx q[3];
rz(-2.1750952) q[3];
sx q[3];
rz(1.7261214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
