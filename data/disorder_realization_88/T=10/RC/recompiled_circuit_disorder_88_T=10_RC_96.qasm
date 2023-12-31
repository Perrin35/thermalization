OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.24703439) q[0];
sx q[0];
rz(-1.8859099) q[0];
sx q[0];
rz(0.32796252) q[0];
rz(2.9070931) q[1];
sx q[1];
rz(-0.20107888) q[1];
sx q[1];
rz(-3.0501563) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5782195) q[0];
sx q[0];
rz(-1.9679929) q[0];
sx q[0];
rz(0.63011516) q[0];
rz(-pi) q[1];
rz(0.25912164) q[2];
sx q[2];
rz(-1.7012351) q[2];
sx q[2];
rz(2.5082617) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.92802231) q[1];
sx q[1];
rz(-1.808128) q[1];
sx q[1];
rz(3.0568394) q[1];
rz(-pi) q[2];
rz(-1.2535291) q[3];
sx q[3];
rz(-2.5228365) q[3];
sx q[3];
rz(-1.3878824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4771007) q[2];
sx q[2];
rz(-2.1680919) q[2];
sx q[2];
rz(2.0155902) q[2];
rz(2.8664355) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(-0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1317516) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(2.4480208) q[0];
rz(2.0454848) q[1];
sx q[1];
rz(-0.98384905) q[1];
sx q[1];
rz(-0.19031659) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2454525) q[0];
sx q[0];
rz(-1.3724694) q[0];
sx q[0];
rz(-1.1597) q[0];
rz(-pi) q[1];
rz(0.53493494) q[2];
sx q[2];
rz(-1.2282279) q[2];
sx q[2];
rz(1.554622) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7029876) q[1];
sx q[1];
rz(-2.4191796) q[1];
sx q[1];
rz(-1.5370675) q[1];
rz(-pi) q[2];
rz(2.042949) q[3];
sx q[3];
rz(-0.90001366) q[3];
sx q[3];
rz(-2.9428279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
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
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.398657) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(-0.8202585) q[0];
rz(2.8495158) q[1];
sx q[1];
rz(-2.0675817) q[1];
sx q[1];
rz(1.2480199) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67278636) q[0];
sx q[0];
rz(-2.5350223) q[0];
sx q[0];
rz(-2.9787105) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6349995) q[2];
sx q[2];
rz(-1.3639796) q[2];
sx q[2];
rz(1.9931672) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6229912) q[1];
sx q[1];
rz(-2.2366183) q[1];
sx q[1];
rz(-1.6559421) q[1];
x q[2];
rz(-2.2218024) q[3];
sx q[3];
rz(-1.3295104) q[3];
sx q[3];
rz(-1.1834708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.32039207) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(-1.9281663) q[2];
rz(-2.976867) q[3];
sx q[3];
rz(-0.89403331) q[3];
sx q[3];
rz(-0.081610672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.7320025) q[0];
sx q[0];
rz(-1.2816757) q[0];
sx q[0];
rz(0.18606342) q[0];
rz(-0.20446725) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(-1.2971372) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064917795) q[0];
sx q[0];
rz(-3.0041822) q[0];
sx q[0];
rz(-1.8924367) q[0];
rz(-1.859971) q[2];
sx q[2];
rz(-1.7316069) q[2];
sx q[2];
rz(-0.5493872) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.101673) q[1];
sx q[1];
rz(-0.97517255) q[1];
sx q[1];
rz(-0.69570978) q[1];
x q[2];
rz(-1.6468871) q[3];
sx q[3];
rz(-0.18008672) q[3];
sx q[3];
rz(-2.0538581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90594784) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(0.91919351) q[2];
rz(-2.8202608) q[3];
sx q[3];
rz(-2.0733757) q[3];
sx q[3];
rz(1.2478158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17386757) q[0];
sx q[0];
rz(-2.6350832) q[0];
sx q[0];
rz(2.2633973) q[0];
rz(-1.8163266) q[1];
sx q[1];
rz(-1.7293431) q[1];
sx q[1];
rz(1.7153046) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8844937) q[0];
sx q[0];
rz(-1.6528659) q[0];
sx q[0];
rz(-0.9135855) q[0];
rz(-pi) q[1];
rz(-2.2141453) q[2];
sx q[2];
rz(-2.7577835) q[2];
sx q[2];
rz(-0.62188934) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.95477415) q[1];
sx q[1];
rz(-0.84328077) q[1];
sx q[1];
rz(0.14528841) q[1];
rz(-1.250508) q[3];
sx q[3];
rz(-0.8257782) q[3];
sx q[3];
rz(0.68009963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8695801) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(-3.0997979) q[2];
rz(-0.061491866) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(0.4367691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3549266) q[0];
sx q[0];
rz(-3.0556398) q[0];
sx q[0];
rz(1.0304931) q[0];
rz(-2.4018535) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(0.57156634) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58493462) q[0];
sx q[0];
rz(-2.1437763) q[0];
sx q[0];
rz(-0.95227382) q[0];
rz(-3.0145054) q[2];
sx q[2];
rz(-2.2327144) q[2];
sx q[2];
rz(1.454929) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.39396573) q[1];
sx q[1];
rz(-0.52226258) q[1];
sx q[1];
rz(-2.0678492) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4404293) q[3];
sx q[3];
rz(-0.94651604) q[3];
sx q[3];
rz(-3.0091156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.17343865) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(-2.5773933) q[2];
rz(3.0155904) q[3];
sx q[3];
rz(-1.6832451) q[3];
sx q[3];
rz(-2.8469767) q[3];
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
rz(1.0903704) q[0];
sx q[0];
rz(-2.9416961) q[0];
sx q[0];
rz(0.71227658) q[0];
rz(-0.5258711) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(-2.4760822) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3411361) q[0];
sx q[0];
rz(-1.8130214) q[0];
sx q[0];
rz(-3.1405297) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7475142) q[2];
sx q[2];
rz(-0.70509796) q[2];
sx q[2];
rz(1.2202028) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4588786) q[1];
sx q[1];
rz(-2.6280118) q[1];
sx q[1];
rz(-1.1125803) q[1];
rz(-2.9969278) q[3];
sx q[3];
rz(-1.8105227) q[3];
sx q[3];
rz(-0.02558115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9843288) q[2];
sx q[2];
rz(-0.66005808) q[2];
sx q[2];
rz(-1.4228014) q[2];
rz(0.11519365) q[3];
sx q[3];
rz(-0.51518232) q[3];
sx q[3];
rz(0.19259024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049906235) q[0];
sx q[0];
rz(-2.005907) q[0];
sx q[0];
rz(-0.19083047) q[0];
rz(-0.62675369) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(0.33871067) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3074293) q[0];
sx q[0];
rz(-1.2273664) q[0];
sx q[0];
rz(2.9817581) q[0];
x q[1];
rz(0.45192265) q[2];
sx q[2];
rz(-1.3616614) q[2];
sx q[2];
rz(-2.506633) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.415886) q[1];
sx q[1];
rz(-0.48735122) q[1];
sx q[1];
rz(-1.8094256) q[1];
rz(2.268928) q[3];
sx q[3];
rz(-0.41655311) q[3];
sx q[3];
rz(2.0779028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.64615858) q[2];
sx q[2];
rz(-2.50714) q[2];
sx q[2];
rz(2.4411566) q[2];
rz(-0.8979848) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(0.17351304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53196466) q[0];
sx q[0];
rz(-2.6599929) q[0];
sx q[0];
rz(0.65834808) q[0];
rz(0.61093962) q[1];
sx q[1];
rz(-1.2555723) q[1];
sx q[1];
rz(3.0019965) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0848207) q[0];
sx q[0];
rz(-0.064282566) q[0];
sx q[0];
rz(-1.2540741) q[0];
x q[1];
rz(-1.6330958) q[2];
sx q[2];
rz(-2.3844516) q[2];
sx q[2];
rz(1.3296668) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86233625) q[1];
sx q[1];
rz(-1.6376054) q[1];
sx q[1];
rz(0.52211296) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2152882) q[3];
sx q[3];
rz(-1.5288121) q[3];
sx q[3];
rz(0.55461649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6167986) q[2];
sx q[2];
rz(-2.1308265) q[2];
sx q[2];
rz(0.33111462) q[2];
rz(-0.75774276) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(-0.087879114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963294) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(2.9560126) q[0];
rz(-2.045385) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(-1.4846444) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0508704) q[0];
sx q[0];
rz(-1.8403247) q[0];
sx q[0];
rz(-1.1140633) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8469641) q[2];
sx q[2];
rz(-1.3535) q[2];
sx q[2];
rz(1.5751788) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1440925) q[1];
sx q[1];
rz(-1.1966685) q[1];
sx q[1];
rz(3.0348026) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1366523) q[3];
sx q[3];
rz(-0.86588174) q[3];
sx q[3];
rz(-2.4019394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2075656) q[2];
sx q[2];
rz(-0.88946122) q[2];
sx q[2];
rz(-2.5893842) q[2];
rz(-2.3637555) q[3];
sx q[3];
rz(-2.2838897) q[3];
sx q[3];
rz(-2.623693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(1.4245695) q[2];
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
