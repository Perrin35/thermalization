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
rz(-1.3396076) q[0];
sx q[0];
rz(2.7924502) q[0];
sx q[0];
rz(10.452441) q[0];
rz(-0.034962058) q[1];
sx q[1];
rz(-1.0171913) q[1];
sx q[1];
rz(-1.4453759) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.611269) q[0];
sx q[0];
rz(-1.9113386) q[0];
sx q[0];
rz(-0.61113071) q[0];
x q[1];
rz(-1.6048172) q[2];
sx q[2];
rz(-1.7482702) q[2];
sx q[2];
rz(0.79171514) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45856341) q[1];
sx q[1];
rz(-1.7818319) q[1];
sx q[1];
rz(0.33331916) q[1];
rz(-pi) q[2];
rz(1.9361922) q[3];
sx q[3];
rz(-2.2789848) q[3];
sx q[3];
rz(2.303249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8029636) q[2];
sx q[2];
rz(-2.6821319) q[2];
sx q[2];
rz(2.6098693) q[2];
rz(1.7470597) q[3];
sx q[3];
rz(-1.2732384) q[3];
sx q[3];
rz(1.9664221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26674536) q[0];
sx q[0];
rz(-1.5044455) q[0];
sx q[0];
rz(2.3496085) q[0];
rz(-0.92164552) q[1];
sx q[1];
rz(-1.4896723) q[1];
sx q[1];
rz(-2.0335061) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9757248) q[0];
sx q[0];
rz(-2.168403) q[0];
sx q[0];
rz(3.0063666) q[0];
rz(-pi) q[1];
rz(-2.1943548) q[2];
sx q[2];
rz(-1.531336) q[2];
sx q[2];
rz(0.90496162) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.82906065) q[1];
sx q[1];
rz(-2.2680813) q[1];
sx q[1];
rz(-0.50364772) q[1];
rz(-pi) q[2];
rz(-3.0854598) q[3];
sx q[3];
rz(-1.2548459) q[3];
sx q[3];
rz(-2.454284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7532928) q[2];
sx q[2];
rz(-2.1123977) q[2];
sx q[2];
rz(-1.8886867) q[2];
rz(0.62475723) q[3];
sx q[3];
rz(-1.2478991) q[3];
sx q[3];
rz(0.46318444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4598292) q[0];
sx q[0];
rz(-2.9637931) q[0];
sx q[0];
rz(-2.8681927) q[0];
rz(0.071648486) q[1];
sx q[1];
rz(-1.1337846) q[1];
sx q[1];
rz(-1.515548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0296299) q[0];
sx q[0];
rz(-0.66034895) q[0];
sx q[0];
rz(0.65159722) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.927761) q[2];
sx q[2];
rz(-2.0139815) q[2];
sx q[2];
rz(0.24493327) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.30574894) q[1];
sx q[1];
rz(-1.4121404) q[1];
sx q[1];
rz(0.58802998) q[1];
rz(-pi) q[2];
rz(0.25392308) q[3];
sx q[3];
rz(-0.99231595) q[3];
sx q[3];
rz(-1.8691269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7800954) q[2];
sx q[2];
rz(-2.4269035) q[2];
sx q[2];
rz(1.3698461) q[2];
rz(-1.0684377) q[3];
sx q[3];
rz(-1.7504102) q[3];
sx q[3];
rz(2.1865602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84900981) q[0];
sx q[0];
rz(-0.42790258) q[0];
sx q[0];
rz(-0.12246116) q[0];
rz(-0.017008688) q[1];
sx q[1];
rz(-1.6288792) q[1];
sx q[1];
rz(-2.2564127) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1830782) q[0];
sx q[0];
rz(-1.4349271) q[0];
sx q[0];
rz(2.0625173) q[0];
rz(-pi) q[1];
rz(-0.34220747) q[2];
sx q[2];
rz(-0.86538431) q[2];
sx q[2];
rz(2.9660385) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.95751132) q[1];
sx q[1];
rz(-2.5169417) q[1];
sx q[1];
rz(2.2923325) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70192317) q[3];
sx q[3];
rz(-1.7198017) q[3];
sx q[3];
rz(1.5535136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.433832) q[2];
sx q[2];
rz(-1.8261352) q[2];
sx q[2];
rz(0.07746499) q[2];
rz(-2.7205983) q[3];
sx q[3];
rz(-1.1869895) q[3];
sx q[3];
rz(0.25119701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0635327) q[0];
sx q[0];
rz(-3.1227626) q[0];
sx q[0];
rz(-2.4596762) q[0];
rz(-0.49184999) q[1];
sx q[1];
rz(-1.0268772) q[1];
sx q[1];
rz(-1.9815725) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9043961) q[0];
sx q[0];
rz(-0.86441308) q[0];
sx q[0];
rz(-0.24817962) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5142646) q[2];
sx q[2];
rz(-1.0011359) q[2];
sx q[2];
rz(-2.4084072) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8012538) q[1];
sx q[1];
rz(-1.7427708) q[1];
sx q[1];
rz(-2.1992963) q[1];
x q[2];
rz(0.97461583) q[3];
sx q[3];
rz(-0.49562956) q[3];
sx q[3];
rz(-2.0295124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68449768) q[2];
sx q[2];
rz(-2.2726077) q[2];
sx q[2];
rz(2.1333466) q[2];
rz(1.6393939) q[3];
sx q[3];
rz(-2.0366171) q[3];
sx q[3];
rz(2.8777299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098792583) q[0];
sx q[0];
rz(-1.5473939) q[0];
sx q[0];
rz(-2.7368326) q[0];
rz(0.81819355) q[1];
sx q[1];
rz(-1.300756) q[1];
sx q[1];
rz(-2.4249446) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7966864) q[0];
sx q[0];
rz(-0.27936882) q[0];
sx q[0];
rz(-3.1059103) q[0];
x q[1];
rz(2.3776182) q[2];
sx q[2];
rz(-1.9696376) q[2];
sx q[2];
rz(-3.0032436) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67988746) q[1];
sx q[1];
rz(-2.6033834) q[1];
sx q[1];
rz(-2.6160282) q[1];
rz(-pi) q[2];
rz(0.57315196) q[3];
sx q[3];
rz(-2.5226197) q[3];
sx q[3];
rz(3.1123146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3396259) q[2];
sx q[2];
rz(-0.64029396) q[2];
sx q[2];
rz(0.65529811) q[2];
rz(-0.35704923) q[3];
sx q[3];
rz(-1.7506295) q[3];
sx q[3];
rz(-2.6220139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94012564) q[0];
sx q[0];
rz(-2.8253912) q[0];
sx q[0];
rz(-1.0200208) q[0];
rz(-2.5992498) q[1];
sx q[1];
rz(-2.0261363) q[1];
sx q[1];
rz(1.7592336) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4750299) q[0];
sx q[0];
rz(-2.2336166) q[0];
sx q[0];
rz(-0.7866089) q[0];
rz(2.5711201) q[2];
sx q[2];
rz(-2.8920724) q[2];
sx q[2];
rz(0.0970627) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0497363) q[1];
sx q[1];
rz(-1.1106792) q[1];
sx q[1];
rz(-0.087398296) q[1];
rz(-pi) q[2];
rz(3.1145688) q[3];
sx q[3];
rz(-0.083649717) q[3];
sx q[3];
rz(-0.25318709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.072784) q[2];
sx q[2];
rz(-1.5280318) q[2];
sx q[2];
rz(0.35045785) q[2];
rz(2.6751878) q[3];
sx q[3];
rz(-2.2240708) q[3];
sx q[3];
rz(-2.1980227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35895178) q[0];
sx q[0];
rz(-0.39059165) q[0];
sx q[0];
rz(-2.1851831) q[0];
rz(-1.4495173) q[1];
sx q[1];
rz(-2.3608975) q[1];
sx q[1];
rz(-2.4704959) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8050552) q[0];
sx q[0];
rz(-1.4292246) q[0];
sx q[0];
rz(-0.086351589) q[0];
rz(-pi) q[1];
x q[1];
rz(0.568957) q[2];
sx q[2];
rz(-1.3929741) q[2];
sx q[2];
rz(0.59019719) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.073428) q[1];
sx q[1];
rz(-0.35919093) q[1];
sx q[1];
rz(-1.753162) q[1];
x q[2];
rz(-2.8550034) q[3];
sx q[3];
rz(-2.4169888) q[3];
sx q[3];
rz(2.2514908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.87390071) q[2];
sx q[2];
rz(-2.716422) q[2];
sx q[2];
rz(2.0153866) q[2];
rz(2.8113484) q[3];
sx q[3];
rz(-1.4010022) q[3];
sx q[3];
rz(0.11708524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0432334) q[0];
sx q[0];
rz(-0.57906228) q[0];
sx q[0];
rz(0.057057127) q[0];
rz(-3.0077899) q[1];
sx q[1];
rz(-2.0963142) q[1];
sx q[1];
rz(0.71279508) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063449115) q[0];
sx q[0];
rz(-2.8196555) q[0];
sx q[0];
rz(0.19970317) q[0];
x q[1];
rz(-1.9781238) q[2];
sx q[2];
rz(-1.1962657) q[2];
sx q[2];
rz(-0.56331149) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1755029) q[1];
sx q[1];
rz(-2.0398629) q[1];
sx q[1];
rz(1.1203946) q[1];
x q[2];
rz(2.2650293) q[3];
sx q[3];
rz(-1.6329771) q[3];
sx q[3];
rz(-2.2060195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5130676) q[2];
sx q[2];
rz(-1.8551989) q[2];
sx q[2];
rz(-3.0600424) q[2];
rz(-1.2201803) q[3];
sx q[3];
rz(-0.27118513) q[3];
sx q[3];
rz(1.0903821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14209014) q[0];
sx q[0];
rz(-1.5220078) q[0];
sx q[0];
rz(-1.581544) q[0];
rz(-2.0060495) q[1];
sx q[1];
rz(-2.0027497) q[1];
sx q[1];
rz(1.1467689) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0479779) q[0];
sx q[0];
rz(-0.85913697) q[0];
sx q[0];
rz(-2.467903) q[0];
x q[1];
rz(-0.79305834) q[2];
sx q[2];
rz(-0.93843776) q[2];
sx q[2];
rz(-1.5989433) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4301871) q[1];
sx q[1];
rz(-1.780695) q[1];
sx q[1];
rz(-1.485157) q[1];
rz(-pi) q[2];
rz(-0.087298079) q[3];
sx q[3];
rz(-1.0074769) q[3];
sx q[3];
rz(0.47730744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.59209383) q[2];
sx q[2];
rz(-1.8213976) q[2];
sx q[2];
rz(2.7756694) q[2];
rz(-0.38771114) q[3];
sx q[3];
rz(-2.0230899) q[3];
sx q[3];
rz(-1.6754735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6112919) q[0];
sx q[0];
rz(-0.74321754) q[0];
sx q[0];
rz(2.8331953) q[0];
rz(0.74956924) q[1];
sx q[1];
rz(-2.0105965) q[1];
sx q[1];
rz(-1.2731193) q[1];
rz(2.7992934) q[2];
sx q[2];
rz(-2.2020349) q[2];
sx q[2];
rz(2.0943691) q[2];
rz(2.5091631) q[3];
sx q[3];
rz(-1.2865744) q[3];
sx q[3];
rz(-2.1435973) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
