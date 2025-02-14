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
rz(1.3996742) q[0];
sx q[0];
rz(-1.6784486) q[0];
sx q[0];
rz(-1.4611257) q[0];
rz(2.6598016) q[1];
sx q[1];
rz(-2.3269589) q[1];
sx q[1];
rz(2.2021267) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50892539) q[0];
sx q[0];
rz(-1.5804176) q[0];
sx q[0];
rz(1.6144572) q[0];
rz(-2.6109735) q[2];
sx q[2];
rz(-1.8920027) q[2];
sx q[2];
rz(-1.8663592) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.11096) q[1];
sx q[1];
rz(-3.1000376) q[1];
sx q[1];
rz(-1.3135733) q[1];
rz(-pi) q[2];
rz(-1.788673) q[3];
sx q[3];
rz(-1.1718201) q[3];
sx q[3];
rz(-0.50673317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7751094) q[2];
sx q[2];
rz(-1.6747687) q[2];
sx q[2];
rz(-2.5169306) q[2];
rz(0.96056739) q[3];
sx q[3];
rz(-1.4550236) q[3];
sx q[3];
rz(2.2666986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4082773) q[0];
sx q[0];
rz(-0.69077078) q[0];
sx q[0];
rz(-2.9498003) q[0];
rz(0.57418144) q[1];
sx q[1];
rz(-1.5230813) q[1];
sx q[1];
rz(0.40723732) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2977266) q[0];
sx q[0];
rz(-2.9653984) q[0];
sx q[0];
rz(1.6398482) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7474465) q[2];
sx q[2];
rz(-1.6345858) q[2];
sx q[2];
rz(-1.072536) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5246689) q[1];
sx q[1];
rz(-1.3868995) q[1];
sx q[1];
rz(-2.6016629) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.408618) q[3];
sx q[3];
rz(-2.6335178) q[3];
sx q[3];
rz(1.8014419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9521956) q[2];
sx q[2];
rz(-1.1789097) q[2];
sx q[2];
rz(0.65323812) q[2];
rz(0.87853986) q[3];
sx q[3];
rz(-2.0432751) q[3];
sx q[3];
rz(-0.11305913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1427796) q[0];
sx q[0];
rz(-1.7072562) q[0];
sx q[0];
rz(-2.0913731) q[0];
rz(0.6160008) q[1];
sx q[1];
rz(-1.417825) q[1];
sx q[1];
rz(-2.2284257) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2722867) q[0];
sx q[0];
rz(-1.7401631) q[0];
sx q[0];
rz(0.10793751) q[0];
x q[1];
rz(1.9550858) q[2];
sx q[2];
rz(-2.8336513) q[2];
sx q[2];
rz(2.4787087) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.84308594) q[1];
sx q[1];
rz(-1.8131327) q[1];
sx q[1];
rz(-2.0451115) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88190474) q[3];
sx q[3];
rz(-1.8671452) q[3];
sx q[3];
rz(-1.0605845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9203165) q[2];
sx q[2];
rz(-2.2344799) q[2];
sx q[2];
rz(-0.059180666) q[2];
rz(-2.3024998) q[3];
sx q[3];
rz(-1.9234761) q[3];
sx q[3];
rz(2.2742417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3686104) q[0];
sx q[0];
rz(-1.0091857) q[0];
sx q[0];
rz(0.34977812) q[0];
rz(-1.3444258) q[1];
sx q[1];
rz(-0.58660048) q[1];
sx q[1];
rz(-3.0288568) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3679593) q[0];
sx q[0];
rz(-1.2293574) q[0];
sx q[0];
rz(-0.48078791) q[0];
rz(-pi) q[1];
rz(-2.2305528) q[2];
sx q[2];
rz(-1.8000812) q[2];
sx q[2];
rz(-2.2509172) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.73600125) q[1];
sx q[1];
rz(-1.8718157) q[1];
sx q[1];
rz(-0.89717866) q[1];
rz(0.68513667) q[3];
sx q[3];
rz(-2.1129012) q[3];
sx q[3];
rz(0.66451162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.747644) q[2];
sx q[2];
rz(-1.5851721) q[2];
sx q[2];
rz(3.0843688) q[2];
rz(1.8852437) q[3];
sx q[3];
rz(-2.5689954) q[3];
sx q[3];
rz(2.4109139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2644065) q[0];
sx q[0];
rz(-1.974396) q[0];
sx q[0];
rz(-0.71006829) q[0];
rz(-2.8801628) q[1];
sx q[1];
rz(-1.5855887) q[1];
sx q[1];
rz(0.95103055) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3426845) q[0];
sx q[0];
rz(-2.5747188) q[0];
sx q[0];
rz(2.5004412) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0847974) q[2];
sx q[2];
rz(-0.99691454) q[2];
sx q[2];
rz(1.4950372) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5366423) q[1];
sx q[1];
rz(-0.8785696) q[1];
sx q[1];
rz(1.3068024) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.099154648) q[3];
sx q[3];
rz(-0.94010779) q[3];
sx q[3];
rz(0.48336467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6422673) q[2];
sx q[2];
rz(-0.32411164) q[2];
sx q[2];
rz(-2.8322241) q[2];
rz(-1.0950836) q[3];
sx q[3];
rz(-2.0131922) q[3];
sx q[3];
rz(-0.40422082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22769044) q[0];
sx q[0];
rz(-2.9588283) q[0];
sx q[0];
rz(0.90232724) q[0];
rz(-0.062571049) q[1];
sx q[1];
rz(-1.2689271) q[1];
sx q[1];
rz(-1.2455469) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0479054) q[0];
sx q[0];
rz(-1.7249247) q[0];
sx q[0];
rz(-2.4612294) q[0];
rz(-pi) q[1];
rz(1.6643413) q[2];
sx q[2];
rz(-0.54160833) q[2];
sx q[2];
rz(0.53047413) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5426297) q[1];
sx q[1];
rz(-1.0815878) q[1];
sx q[1];
rz(-1.0973003) q[1];
rz(-pi) q[2];
rz(2.0537187) q[3];
sx q[3];
rz(-2.1854162) q[3];
sx q[3];
rz(-2.9642005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41783276) q[2];
sx q[2];
rz(-1.2699026) q[2];
sx q[2];
rz(-1.7977686) q[2];
rz(-1.8580681) q[3];
sx q[3];
rz(-0.44461861) q[3];
sx q[3];
rz(-2.4166687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.035324) q[0];
sx q[0];
rz(-3.1043053) q[0];
sx q[0];
rz(-2.7653747) q[0];
rz(2.6694934) q[1];
sx q[1];
rz(-2.4596228) q[1];
sx q[1];
rz(-2.6825486) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7721586) q[0];
sx q[0];
rz(-0.87352814) q[0];
sx q[0];
rz(-2.5791772) q[0];
rz(-1.3960346) q[2];
sx q[2];
rz(-1.1272001) q[2];
sx q[2];
rz(0.52973739) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9487293) q[1];
sx q[1];
rz(-1.4079844) q[1];
sx q[1];
rz(1.639674) q[1];
rz(-pi) q[2];
rz(0.01851933) q[3];
sx q[3];
rz(-0.59773047) q[3];
sx q[3];
rz(1.2883119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.072319357) q[2];
sx q[2];
rz(-2.4299712) q[2];
sx q[2];
rz(-2.7834328) q[2];
rz(0.88851309) q[3];
sx q[3];
rz(-2.1202959) q[3];
sx q[3];
rz(-2.2108938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096980378) q[0];
sx q[0];
rz(-0.29878578) q[0];
sx q[0];
rz(-0.97691798) q[0];
rz(0.75374976) q[1];
sx q[1];
rz(-1.6981373) q[1];
sx q[1];
rz(-1.6389729) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80874211) q[0];
sx q[0];
rz(-0.70510222) q[0];
sx q[0];
rz(-0.1211959) q[0];
x q[1];
rz(2.4719878) q[2];
sx q[2];
rz(-1.9960224) q[2];
sx q[2];
rz(2.2086672) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5411631) q[1];
sx q[1];
rz(-0.90879295) q[1];
sx q[1];
rz(1.1837436) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94705233) q[3];
sx q[3];
rz(-1.3823518) q[3];
sx q[3];
rz(-0.59977967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.94540191) q[2];
sx q[2];
rz(-1.6795029) q[2];
sx q[2];
rz(0.84575829) q[2];
rz(-0.78019199) q[3];
sx q[3];
rz(-1.3581685) q[3];
sx q[3];
rz(0.63854027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82275259) q[0];
sx q[0];
rz(-1.0718811) q[0];
sx q[0];
rz(0.34291294) q[0];
rz(-2.137939) q[1];
sx q[1];
rz(-1.9099078) q[1];
sx q[1];
rz(0.71663219) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1176004) q[0];
sx q[0];
rz(-2.2870734) q[0];
sx q[0];
rz(2.0407659) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54729531) q[2];
sx q[2];
rz(-2.5792053) q[2];
sx q[2];
rz(2.9862816) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9327859) q[1];
sx q[1];
rz(-0.98742551) q[1];
sx q[1];
rz(1.2482743) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9699205) q[3];
sx q[3];
rz(-2.2033327) q[3];
sx q[3];
rz(0.61033953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0493244) q[2];
sx q[2];
rz(-1.6246139) q[2];
sx q[2];
rz(-2.1410904) q[2];
rz(2.3246824) q[3];
sx q[3];
rz(-1.2714081) q[3];
sx q[3];
rz(0.36417714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3889121) q[0];
sx q[0];
rz(-2.2081544) q[0];
sx q[0];
rz(0.77891427) q[0];
rz(0.79175103) q[1];
sx q[1];
rz(-1.2261032) q[1];
sx q[1];
rz(-2.7955999) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2591922) q[0];
sx q[0];
rz(-1.1244052) q[0];
sx q[0];
rz(1.4784228) q[0];
rz(1.5322288) q[2];
sx q[2];
rz(-0.30353217) q[2];
sx q[2];
rz(0.26798778) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.78948632) q[1];
sx q[1];
rz(-2.5946684) q[1];
sx q[1];
rz(-1.5032306) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84684144) q[3];
sx q[3];
rz(-1.2961642) q[3];
sx q[3];
rz(1.3534897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.148823) q[2];
sx q[2];
rz(-1.6042446) q[2];
sx q[2];
rz(-2.0796622) q[2];
rz(-2.3289833) q[3];
sx q[3];
rz(-1.142623) q[3];
sx q[3];
rz(-0.44212166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.84520311) q[0];
sx q[0];
rz(-0.82798216) q[0];
sx q[0];
rz(1.6430586) q[0];
rz(-2.8614112) q[1];
sx q[1];
rz(-1.9277086) q[1];
sx q[1];
rz(2.3452506) q[1];
rz(1.5308357) q[2];
sx q[2];
rz(-1.162095) q[2];
sx q[2];
rz(3.0065995) q[2];
rz(-0.82874684) q[3];
sx q[3];
rz(-1.256905) q[3];
sx q[3];
rz(-1.4245413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
