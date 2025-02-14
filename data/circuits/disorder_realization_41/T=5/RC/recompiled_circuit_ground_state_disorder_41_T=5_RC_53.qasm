OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5858894) q[0];
sx q[0];
rz(-1.1893505) q[0];
sx q[0];
rz(-1.1638292) q[0];
rz(2.3948506) q[1];
sx q[1];
rz(-2.876694) q[1];
sx q[1];
rz(2.9704111) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25397444) q[0];
sx q[0];
rz(-1.3518392) q[0];
sx q[0];
rz(2.9848214) q[0];
rz(-pi) q[1];
rz(0.88340448) q[2];
sx q[2];
rz(-0.92057487) q[2];
sx q[2];
rz(-1.5387077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.65003465) q[1];
sx q[1];
rz(-1.2753873) q[1];
sx q[1];
rz(-0.00062382767) q[1];
x q[2];
rz(1.5653651) q[3];
sx q[3];
rz(-1.587638) q[3];
sx q[3];
rz(1.4726437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.050194) q[2];
sx q[2];
rz(-1.5768496) q[2];
sx q[2];
rz(2.0719299) q[2];
rz(-1.4403249) q[3];
sx q[3];
rz(-1.0245208) q[3];
sx q[3];
rz(1.9694156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1759724) q[0];
sx q[0];
rz(-1.4485285) q[0];
sx q[0];
rz(-2.5536221) q[0];
rz(-1.8486456) q[1];
sx q[1];
rz(-0.99423948) q[1];
sx q[1];
rz(2.9867244) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17093824) q[0];
sx q[0];
rz(-2.0679255) q[0];
sx q[0];
rz(1.6927822) q[0];
rz(1.7172408) q[2];
sx q[2];
rz(-1.3322209) q[2];
sx q[2];
rz(1.6171136) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5649287) q[1];
sx q[1];
rz(-1.1047078) q[1];
sx q[1];
rz(-2.7647777) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9294176) q[3];
sx q[3];
rz(-1.4373693) q[3];
sx q[3];
rz(-1.9047036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.64314848) q[2];
sx q[2];
rz(-0.73068205) q[2];
sx q[2];
rz(-0.11239642) q[2];
rz(-2.3183838) q[3];
sx q[3];
rz(-1.9197437) q[3];
sx q[3];
rz(2.6196151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4560029) q[0];
sx q[0];
rz(-2.0528448) q[0];
sx q[0];
rz(-2.1734557) q[0];
rz(2.0227506) q[1];
sx q[1];
rz(-1.6066931) q[1];
sx q[1];
rz(1.3333837) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4844167) q[0];
sx q[0];
rz(-1.4098995) q[0];
sx q[0];
rz(-0.96477525) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5381375) q[2];
sx q[2];
rz(-1.1759182) q[2];
sx q[2];
rz(-1.6445356) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3231877) q[1];
sx q[1];
rz(-1.0099704) q[1];
sx q[1];
rz(-3.0238999) q[1];
rz(-pi) q[2];
rz(1.2127688) q[3];
sx q[3];
rz(-1.6505401) q[3];
sx q[3];
rz(0.84313508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.94990388) q[2];
sx q[2];
rz(-2.1046941) q[2];
sx q[2];
rz(0.26491586) q[2];
rz(-0.31629899) q[3];
sx q[3];
rz(-2.8079872) q[3];
sx q[3];
rz(1.3914289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.815627) q[0];
sx q[0];
rz(-0.2206603) q[0];
sx q[0];
rz(1.6478446) q[0];
rz(-1.7296467) q[1];
sx q[1];
rz(-1.0271881) q[1];
sx q[1];
rz(0.6573917) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7058168) q[0];
sx q[0];
rz(-0.12744337) q[0];
sx q[0];
rz(-0.61690046) q[0];
rz(0.98720647) q[2];
sx q[2];
rz(-2.2650717) q[2];
sx q[2];
rz(-2.2615446) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.58739793) q[1];
sx q[1];
rz(-0.19853354) q[1];
sx q[1];
rz(1.3819329) q[1];
x q[2];
rz(-0.811565) q[3];
sx q[3];
rz(-0.92716427) q[3];
sx q[3];
rz(-2.0763458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.21939453) q[2];
sx q[2];
rz(-1.9191091) q[2];
sx q[2];
rz(-2.738319) q[2];
rz(-1.9379617) q[3];
sx q[3];
rz(-1.6324685) q[3];
sx q[3];
rz(-0.50886124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68172425) q[0];
sx q[0];
rz(-0.42400703) q[0];
sx q[0];
rz(-0.62393171) q[0];
rz(0.722018) q[1];
sx q[1];
rz(-1.7941509) q[1];
sx q[1];
rz(-3.0375979) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.390966) q[0];
sx q[0];
rz(-1.4292882) q[0];
sx q[0];
rz(0.69049375) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7816167) q[2];
sx q[2];
rz(-1.2449387) q[2];
sx q[2];
rz(2.6687372) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5370868) q[1];
sx q[1];
rz(-2.1963122) q[1];
sx q[1];
rz(0.64670678) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4173043) q[3];
sx q[3];
rz(-2.2632627) q[3];
sx q[3];
rz(2.9812733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35088745) q[2];
sx q[2];
rz(-0.71983379) q[2];
sx q[2];
rz(0.78901115) q[2];
rz(-1.8770494) q[3];
sx q[3];
rz(-1.0435373) q[3];
sx q[3];
rz(0.98715442) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3037381) q[0];
sx q[0];
rz(-0.721295) q[0];
sx q[0];
rz(-0.20027941) q[0];
rz(1.4664949) q[1];
sx q[1];
rz(-1.8148345) q[1];
sx q[1];
rz(-1.6237578) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12833327) q[0];
sx q[0];
rz(-2.4966514) q[0];
sx q[0];
rz(-0.63887424) q[0];
x q[1];
rz(-2.0744223) q[2];
sx q[2];
rz(-1.6852411) q[2];
sx q[2];
rz(-1.9527854) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7519546) q[1];
sx q[1];
rz(-1.5742969) q[1];
sx q[1];
rz(0.59094999) q[1];
x q[2];
rz(0.79849859) q[3];
sx q[3];
rz(-1.7424462) q[3];
sx q[3];
rz(0.6397748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.62310654) q[2];
sx q[2];
rz(-0.97272626) q[2];
sx q[2];
rz(2.8002807) q[2];
rz(0.85706472) q[3];
sx q[3];
rz(-1.2271481) q[3];
sx q[3];
rz(-2.6763693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.051006) q[0];
sx q[0];
rz(-2.107928) q[0];
sx q[0];
rz(1.2475913) q[0];
rz(-0.11820758) q[1];
sx q[1];
rz(-1.2756313) q[1];
sx q[1];
rz(-0.035621312) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8310332) q[0];
sx q[0];
rz(-0.11971029) q[0];
sx q[0];
rz(-1.381078) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9526677) q[2];
sx q[2];
rz(-1.1729203) q[2];
sx q[2];
rz(0.16824761) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.25917398) q[1];
sx q[1];
rz(-0.66773623) q[1];
sx q[1];
rz(-2.7386433) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3603975) q[3];
sx q[3];
rz(-1.6577814) q[3];
sx q[3];
rz(1.1802105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9901765) q[2];
sx q[2];
rz(-1.6665062) q[2];
sx q[2];
rz(2.0659633) q[2];
rz(-2.6208124) q[3];
sx q[3];
rz(-2.5213089) q[3];
sx q[3];
rz(-1.552593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2369279) q[0];
sx q[0];
rz(-1.0228782) q[0];
sx q[0];
rz(0.9032332) q[0];
rz(1.6385993) q[1];
sx q[1];
rz(-1.4475977) q[1];
sx q[1];
rz(-1.2618056) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5905268) q[0];
sx q[0];
rz(-1.3414978) q[0];
sx q[0];
rz(0.14044827) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4440358) q[2];
sx q[2];
rz(-0.64854014) q[2];
sx q[2];
rz(1.7096138) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3751531) q[1];
sx q[1];
rz(-0.10442153) q[1];
sx q[1];
rz(-0.5071284) q[1];
rz(-pi) q[2];
rz(0.60537548) q[3];
sx q[3];
rz(-1.9937111) q[3];
sx q[3];
rz(1.133906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.5791851) q[2];
sx q[2];
rz(-2.4895442) q[2];
sx q[2];
rz(-0.11422608) q[2];
rz(1.7646029) q[3];
sx q[3];
rz(-1.1403133) q[3];
sx q[3];
rz(2.637114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1219516) q[0];
sx q[0];
rz(-0.99948245) q[0];
sx q[0];
rz(-1.1720538) q[0];
rz(-1.3837586) q[1];
sx q[1];
rz(-1.9969321) q[1];
sx q[1];
rz(-3.0269472) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8269236) q[0];
sx q[0];
rz(-2.0631587) q[0];
sx q[0];
rz(-0.63757245) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4769606) q[2];
sx q[2];
rz(-2.0619446) q[2];
sx q[2];
rz(1.067184) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.35797) q[1];
sx q[1];
rz(-1.8183892) q[1];
sx q[1];
rz(-0.072958306) q[1];
rz(-pi) q[2];
rz(0.10520868) q[3];
sx q[3];
rz(-0.79920879) q[3];
sx q[3];
rz(-1.5463643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5961479) q[2];
sx q[2];
rz(-2.9344276) q[2];
sx q[2];
rz(-2.2607415) q[2];
rz(-2.7018069) q[3];
sx q[3];
rz(-1.9763016) q[3];
sx q[3];
rz(-2.0161276) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83774829) q[0];
sx q[0];
rz(-1.4864018) q[0];
sx q[0];
rz(-3.0434171) q[0];
rz(2.5671666) q[1];
sx q[1];
rz(-1.6822633) q[1];
sx q[1];
rz(0.59304768) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3705854) q[0];
sx q[0];
rz(-2.2789776) q[0];
sx q[0];
rz(-2.7561491) q[0];
x q[1];
rz(0.93043987) q[2];
sx q[2];
rz(-2.499492) q[2];
sx q[2];
rz(3.0885027) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9657085) q[1];
sx q[1];
rz(-1.5595601) q[1];
sx q[1];
rz(1.5574993) q[1];
x q[2];
rz(0.63238588) q[3];
sx q[3];
rz(-1.8416071) q[3];
sx q[3];
rz(2.2512844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1526996) q[2];
sx q[2];
rz(-1.1836735) q[2];
sx q[2];
rz(2.2031671) q[2];
rz(1.3931795) q[3];
sx q[3];
rz(-1.3104855) q[3];
sx q[3];
rz(-0.2383298) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68168454) q[0];
sx q[0];
rz(-0.60612283) q[0];
sx q[0];
rz(-1.7932307) q[0];
rz(-1.0472736) q[1];
sx q[1];
rz(-2.2127163) q[1];
sx q[1];
rz(-0.97102078) q[1];
rz(0.045108724) q[2];
sx q[2];
rz(-1.6354297) q[2];
sx q[2];
rz(2.9387504) q[2];
rz(0.18235124) q[3];
sx q[3];
rz(-1.3538881) q[3];
sx q[3];
rz(2.8543579) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
