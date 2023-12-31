OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(-3.0598109) q[0];
sx q[0];
rz(2.6401289) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(0.3224386) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48243385) q[0];
sx q[0];
rz(-1.6560439) q[0];
sx q[0];
rz(-1.2486588) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2519977) q[2];
sx q[2];
rz(-2.0686364) q[2];
sx q[2];
rz(1.3537458) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0281801) q[1];
sx q[1];
rz(-2.1537188) q[1];
sx q[1];
rz(0.41863538) q[1];
x q[2];
rz(-0.93482165) q[3];
sx q[3];
rz(-1.9536195) q[3];
sx q[3];
rz(-2.8179907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.50513187) q[2];
sx q[2];
rz(-0.59288609) q[2];
sx q[2];
rz(-0.55603975) q[2];
rz(0.83267823) q[3];
sx q[3];
rz(-1.6502389) q[3];
sx q[3];
rz(0.94579831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6933724) q[0];
sx q[0];
rz(-1.6813261) q[0];
sx q[0];
rz(-0.15727501) q[0];
rz(-2.8804624) q[1];
sx q[1];
rz(-1.3477247) q[1];
sx q[1];
rz(0.10903407) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5058274) q[0];
sx q[0];
rz(-2.8553537) q[0];
sx q[0];
rz(2.2155227) q[0];
x q[1];
rz(0.26308665) q[2];
sx q[2];
rz(-1.7619942) q[2];
sx q[2];
rz(0.69743246) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.82169689) q[1];
sx q[1];
rz(-2.3488455) q[1];
sx q[1];
rz(-0.68383033) q[1];
rz(-pi) q[2];
rz(2.845329) q[3];
sx q[3];
rz(-2.6014572) q[3];
sx q[3];
rz(1.6903282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9033501) q[2];
sx q[2];
rz(-1.976333) q[2];
sx q[2];
rz(-1.8781352) q[2];
rz(0.3271099) q[3];
sx q[3];
rz(-1.5771022) q[3];
sx q[3];
rz(1.2143149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064421244) q[0];
sx q[0];
rz(-3.0922958) q[0];
sx q[0];
rz(1.7984614) q[0];
rz(-2.893977) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(-2.6599191) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0945064) q[0];
sx q[0];
rz(-1.1990093) q[0];
sx q[0];
rz(-2.0949754) q[0];
rz(-pi) q[1];
rz(0.9573612) q[2];
sx q[2];
rz(-1.8463677) q[2];
sx q[2];
rz(-2.7797109) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.958365) q[1];
sx q[1];
rz(-2.1124358) q[1];
sx q[1];
rz(0.72802131) q[1];
x q[2];
rz(-2.5211469) q[3];
sx q[3];
rz(-2.1944322) q[3];
sx q[3];
rz(-0.40943957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3383011) q[2];
sx q[2];
rz(-0.81739601) q[2];
sx q[2];
rz(-0.49989191) q[2];
rz(0.56097427) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(-1.4612173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-2.9445779) q[0];
sx q[0];
rz(-0.16600969) q[0];
sx q[0];
rz(2.5894077) q[0];
rz(1.588297) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(1.8968556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9024076) q[0];
sx q[0];
rz(-0.33948487) q[0];
sx q[0];
rz(-3.1367338) q[0];
rz(-1.3696026) q[2];
sx q[2];
rz(-1.6987213) q[2];
sx q[2];
rz(0.97845562) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1340449) q[1];
sx q[1];
rz(-1.782685) q[1];
sx q[1];
rz(2.120554) q[1];
rz(-pi) q[2];
rz(2.7987715) q[3];
sx q[3];
rz(-0.59844136) q[3];
sx q[3];
rz(-3.0330021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.84919471) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(1.1506895) q[2];
rz(-1.6644647) q[3];
sx q[3];
rz(-1.632558) q[3];
sx q[3];
rz(-0.48294827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078159049) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(-0.081469014) q[0];
rz(0.062462656) q[1];
sx q[1];
rz(-1.1413347) q[1];
sx q[1];
rz(-1.5030456) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1254079) q[0];
sx q[0];
rz(-2.5446919) q[0];
sx q[0];
rz(-1.3422658) q[0];
rz(0.083085255) q[2];
sx q[2];
rz(-1.6787046) q[2];
sx q[2];
rz(-1.8748869) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0283708) q[1];
sx q[1];
rz(-1.8824667) q[1];
sx q[1];
rz(0.10722864) q[1];
rz(-pi) q[2];
rz(1.8893858) q[3];
sx q[3];
rz(-1.0831523) q[3];
sx q[3];
rz(2.3209751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6333255) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(-1.903669) q[2];
rz(-2.0189019) q[3];
sx q[3];
rz(-2.4653547) q[3];
sx q[3];
rz(-0.52156633) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6102819) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(2.8748728) q[0];
rz(-2.5807014) q[1];
sx q[1];
rz(-1.2979049) q[1];
sx q[1];
rz(2.3430603) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079433867) q[0];
sx q[0];
rz(-1.845813) q[0];
sx q[0];
rz(0.14607231) q[0];
rz(-pi) q[1];
x q[1];
rz(2.550225) q[2];
sx q[2];
rz(-0.57833507) q[2];
sx q[2];
rz(-0.27507281) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27998566) q[1];
sx q[1];
rz(-1.0000739) q[1];
sx q[1];
rz(3.0793889) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60243209) q[3];
sx q[3];
rz(-0.89655399) q[3];
sx q[3];
rz(2.7979421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9810527) q[2];
sx q[2];
rz(-1.2489698) q[2];
sx q[2];
rz(-0.36671656) q[2];
rz(1.8803053) q[3];
sx q[3];
rz(-1.4533318) q[3];
sx q[3];
rz(-3.0453851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325608) q[0];
sx q[0];
rz(-0.92027396) q[0];
sx q[0];
rz(-2.5352056) q[0];
rz(2.9442893) q[1];
sx q[1];
rz(-2.0154672) q[1];
sx q[1];
rz(2.6775449) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1394135) q[0];
sx q[0];
rz(-2.3521949) q[0];
sx q[0];
rz(2.1887357) q[0];
rz(3.1214141) q[2];
sx q[2];
rz(-1.8115461) q[2];
sx q[2];
rz(1.9764331) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11394994) q[1];
sx q[1];
rz(-1.0877891) q[1];
sx q[1];
rz(0.34160683) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8125698) q[3];
sx q[3];
rz(-2.892422) q[3];
sx q[3];
rz(0.84013018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2074034) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(-0.25804538) q[2];
rz(-1.1856273) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(-3.1380838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.946452) q[0];
sx q[0];
rz(-1.2807245) q[0];
sx q[0];
rz(-0.38129693) q[0];
rz(-0.095245846) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(-1.7000748) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38948108) q[0];
sx q[0];
rz(-1.507326) q[0];
sx q[0];
rz(-1.5556637) q[0];
rz(-pi) q[1];
rz(-1.5324138) q[2];
sx q[2];
rz(-2.1830609) q[2];
sx q[2];
rz(1.6910764) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0563593) q[1];
sx q[1];
rz(-1.4409522) q[1];
sx q[1];
rz(1.6972046) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0184047) q[3];
sx q[3];
rz(-2.8052969) q[3];
sx q[3];
rz(1.5817225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9986481) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(2.9150035) q[2];
rz(2.6930124) q[3];
sx q[3];
rz(-1.5357163) q[3];
sx q[3];
rz(-2.3118238) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7609693) q[0];
sx q[0];
rz(-2.3801104) q[0];
sx q[0];
rz(1.7425849) q[0];
rz(-2.8245068) q[1];
sx q[1];
rz(-1.6665019) q[1];
sx q[1];
rz(0.98659602) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6331659) q[0];
sx q[0];
rz(-1.0796483) q[0];
sx q[0];
rz(-1.7500061) q[0];
rz(-pi) q[1];
rz(-0.53194745) q[2];
sx q[2];
rz(-0.94594687) q[2];
sx q[2];
rz(2.9192386) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0341067) q[1];
sx q[1];
rz(-1.7120546) q[1];
sx q[1];
rz(-0.60955255) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9429728) q[3];
sx q[3];
rz(-2.6937006) q[3];
sx q[3];
rz(0.4263634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2150779) q[2];
sx q[2];
rz(-2.4145917) q[2];
sx q[2];
rz(-2.731936) q[2];
rz(2.8783197) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(0.51945654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39919329) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(1.7364527) q[0];
rz(2.3204904) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(-1.4155037) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1499407) q[0];
sx q[0];
rz(-1.2399925) q[0];
sx q[0];
rz(-1.8206235) q[0];
rz(-pi) q[1];
rz(0.24869463) q[2];
sx q[2];
rz(-1.9242312) q[2];
sx q[2];
rz(-2.8949708) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7064757) q[1];
sx q[1];
rz(-2.8787328) q[1];
sx q[1];
rz(0.95107066) q[1];
rz(-pi) q[2];
rz(-2.0831765) q[3];
sx q[3];
rz(-1.7184162) q[3];
sx q[3];
rz(-2.3073713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4225509) q[2];
sx q[2];
rz(-0.30095235) q[2];
sx q[2];
rz(0.12410513) q[2];
rz(-2.1758046) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(2.4805099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.538095) q[0];
sx q[0];
rz(-2.8932543) q[0];
sx q[0];
rz(2.2809991) q[0];
rz(2.8339236) q[1];
sx q[1];
rz(-1.2528906) q[1];
sx q[1];
rz(1.2045592) q[1];
rz(2.8812257) q[2];
sx q[2];
rz(-1.7454864) q[2];
sx q[2];
rz(-2.7228552) q[2];
rz(0.55660558) q[3];
sx q[3];
rz(-1.6986596) q[3];
sx q[3];
rz(-1.2996775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
