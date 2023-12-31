OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2330115) q[0];
sx q[0];
rz(-1.1865948) q[0];
sx q[0];
rz(-1.9957805) q[0];
rz(0.97283483) q[1];
sx q[1];
rz(-1.4714779) q[1];
sx q[1];
rz(0.29247984) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4300633) q[0];
sx q[0];
rz(-1.9646514) q[0];
sx q[0];
rz(-0.060056134) q[0];
rz(-pi) q[1];
rz(-3.0897065) q[2];
sx q[2];
rz(-2.3615026) q[2];
sx q[2];
rz(0.29733301) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3136183) q[1];
sx q[1];
rz(-2.5458114) q[1];
sx q[1];
rz(1.4890563) q[1];
rz(-pi) q[2];
rz(-1.0969767) q[3];
sx q[3];
rz(-1.5690656) q[3];
sx q[3];
rz(-0.8918106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4108489) q[2];
sx q[2];
rz(-2.0480672) q[2];
sx q[2];
rz(0.4494108) q[2];
rz(0.64569008) q[3];
sx q[3];
rz(-2.460545) q[3];
sx q[3];
rz(1.2908363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21355024) q[0];
sx q[0];
rz(-0.88590652) q[0];
sx q[0];
rz(-2.399562) q[0];
rz(1.4713326) q[1];
sx q[1];
rz(-0.69683087) q[1];
sx q[1];
rz(-1.3630294) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8105158) q[0];
sx q[0];
rz(-1.1805981) q[0];
sx q[0];
rz(-0.66731989) q[0];
rz(-0.46911932) q[2];
sx q[2];
rz(-1.763478) q[2];
sx q[2];
rz(2.0872781) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0568309) q[1];
sx q[1];
rz(-2.8014604) q[1];
sx q[1];
rz(0.057573307) q[1];
rz(-pi) q[2];
rz(0.19248776) q[3];
sx q[3];
rz(-1.3093595) q[3];
sx q[3];
rz(0.9769494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8403975) q[2];
sx q[2];
rz(-0.49673721) q[2];
sx q[2];
rz(1.2724686) q[2];
rz(2.8043591) q[3];
sx q[3];
rz(-1.971258) q[3];
sx q[3];
rz(-3.1315394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067588016) q[0];
sx q[0];
rz(-2.8019866) q[0];
sx q[0];
rz(-0.2628251) q[0];
rz(2.7858531) q[1];
sx q[1];
rz(-1.4590615) q[1];
sx q[1];
rz(1.9062818) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4848029) q[0];
sx q[0];
rz(-0.018253837) q[0];
sx q[0];
rz(2.546893) q[0];
rz(1.8112196) q[2];
sx q[2];
rz(-1.5449107) q[2];
sx q[2];
rz(-2.7607069) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0318109) q[1];
sx q[1];
rz(-2.8019252) q[1];
sx q[1];
rz(-0.32082816) q[1];
rz(-pi) q[2];
rz(1.1615109) q[3];
sx q[3];
rz(-0.73771362) q[3];
sx q[3];
rz(-1.1324594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1245023) q[2];
sx q[2];
rz(-1.7586781) q[2];
sx q[2];
rz(0.18033218) q[2];
rz(-0.63594937) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(2.7699871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3746049) q[0];
sx q[0];
rz(-0.42414442) q[0];
sx q[0];
rz(2.4225127) q[0];
rz(-2.6285697) q[1];
sx q[1];
rz(-1.2027556) q[1];
sx q[1];
rz(1.0346574) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9229139) q[0];
sx q[0];
rz(-1.551034) q[0];
sx q[0];
rz(-3.0826871) q[0];
x q[1];
rz(-1.2345418) q[2];
sx q[2];
rz(-0.14306919) q[2];
sx q[2];
rz(-2.2640995) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3116341) q[1];
sx q[1];
rz(-1.8533851) q[1];
sx q[1];
rz(-2.4294873) q[1];
rz(1.3961117) q[3];
sx q[3];
rz(-1.1038728) q[3];
sx q[3];
rz(-0.13388453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.93306142) q[2];
sx q[2];
rz(-2.8716817) q[2];
sx q[2];
rz(2.9053524) q[2];
rz(-1.158372) q[3];
sx q[3];
rz(-1.002243) q[3];
sx q[3];
rz(-0.85197824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(3.1119969) q[0];
sx q[0];
rz(-1.010226) q[0];
sx q[0];
rz(0.90233666) q[0];
rz(-0.92102712) q[1];
sx q[1];
rz(-0.61704707) q[1];
sx q[1];
rz(-2.5193118) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.92684) q[0];
sx q[0];
rz(-1.2878969) q[0];
sx q[0];
rz(-2.9515285) q[0];
rz(-1.6629409) q[2];
sx q[2];
rz(-2.2170057) q[2];
sx q[2];
rz(2.9421633) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9363074) q[1];
sx q[1];
rz(-2.949932) q[1];
sx q[1];
rz(-1.081341) q[1];
rz(-0.47982521) q[3];
sx q[3];
rz(-1.2180426) q[3];
sx q[3];
rz(-1.6640116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92419147) q[2];
sx q[2];
rz(-1.3240341) q[2];
sx q[2];
rz(3.1203111) q[2];
rz(-1.4422669) q[3];
sx q[3];
rz(-0.71443021) q[3];
sx q[3];
rz(0.91059476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845881) q[0];
sx q[0];
rz(-0.5963043) q[0];
sx q[0];
rz(0.044629991) q[0];
rz(2.1394829) q[1];
sx q[1];
rz(-1.0005181) q[1];
sx q[1];
rz(0.034428509) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39715595) q[0];
sx q[0];
rz(-1.7560871) q[0];
sx q[0];
rz(0.013884355) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0677412) q[2];
sx q[2];
rz(-1.4639374) q[2];
sx q[2];
rz(1.6018794) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.039636314) q[1];
sx q[1];
rz(-0.64412457) q[1];
sx q[1];
rz(-2.1307751) q[1];
rz(-pi) q[2];
rz(0.19272007) q[3];
sx q[3];
rz(-0.7145213) q[3];
sx q[3];
rz(-2.2477828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.78952152) q[2];
sx q[2];
rz(-1.8692724) q[2];
sx q[2];
rz(-0.96674031) q[2];
rz(0.012185193) q[3];
sx q[3];
rz(-0.92958486) q[3];
sx q[3];
rz(-0.10908443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1352585) q[0];
sx q[0];
rz(-2.874458) q[0];
sx q[0];
rz(2.1475041) q[0];
rz(-3.0864691) q[1];
sx q[1];
rz(-1.2722641) q[1];
sx q[1];
rz(-0.73928839) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2544884) q[0];
sx q[0];
rz(-1.6657889) q[0];
sx q[0];
rz(1.4332921) q[0];
x q[1];
rz(2.5653966) q[2];
sx q[2];
rz(-1.9857166) q[2];
sx q[2];
rz(2.8232099) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1623605) q[1];
sx q[1];
rz(-1.3812961) q[1];
sx q[1];
rz(-2.6081671) q[1];
rz(-pi) q[2];
rz(-1.2823295) q[3];
sx q[3];
rz(-0.97086421) q[3];
sx q[3];
rz(-2.695042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.58549515) q[2];
sx q[2];
rz(-1.7332417) q[2];
sx q[2];
rz(0.56376702) q[2];
rz(0.24027763) q[3];
sx q[3];
rz(-0.53301817) q[3];
sx q[3];
rz(-1.7668004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69650841) q[0];
sx q[0];
rz(-2.962528) q[0];
sx q[0];
rz(-0.58404303) q[0];
rz(-0.42075992) q[1];
sx q[1];
rz(-1.0412443) q[1];
sx q[1];
rz(-2.8631794) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7137372) q[0];
sx q[0];
rz(-2.5579778) q[0];
sx q[0];
rz(-1.7996656) q[0];
rz(0.39324795) q[2];
sx q[2];
rz(-1.3967782) q[2];
sx q[2];
rz(-2.4207123) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0502187) q[1];
sx q[1];
rz(-2.723657) q[1];
sx q[1];
rz(-2.1347079) q[1];
rz(-pi) q[2];
rz(-1.1823468) q[3];
sx q[3];
rz(-2.5591345) q[3];
sx q[3];
rz(2.9782481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.86928308) q[2];
sx q[2];
rz(-0.81988207) q[2];
sx q[2];
rz(0.38796866) q[2];
rz(-1.3502454) q[3];
sx q[3];
rz(-0.46348059) q[3];
sx q[3];
rz(0.2750245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824654) q[0];
sx q[0];
rz(-0.21331856) q[0];
sx q[0];
rz(-0.7318837) q[0];
rz(-0.16648079) q[1];
sx q[1];
rz(-0.44547588) q[1];
sx q[1];
rz(3.0678715) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6442935) q[0];
sx q[0];
rz(-1.6719581) q[0];
sx q[0];
rz(-2.0072719) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97754064) q[2];
sx q[2];
rz(-2.9974555) q[2];
sx q[2];
rz(1.6050715) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8718865) q[1];
sx q[1];
rz(-1.6082113) q[1];
sx q[1];
rz(-0.58362959) q[1];
rz(-pi) q[2];
rz(2.5684486) q[3];
sx q[3];
rz(-1.4775652) q[3];
sx q[3];
rz(-1.7236818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5471197) q[2];
sx q[2];
rz(-0.23590817) q[2];
sx q[2];
rz(2.7129042) q[2];
rz(1.3171014) q[3];
sx q[3];
rz(-1.3573815) q[3];
sx q[3];
rz(-1.930442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4093032) q[0];
sx q[0];
rz(-1.4855054) q[0];
sx q[0];
rz(2.0987341) q[0];
rz(1.5111142) q[1];
sx q[1];
rz(-1.7621367) q[1];
sx q[1];
rz(1.7932549) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74603426) q[0];
sx q[0];
rz(-1.6840044) q[0];
sx q[0];
rz(3.1157225) q[0];
x q[1];
rz(1.4820319) q[2];
sx q[2];
rz(-2.4030493) q[2];
sx q[2];
rz(-0.86631394) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.18689449) q[1];
sx q[1];
rz(-2.8098626) q[1];
sx q[1];
rz(1.7897254) q[1];
rz(-2.9006349) q[3];
sx q[3];
rz(-1.3853419) q[3];
sx q[3];
rz(-2.0897739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3412791) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(3.0467765) q[2];
rz(-1.9684277) q[3];
sx q[3];
rz(-1.4011819) q[3];
sx q[3];
rz(-0.15869424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6288347) q[0];
sx q[0];
rz(-0.47369581) q[0];
sx q[0];
rz(1.0439903) q[0];
rz(-1.5402773) q[1];
sx q[1];
rz(-1.5834783) q[1];
sx q[1];
rz(-1.6316354) q[1];
rz(1.9376471) q[2];
sx q[2];
rz(-0.22782142) q[2];
sx q[2];
rz(0.08624764) q[2];
rz(0.52900984) q[3];
sx q[3];
rz(-0.88376868) q[3];
sx q[3];
rz(0.39467011) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
