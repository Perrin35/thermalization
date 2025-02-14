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
rz(1.709197) q[0];
sx q[0];
rz(3.4442918) q[0];
sx q[0];
rz(11.533307) q[0];
rz(1.9167702) q[1];
sx q[1];
rz(-1.6515825) q[1];
sx q[1];
rz(0.19436793) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2107735) q[0];
sx q[0];
rz(-1.418485) q[0];
sx q[0];
rz(2.5126626) q[0];
rz(-pi) q[1];
rz(-0.35440234) q[2];
sx q[2];
rz(-0.72847073) q[2];
sx q[2];
rz(-0.013106339) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3981975) q[1];
sx q[1];
rz(-1.6170073) q[1];
sx q[1];
rz(0.034548918) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0048054) q[3];
sx q[3];
rz(-0.99054407) q[3];
sx q[3];
rz(-2.2541719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4091461) q[2];
sx q[2];
rz(-2.0398085) q[2];
sx q[2];
rz(-0.5644325) q[2];
rz(1.8730646) q[3];
sx q[3];
rz(-3.1334435) q[3];
sx q[3];
rz(2.234999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0821575) q[0];
sx q[0];
rz(-3.1340288) q[0];
sx q[0];
rz(0.91307688) q[0];
rz(-0.68436855) q[1];
sx q[1];
rz(-0.00033907779) q[1];
sx q[1];
rz(-2.2025462) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8565687) q[0];
sx q[0];
rz(-1.9459581) q[0];
sx q[0];
rz(0.60174353) q[0];
x q[1];
rz(0.46397658) q[2];
sx q[2];
rz(-1.5590406) q[2];
sx q[2];
rz(-1.3323615) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9833656) q[1];
sx q[1];
rz(-1.2364796) q[1];
sx q[1];
rz(2.0086221) q[1];
x q[2];
rz(-3.0747218) q[3];
sx q[3];
rz(-1.9767728) q[3];
sx q[3];
rz(1.2479046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3188476) q[2];
sx q[2];
rz(-0.007195909) q[2];
sx q[2];
rz(2.2491271) q[2];
rz(1.578963) q[3];
sx q[3];
rz(-0.024450863) q[3];
sx q[3];
rz(-1.8306556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40176323) q[0];
sx q[0];
rz(-0.50245291) q[0];
sx q[0];
rz(-2.7318562) q[0];
rz(-0.01092625) q[1];
sx q[1];
rz(-2.9210563) q[1];
sx q[1];
rz(-1.797537) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4554268) q[0];
sx q[0];
rz(-0.51138216) q[0];
sx q[0];
rz(-1.6126942) q[0];
x q[1];
rz(0.063882752) q[2];
sx q[2];
rz(-1.4555571) q[2];
sx q[2];
rz(1.921953) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2933767) q[1];
sx q[1];
rz(-2.8753706) q[1];
sx q[1];
rz(-1.3220399) q[1];
x q[2];
rz(1.4555172) q[3];
sx q[3];
rz(-1.56101) q[3];
sx q[3];
rz(-0.20231314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.672013) q[2];
sx q[2];
rz(-1.6894222) q[2];
sx q[2];
rz(-0.043188485) q[2];
rz(2.0016661) q[3];
sx q[3];
rz(-2.98525) q[3];
sx q[3];
rz(1.202701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4129055) q[0];
sx q[0];
rz(-0.40189704) q[0];
sx q[0];
rz(-1.3380949) q[0];
rz(-0.69522229) q[1];
sx q[1];
rz(-3.0137364) q[1];
sx q[1];
rz(-2.7241838) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.886717) q[0];
sx q[0];
rz(-2.4435197) q[0];
sx q[0];
rz(1.8398192) q[0];
rz(-pi) q[1];
rz(-1.2222671) q[2];
sx q[2];
rz(-2.4651595) q[2];
sx q[2];
rz(1.8944905) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1811235) q[1];
sx q[1];
rz(-2.4669018) q[1];
sx q[1];
rz(1.2441649) q[1];
rz(-pi) q[2];
rz(2.999971) q[3];
sx q[3];
rz(-1.1417665) q[3];
sx q[3];
rz(0.98875586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7030316) q[2];
sx q[2];
rz(-2.265354) q[2];
sx q[2];
rz(1.38928) q[2];
rz(-2.321068) q[3];
sx q[3];
rz(-1.5503784) q[3];
sx q[3];
rz(0.70782053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97838068) q[0];
sx q[0];
rz(-3.1040525) q[0];
sx q[0];
rz(-2.1782844) q[0];
rz(-2.9523201) q[1];
sx q[1];
rz(-3.1258588) q[1];
sx q[1];
rz(0.16545573) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67901299) q[0];
sx q[0];
rz(-1.5272584) q[0];
sx q[0];
rz(-1.5444368) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2654183) q[2];
sx q[2];
rz(-1.3378222) q[2];
sx q[2];
rz(0.35859749) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.206659) q[1];
sx q[1];
rz(-2.3689785) q[1];
sx q[1];
rz(1.7375224) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1894375) q[3];
sx q[3];
rz(-0.51636558) q[3];
sx q[3];
rz(-1.0960033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8339707) q[2];
sx q[2];
rz(-0.51887363) q[2];
sx q[2];
rz(-2.0818254) q[2];
rz(2.4506532) q[3];
sx q[3];
rz(-0.86425224) q[3];
sx q[3];
rz(1.8701657) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6178013) q[0];
sx q[0];
rz(-0.093136223) q[0];
sx q[0];
rz(1.5366489) q[0];
rz(-0.45283428) q[1];
sx q[1];
rz(-3.1335399) q[1];
sx q[1];
rz(-1.7013928) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0217845) q[0];
sx q[0];
rz(-0.23384604) q[0];
sx q[0];
rz(0.77063693) q[0];
rz(-pi) q[1];
rz(0.78078713) q[2];
sx q[2];
rz(-2.8344354) q[2];
sx q[2];
rz(3.0083092) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3611412) q[1];
sx q[1];
rz(-1.5634426) q[1];
sx q[1];
rz(0.15463923) q[1];
rz(-pi) q[2];
rz(-3.0423445) q[3];
sx q[3];
rz(-1.7407047) q[3];
sx q[3];
rz(-0.75796321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77233934) q[2];
sx q[2];
rz(-2.1489096) q[2];
sx q[2];
rz(-0.079027979) q[2];
rz(2.2760271) q[3];
sx q[3];
rz(-0.95439684) q[3];
sx q[3];
rz(1.0237833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2257776) q[0];
sx q[0];
rz(-0.0053891698) q[0];
sx q[0];
rz(0.22501568) q[0];
rz(2.8354538) q[1];
sx q[1];
rz(-3.1248416) q[1];
sx q[1];
rz(-0.87047815) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9266548) q[0];
sx q[0];
rz(-1.4995725) q[0];
sx q[0];
rz(-1.5230973) q[0];
x q[1];
rz(-0.71484675) q[2];
sx q[2];
rz(-2.0185592) q[2];
sx q[2];
rz(0.078850672) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6806769) q[1];
sx q[1];
rz(-0.80246882) q[1];
sx q[1];
rz(-1.4867307) q[1];
rz(-pi) q[2];
rz(2.7288658) q[3];
sx q[3];
rz(-2.0214635) q[3];
sx q[3];
rz(-0.13201763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47591448) q[2];
sx q[2];
rz(-3.0147538) q[2];
sx q[2];
rz(1.0352943) q[2];
rz(0.78694844) q[3];
sx q[3];
rz(-3.0988155) q[3];
sx q[3];
rz(-0.27534819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.03054522) q[0];
sx q[0];
rz(-3.1304517) q[0];
sx q[0];
rz(0.025644843) q[0];
rz(-1.8772839) q[1];
sx q[1];
rz(-0.023921078) q[1];
sx q[1];
rz(-2.4425676) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5830987) q[0];
sx q[0];
rz(-1.8374658) q[0];
sx q[0];
rz(1.3400274) q[0];
rz(-pi) q[1];
rz(0.46866663) q[2];
sx q[2];
rz(-1.7231517) q[2];
sx q[2];
rz(3.0909757) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0009202) q[1];
sx q[1];
rz(-0.39127054) q[1];
sx q[1];
rz(-2.9468582) q[1];
rz(-2.5949536) q[3];
sx q[3];
rz(-2.1501503) q[3];
sx q[3];
rz(2.1648615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.25253025) q[2];
sx q[2];
rz(-2.3928596) q[2];
sx q[2];
rz(-2.8177281) q[2];
rz(0.1570541) q[3];
sx q[3];
rz(-1.2374977) q[3];
sx q[3];
rz(-2.9878555) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57482982) q[0];
sx q[0];
rz(-0.014641849) q[0];
sx q[0];
rz(-0.59004849) q[0];
rz(-0.7198965) q[1];
sx q[1];
rz(-3.0820334) q[1];
sx q[1];
rz(0.86729008) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27408263) q[0];
sx q[0];
rz(-1.6757586) q[0];
sx q[0];
rz(-0.63469736) q[0];
rz(-1.2008823) q[2];
sx q[2];
rz(-0.6354161) q[2];
sx q[2];
rz(-0.70178343) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0194975) q[1];
sx q[1];
rz(-1.6830793) q[1];
sx q[1];
rz(-3.0730559) q[1];
rz(-pi) q[2];
rz(-1.7884364) q[3];
sx q[3];
rz(-3*pi/13) q[3];
sx q[3];
rz(-2.9293037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94816339) q[2];
sx q[2];
rz(-2.2488504) q[2];
sx q[2];
rz(-0.11290045) q[2];
rz(-2.0924977) q[3];
sx q[3];
rz(-1.5821404) q[3];
sx q[3];
rz(-0.73672867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.076544) q[0];
sx q[0];
rz(-1.5815409) q[0];
sx q[0];
rz(1.5034058) q[0];
rz(-0.63649559) q[1];
sx q[1];
rz(-0.46000767) q[1];
sx q[1];
rz(1.5719302) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8419452) q[0];
sx q[0];
rz(-1.491786) q[0];
sx q[0];
rz(-1.6539198) q[0];
rz(-pi) q[1];
rz(3.1384597) q[2];
sx q[2];
rz(-1.5701446) q[2];
sx q[2];
rz(2.3359483) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1366982) q[1];
sx q[1];
rz(-1.5678798) q[1];
sx q[1];
rz(-1.5705622) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7336373) q[3];
sx q[3];
rz(-1.0615088) q[3];
sx q[3];
rz(-0.41679888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2703209) q[2];
sx q[2];
rz(-1.511829) q[2];
sx q[2];
rz(2.9809269) q[2];
rz(-1.2284944) q[3];
sx q[3];
rz(-0.03385032) q[3];
sx q[3];
rz(-0.20139774) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0495618) q[0];
sx q[0];
rz(-1.1783538) q[0];
sx q[0];
rz(0.023068064) q[0];
rz(-1.5184215) q[1];
sx q[1];
rz(-0.36133125) q[1];
sx q[1];
rz(0.28892118) q[1];
rz(-3.1186947) q[2];
sx q[2];
rz(-1.2780634) q[2];
sx q[2];
rz(-2.7414049) q[2];
rz(2.0102228) q[3];
sx q[3];
rz(-1.4356339) q[3];
sx q[3];
rz(2.1826134) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
