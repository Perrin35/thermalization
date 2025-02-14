OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.159654) q[0];
sx q[0];
rz(-1.5647793) q[0];
sx q[0];
rz(1.2793581) q[0];
rz(0.98969412) q[1];
sx q[1];
rz(5.0862105) q[1];
sx q[1];
rz(12.401019) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26583126) q[0];
sx q[0];
rz(-2.1935757) q[0];
sx q[0];
rz(2.7490007) q[0];
rz(-pi) q[1];
rz(-0.44423361) q[2];
sx q[2];
rz(-1.9529251) q[2];
sx q[2];
rz(1.5096779) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0667116) q[1];
sx q[1];
rz(-2.0600658) q[1];
sx q[1];
rz(2.679945) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8836423) q[3];
sx q[3];
rz(-1.2688302) q[3];
sx q[3];
rz(-2.0998968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2692261) q[2];
sx q[2];
rz(-1.9908315) q[2];
sx q[2];
rz(-1.9744138) q[2];
rz(0.40927467) q[3];
sx q[3];
rz(-0.27547488) q[3];
sx q[3];
rz(-2.1854713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.069365) q[0];
sx q[0];
rz(-0.15412155) q[0];
sx q[0];
rz(-1.7896205) q[0];
rz(-1.4714454) q[1];
sx q[1];
rz(-2.2189326) q[1];
sx q[1];
rz(-2.2460489) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.937718) q[0];
sx q[0];
rz(-0.26709891) q[0];
sx q[0];
rz(-0.0077203053) q[0];
rz(1.6487213) q[2];
sx q[2];
rz(-1.0977626) q[2];
sx q[2];
rz(0.89240197) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9546763) q[1];
sx q[1];
rz(-0.82789956) q[1];
sx q[1];
rz(2.1174341) q[1];
x q[2];
rz(1.9905521) q[3];
sx q[3];
rz(-1.7872056) q[3];
sx q[3];
rz(-1.2866502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6764549) q[2];
sx q[2];
rz(-0.82200161) q[2];
sx q[2];
rz(1.4466393) q[2];
rz(1.0195352) q[3];
sx q[3];
rz(-1.8229702) q[3];
sx q[3];
rz(2.8442966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4296221) q[0];
sx q[0];
rz(-1.1638887) q[0];
sx q[0];
rz(-0.0042560552) q[0];
rz(0.70107067) q[1];
sx q[1];
rz(-0.509976) q[1];
sx q[1];
rz(-0.027499011) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51873678) q[0];
sx q[0];
rz(-2.7122926) q[0];
sx q[0];
rz(1.1909498) q[0];
rz(-pi) q[1];
rz(-0.84553366) q[2];
sx q[2];
rz(-1.4730244) q[2];
sx q[2];
rz(0.95319437) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.356952) q[1];
sx q[1];
rz(-0.36423794) q[1];
sx q[1];
rz(-2.0924773) q[1];
rz(-pi) q[2];
rz(-1.6106538) q[3];
sx q[3];
rz(-0.33646001) q[3];
sx q[3];
rz(-2.5655059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6444401) q[2];
sx q[2];
rz(-1.7915598) q[2];
sx q[2];
rz(0.39786097) q[2];
rz(1.0128939) q[3];
sx q[3];
rz(-2.1161067) q[3];
sx q[3];
rz(-0.40867543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6560646) q[0];
sx q[0];
rz(-1.8714454) q[0];
sx q[0];
rz(-0.84554607) q[0];
rz(-0.0088508765) q[1];
sx q[1];
rz(-0.84819853) q[1];
sx q[1];
rz(1.1525851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6229079) q[0];
sx q[0];
rz(-0.9810514) q[0];
sx q[0];
rz(2.6481762) q[0];
rz(-pi) q[1];
rz(0.25817008) q[2];
sx q[2];
rz(-1.1328146) q[2];
sx q[2];
rz(-0.5817619) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2611672) q[1];
sx q[1];
rz(-0.5929872) q[1];
sx q[1];
rz(3.0914607) q[1];
x q[2];
rz(-2.827876) q[3];
sx q[3];
rz(-1.507618) q[3];
sx q[3];
rz(-2.3977195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3919966) q[2];
sx q[2];
rz(-2.1139202) q[2];
sx q[2];
rz(1.2006302) q[2];
rz(0.47731733) q[3];
sx q[3];
rz(-0.47003191) q[3];
sx q[3];
rz(-2.4376455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93477997) q[0];
sx q[0];
rz(-2.2581357) q[0];
sx q[0];
rz(-2.2373037) q[0];
rz(0.9710871) q[1];
sx q[1];
rz(-1.9458408) q[1];
sx q[1];
rz(2.8579874) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0440031) q[0];
sx q[0];
rz(-1.4999125) q[0];
sx q[0];
rz(-1.5547836) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7104811) q[2];
sx q[2];
rz(-1.1656802) q[2];
sx q[2];
rz(0.27027915) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7452104) q[1];
sx q[1];
rz(-1.3644427) q[1];
sx q[1];
rz(0.8861122) q[1];
rz(-2.141385) q[3];
sx q[3];
rz(-2.3688487) q[3];
sx q[3];
rz(-0.34602133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7859555) q[2];
sx q[2];
rz(-1.2355618) q[2];
sx q[2];
rz(2.9034485) q[2];
rz(-2.1613878) q[3];
sx q[3];
rz(-1.9832059) q[3];
sx q[3];
rz(2.3908206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7439483) q[0];
sx q[0];
rz(-1.605796) q[0];
sx q[0];
rz(-2.7342947) q[0];
rz(0.40335718) q[1];
sx q[1];
rz(-2.401001) q[1];
sx q[1];
rz(-2.2743646) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3639528) q[0];
sx q[0];
rz(-2.1356076) q[0];
sx q[0];
rz(3.0710247) q[0];
rz(2.5851493) q[2];
sx q[2];
rz(-1.3927937) q[2];
sx q[2];
rz(1.0484421) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7379511) q[1];
sx q[1];
rz(-0.80437901) q[1];
sx q[1];
rz(-2.7405095) q[1];
rz(1.9656455) q[3];
sx q[3];
rz(-0.32002745) q[3];
sx q[3];
rz(1.4744028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.12157089) q[2];
sx q[2];
rz(-2.3131804) q[2];
sx q[2];
rz(-1.7106445) q[2];
rz(-0.53267789) q[3];
sx q[3];
rz(-1.0360274) q[3];
sx q[3];
rz(1.7238341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
rz(1.1581887) q[0];
sx q[0];
rz(-0.15022763) q[0];
sx q[0];
rz(2.0576553) q[0];
rz(-0.051941959) q[1];
sx q[1];
rz(-0.40819326) q[1];
sx q[1];
rz(-3.1165677) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0753703) q[0];
sx q[0];
rz(-2.0522723) q[0];
sx q[0];
rz(-1.6027149) q[0];
rz(-1.0919326) q[2];
sx q[2];
rz(-2.3298878) q[2];
sx q[2];
rz(1.020249) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96341771) q[1];
sx q[1];
rz(-1.9994078) q[1];
sx q[1];
rz(-2.9531698) q[1];
rz(2.7152578) q[3];
sx q[3];
rz(-1.4407743) q[3];
sx q[3];
rz(-0.015817964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0673125) q[2];
sx q[2];
rz(-1.68579) q[2];
sx q[2];
rz(2.32453) q[2];
rz(2.1909511) q[3];
sx q[3];
rz(-1.0167511) q[3];
sx q[3];
rz(1.1869717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0364712) q[0];
sx q[0];
rz(-2.8558185) q[0];
sx q[0];
rz(-0.60086077) q[0];
rz(-0.034424456) q[1];
sx q[1];
rz(-1.333678) q[1];
sx q[1];
rz(1.6011802) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87266028) q[0];
sx q[0];
rz(-1.0300228) q[0];
sx q[0];
rz(-0.85178225) q[0];
rz(-2.41018) q[2];
sx q[2];
rz(-0.81899511) q[2];
sx q[2];
rz(2.7186944) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9084165) q[1];
sx q[1];
rz(-1.5926835) q[1];
sx q[1];
rz(-0.43969701) q[1];
rz(-pi) q[2];
rz(0.54564666) q[3];
sx q[3];
rz(-1.9522523) q[3];
sx q[3];
rz(1.1798525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2711266) q[2];
sx q[2];
rz(-1.8693962) q[2];
sx q[2];
rz(-1.4938483) q[2];
rz(2.3624524) q[3];
sx q[3];
rz(-1.6611764) q[3];
sx q[3];
rz(2.486865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43799022) q[0];
sx q[0];
rz(-2.0605189) q[0];
sx q[0];
rz(-0.78824747) q[0];
rz(-1.8986656) q[1];
sx q[1];
rz(-1.5595167) q[1];
sx q[1];
rz(2.8020614) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6166016) q[0];
sx q[0];
rz(-0.22540936) q[0];
sx q[0];
rz(-1.0269643) q[0];
rz(-pi) q[1];
rz(-1.2844769) q[2];
sx q[2];
rz(-1.0272046) q[2];
sx q[2];
rz(1.408615) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0362944) q[1];
sx q[1];
rz(-0.52043167) q[1];
sx q[1];
rz(-0.6639414) q[1];
rz(-pi) q[2];
rz(0.65776627) q[3];
sx q[3];
rz(-0.89221803) q[3];
sx q[3];
rz(2.984131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1303611) q[2];
sx q[2];
rz(-1.9731382) q[2];
sx q[2];
rz(1.2569024) q[2];
rz(-1.0907178) q[3];
sx q[3];
rz(-1.2033477) q[3];
sx q[3];
rz(-0.91712657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85482368) q[0];
sx q[0];
rz(-1.6463065) q[0];
sx q[0];
rz(2.0538034) q[0];
rz(-0.65025672) q[1];
sx q[1];
rz(-1.4573263) q[1];
sx q[1];
rz(0.021171721) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0567704) q[0];
sx q[0];
rz(-2.6426043) q[0];
sx q[0];
rz(-0.4468687) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6708605) q[2];
sx q[2];
rz(-1.7463049) q[2];
sx q[2];
rz(0.99711768) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6172006) q[1];
sx q[1];
rz(-0.78689303) q[1];
sx q[1];
rz(1.2715497) q[1];
rz(-pi) q[2];
rz(1.0383426) q[3];
sx q[3];
rz(-1.5335011) q[3];
sx q[3];
rz(-1.9017912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1374958) q[2];
sx q[2];
rz(-1.3429514) q[2];
sx q[2];
rz(-2.1583978) q[2];
rz(2.148597) q[3];
sx q[3];
rz(-2.0296622) q[3];
sx q[3];
rz(2.907311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8773593) q[0];
sx q[0];
rz(-2.3491884) q[0];
sx q[0];
rz(1.0607251) q[0];
rz(2.0296774) q[1];
sx q[1];
rz(-2.8363375) q[1];
sx q[1];
rz(-3.0189966) q[1];
rz(-1.1373043) q[2];
sx q[2];
rz(-0.52327427) q[2];
sx q[2];
rz(0.56868194) q[2];
rz(-2.6186416) q[3];
sx q[3];
rz(-1.5613489) q[3];
sx q[3];
rz(-0.47351826) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
