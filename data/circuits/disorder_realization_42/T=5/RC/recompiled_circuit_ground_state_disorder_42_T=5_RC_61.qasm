OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8920933) q[0];
sx q[0];
rz(-1.1846932) q[0];
sx q[0];
rz(2.0062334) q[0];
rz(2.7197977) q[1];
sx q[1];
rz(-0.78039688) q[1];
sx q[1];
rz(0.62224046) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6506158) q[0];
sx q[0];
rz(-0.19189206) q[0];
sx q[0];
rz(-2.2723115) q[0];
rz(-pi) q[1];
x q[1];
rz(0.038921629) q[2];
sx q[2];
rz(-2.0177671) q[2];
sx q[2];
rz(1.0124504) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.92889228) q[1];
sx q[1];
rz(-1.5824741) q[1];
sx q[1];
rz(2.547193) q[1];
rz(-pi) q[2];
rz(-2.9184266) q[3];
sx q[3];
rz(-1.0112178) q[3];
sx q[3];
rz(2.1055438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8047831) q[2];
sx q[2];
rz(-1.6853036) q[2];
sx q[2];
rz(-0.08610227) q[2];
rz(-0.51186776) q[3];
sx q[3];
rz(-2.0972926) q[3];
sx q[3];
rz(2.4477203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7030199) q[0];
sx q[0];
rz(-0.95412552) q[0];
sx q[0];
rz(2.2254206) q[0];
rz(-2.84962) q[1];
sx q[1];
rz(-2.5501854) q[1];
sx q[1];
rz(0.77436647) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.403023) q[0];
sx q[0];
rz(-2.6576808) q[0];
sx q[0];
rz(-0.83173521) q[0];
x q[1];
rz(-1.5874812) q[2];
sx q[2];
rz(-2.7962001) q[2];
sx q[2];
rz(0.85341233) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9271597) q[1];
sx q[1];
rz(-1.4478874) q[1];
sx q[1];
rz(1.4686119) q[1];
x q[2];
rz(-0.99454576) q[3];
sx q[3];
rz(-2.3062097) q[3];
sx q[3];
rz(1.6540838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.795221) q[2];
sx q[2];
rz(-0.56196153) q[2];
sx q[2];
rz(-0.58491659) q[2];
rz(1.6631205) q[3];
sx q[3];
rz(-1.9641179) q[3];
sx q[3];
rz(1.4896721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65603489) q[0];
sx q[0];
rz(-1.0233044) q[0];
sx q[0];
rz(-2.5265332) q[0];
rz(-2.8133605) q[1];
sx q[1];
rz(-2.131772) q[1];
sx q[1];
rz(-2.097791) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39343483) q[0];
sx q[0];
rz(-1.6232449) q[0];
sx q[0];
rz(1.5657052) q[0];
x q[1];
rz(-1.5756025) q[2];
sx q[2];
rz(-1.9991181) q[2];
sx q[2];
rz(2.7132963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.17576007) q[1];
sx q[1];
rz(-1.1783667) q[1];
sx q[1];
rz(-2.2821941) q[1];
rz(-pi) q[2];
rz(1.2112593) q[3];
sx q[3];
rz(-1.94424) q[3];
sx q[3];
rz(-1.5869753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.388776) q[2];
sx q[2];
rz(-2.7478168) q[2];
sx q[2];
rz(-2.8326995) q[2];
rz(-2.654352) q[3];
sx q[3];
rz(-1.4332708) q[3];
sx q[3];
rz(-0.88700956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3798645) q[0];
sx q[0];
rz(-0.74493113) q[0];
sx q[0];
rz(-2.6640025) q[0];
rz(1.8567122) q[1];
sx q[1];
rz(-1.4297337) q[1];
sx q[1];
rz(0.46599785) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6687688) q[0];
sx q[0];
rz(-1.6472938) q[0];
sx q[0];
rz(-0.55474218) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1144756) q[2];
sx q[2];
rz(-1.9812036) q[2];
sx q[2];
rz(2.9027651) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6324959) q[1];
sx q[1];
rz(-1.5358616) q[1];
sx q[1];
rz(0.81991244) q[1];
rz(2.4336045) q[3];
sx q[3];
rz(-2.3205415) q[3];
sx q[3];
rz(-0.94828965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3232702) q[2];
sx q[2];
rz(-2.6939836) q[2];
sx q[2];
rz(1.4480048) q[2];
rz(-1.4130392) q[3];
sx q[3];
rz(-1.5331242) q[3];
sx q[3];
rz(-1.9938699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84561658) q[0];
sx q[0];
rz(-2.4281261) q[0];
sx q[0];
rz(0.071685858) q[0];
rz(-1.4777615) q[1];
sx q[1];
rz(-1.3641337) q[1];
sx q[1];
rz(-0.62612265) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6090995) q[0];
sx q[0];
rz(-0.83864826) q[0];
sx q[0];
rz(-3.0200543) q[0];
rz(-1.4914054) q[2];
sx q[2];
rz(-2.491386) q[2];
sx q[2];
rz(-2.2967867) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4461527) q[1];
sx q[1];
rz(-1.3976516) q[1];
sx q[1];
rz(-2.3320564) q[1];
x q[2];
rz(-1.7528698) q[3];
sx q[3];
rz(-1.9183049) q[3];
sx q[3];
rz(3.0680498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6757297) q[2];
sx q[2];
rz(-0.53368038) q[2];
sx q[2];
rz(1.6443058) q[2];
rz(-0.40694445) q[3];
sx q[3];
rz(-0.99138433) q[3];
sx q[3];
rz(-2.2814894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1740455) q[0];
sx q[0];
rz(-1.8438735) q[0];
sx q[0];
rz(-2.8795854) q[0];
rz(1.920248) q[1];
sx q[1];
rz(-1.5120993) q[1];
sx q[1];
rz(1.2006203) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9438433) q[0];
sx q[0];
rz(-2.245994) q[0];
sx q[0];
rz(-0.64267107) q[0];
rz(-pi) q[1];
rz(-0.85882218) q[2];
sx q[2];
rz(-1.1505652) q[2];
sx q[2];
rz(2.4855297) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9900254) q[1];
sx q[1];
rz(-2.5737816) q[1];
sx q[1];
rz(2.0358596) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2030199) q[3];
sx q[3];
rz(-3.0431755) q[3];
sx q[3];
rz(-0.84151387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0210586) q[2];
sx q[2];
rz(-2.6594682) q[2];
sx q[2];
rz(-0.47312197) q[2];
rz(1.037723) q[3];
sx q[3];
rz(-0.40268746) q[3];
sx q[3];
rz(-1.13824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6556019) q[0];
sx q[0];
rz(-2.7779873) q[0];
sx q[0];
rz(2.3714491) q[0];
rz(-0.4153525) q[1];
sx q[1];
rz(-1.5130006) q[1];
sx q[1];
rz(1.8619246) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3949385) q[0];
sx q[0];
rz(-1.5956559) q[0];
sx q[0];
rz(-0.8381054) q[0];
x q[1];
rz(1.1392904) q[2];
sx q[2];
rz(-1.4551292) q[2];
sx q[2];
rz(2.3097599) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.61381147) q[1];
sx q[1];
rz(-1.4134058) q[1];
sx q[1];
rz(1.792683) q[1];
rz(3.123056) q[3];
sx q[3];
rz(-2.662475) q[3];
sx q[3];
rz(-1.9974513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8659849) q[2];
sx q[2];
rz(-1.7476247) q[2];
sx q[2];
rz(2.7371791) q[2];
rz(1.154493) q[3];
sx q[3];
rz(-1.6136074) q[3];
sx q[3];
rz(-1.9563458) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6605717) q[0];
sx q[0];
rz(-2.4206929) q[0];
sx q[0];
rz(0.50203669) q[0];
rz(-2.1227116) q[1];
sx q[1];
rz(-1.6682245) q[1];
sx q[1];
rz(-0.92207164) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4785504) q[0];
sx q[0];
rz(-1.6660569) q[0];
sx q[0];
rz(0.24597286) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8174632) q[2];
sx q[2];
rz(-0.96117678) q[2];
sx q[2];
rz(2.7423046) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3193911) q[1];
sx q[1];
rz(-1.5272045) q[1];
sx q[1];
rz(1.6258214) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4430674) q[3];
sx q[3];
rz(-1.5683577) q[3];
sx q[3];
rz(-1.7170441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.33793882) q[2];
sx q[2];
rz(-1.2793845) q[2];
sx q[2];
rz(-3.0214018) q[2];
rz(3.1386612) q[3];
sx q[3];
rz(-2.7836697) q[3];
sx q[3];
rz(0.005793747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1015162) q[0];
sx q[0];
rz(-0.67779556) q[0];
sx q[0];
rz(-1.6325604) q[0];
rz(2.312233) q[1];
sx q[1];
rz(-0.50193915) q[1];
sx q[1];
rz(2.0109743) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.254758) q[0];
sx q[0];
rz(-1.8811329) q[0];
sx q[0];
rz(-1.929105) q[0];
x q[1];
rz(-2.0810602) q[2];
sx q[2];
rz(-0.9212554) q[2];
sx q[2];
rz(-0.559597) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3106416) q[1];
sx q[1];
rz(-2.0071908) q[1];
sx q[1];
rz(2.2121625) q[1];
rz(-pi) q[2];
rz(-1.0535766) q[3];
sx q[3];
rz(-1.791009) q[3];
sx q[3];
rz(2.6802879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.074038) q[2];
sx q[2];
rz(-0.87104565) q[2];
sx q[2];
rz(1.3646431) q[2];
rz(3.126295) q[3];
sx q[3];
rz(-0.8316032) q[3];
sx q[3];
rz(1.8196222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44052723) q[0];
sx q[0];
rz(-2.4897713) q[0];
sx q[0];
rz(-0.1524674) q[0];
rz(1.7561779) q[1];
sx q[1];
rz(-2.0962174) q[1];
sx q[1];
rz(-0.32858953) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6406031) q[0];
sx q[0];
rz(-3.1251934) q[0];
sx q[0];
rz(2.5522405) q[0];
rz(1.8118906) q[2];
sx q[2];
rz(-0.82359353) q[2];
sx q[2];
rz(-1.4581949) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6870881) q[1];
sx q[1];
rz(-1.3326613) q[1];
sx q[1];
rz(3.1109875) q[1];
rz(-pi) q[2];
rz(2.4963107) q[3];
sx q[3];
rz(-1.8580164) q[3];
sx q[3];
rz(0.20791277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3370328) q[2];
sx q[2];
rz(-2.0064662) q[2];
sx q[2];
rz(1.0767153) q[2];
rz(2.1870901) q[3];
sx q[3];
rz(-1.7007622) q[3];
sx q[3];
rz(-2.8452828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8148707) q[0];
sx q[0];
rz(-2.1537415) q[0];
sx q[0];
rz(1.4065773) q[0];
rz(1.8234491) q[1];
sx q[1];
rz(-1.6271918) q[1];
sx q[1];
rz(0.79961332) q[1];
rz(1.6257269) q[2];
sx q[2];
rz(-1.5274897) q[2];
sx q[2];
rz(-2.4804583) q[2];
rz(1.7794505) q[3];
sx q[3];
rz(-2.2310774) q[3];
sx q[3];
rz(3.0175573) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
