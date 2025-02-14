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
rz(-1.8622346) q[0];
rz(-2.1518985) q[1];
sx q[1];
rz(-1.9446179) q[1];
sx q[1];
rz(-2.9762414) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7900766) q[0];
sx q[0];
rz(-2.4195603) q[0];
sx q[0];
rz(-2.0603098) q[0];
rz(-pi) q[1];
rz(-1.1520391) q[2];
sx q[2];
rz(-1.9809696) q[2];
sx q[2];
rz(2.9048186) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8805485) q[1];
sx q[1];
rz(-2.4820584) q[1];
sx q[1];
rz(0.87415965) q[1];
rz(-pi) q[2];
rz(2.2569916) q[3];
sx q[3];
rz(-0.39456083) q[3];
sx q[3];
rz(1.7673499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87236658) q[2];
sx q[2];
rz(-1.1507611) q[2];
sx q[2];
rz(-1.9744138) q[2];
rz(-2.732318) q[3];
sx q[3];
rz(-0.27547488) q[3];
sx q[3];
rz(0.95612139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072227612) q[0];
sx q[0];
rz(-2.9874711) q[0];
sx q[0];
rz(-1.3519721) q[0];
rz(-1.6701472) q[1];
sx q[1];
rz(-2.2189326) q[1];
sx q[1];
rz(2.2460489) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9297139) q[0];
sx q[0];
rz(-1.3037056) q[0];
sx q[0];
rz(1.5729089) q[0];
x q[1];
rz(-2.9906435) q[2];
sx q[2];
rz(-0.47892919) q[2];
sx q[2];
rz(-0.72266173) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.99341678) q[1];
sx q[1];
rz(-1.9636781) q[1];
sx q[1];
rz(-0.82156374) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0655153) q[3];
sx q[3];
rz(-2.6722994) q[3];
sx q[3];
rz(2.9772907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6764549) q[2];
sx q[2];
rz(-0.82200161) q[2];
sx q[2];
rz(1.6949534) q[2];
rz(-2.1220574) q[3];
sx q[3];
rz(-1.8229702) q[3];
sx q[3];
rz(2.8442966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4296221) q[0];
sx q[0];
rz(-1.977704) q[0];
sx q[0];
rz(3.1373366) q[0];
rz(0.70107067) q[1];
sx q[1];
rz(-2.6316167) q[1];
sx q[1];
rz(0.027499011) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51873678) q[0];
sx q[0];
rz(-0.42930005) q[0];
sx q[0];
rz(1.1909498) q[0];
rz(-pi) q[1];
rz(-0.13032924) q[2];
sx q[2];
rz(-2.2918335) q[2];
sx q[2];
rz(0.70391612) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8054652) q[1];
sx q[1];
rz(-1.8847817) q[1];
sx q[1];
rz(2.9538395) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5309389) q[3];
sx q[3];
rz(-0.33646001) q[3];
sx q[3];
rz(0.57608673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6444401) q[2];
sx q[2];
rz(-1.7915598) q[2];
sx q[2];
rz(-2.7437317) q[2];
rz(-2.1286987) q[3];
sx q[3];
rz(-2.1161067) q[3];
sx q[3];
rz(2.7329172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6560646) q[0];
sx q[0];
rz(-1.2701472) q[0];
sx q[0];
rz(2.2960466) q[0];
rz(-0.0088508765) q[1];
sx q[1];
rz(-0.84819853) q[1];
sx q[1];
rz(1.1525851) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8542957) q[0];
sx q[0];
rz(-2.3919772) q[0];
sx q[0];
rz(-2.1867256) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0218618) q[2];
sx q[2];
rz(-1.3374723) q[2];
sx q[2];
rz(0.87750669) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.2680451) q[1];
sx q[1];
rz(-1.5427886) q[1];
sx q[1];
rz(0.59240474) q[1];
rz(2.939393) q[3];
sx q[3];
rz(-0.31980896) q[3];
sx q[3];
rz(-2.5068738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.74959603) q[2];
sx q[2];
rz(-2.1139202) q[2];
sx q[2];
rz(1.2006302) q[2];
rz(2.6642753) q[3];
sx q[3];
rz(-2.6715607) q[3];
sx q[3];
rz(-2.4376455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93477997) q[0];
sx q[0];
rz(-2.2581357) q[0];
sx q[0];
rz(-0.90428895) q[0];
rz(2.1705056) q[1];
sx q[1];
rz(-1.1957518) q[1];
sx q[1];
rz(-0.28360525) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87521919) q[0];
sx q[0];
rz(-3.0689256) q[0];
sx q[0];
rz(-0.22180264) q[0];
rz(-pi) q[1];
rz(-2.0118158) q[2];
sx q[2];
rz(-1.9649817) q[2];
sx q[2];
rz(-1.6617384) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1329258) q[1];
sx q[1];
rz(-2.2382793) q[1];
sx q[1];
rz(-0.26393601) q[1];
rz(2.6568703) q[3];
sx q[3];
rz(-2.1987763) q[3];
sx q[3];
rz(1.0768277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7859555) q[2];
sx q[2];
rz(-1.9060308) q[2];
sx q[2];
rz(-2.9034485) q[2];
rz(2.1613878) q[3];
sx q[3];
rz(-1.9832059) q[3];
sx q[3];
rz(0.75077209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
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
rz(0.86722803) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2326538) q[0];
sx q[0];
rz(-0.5687269) q[0];
sx q[0];
rz(-1.6816116) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32829702) q[2];
sx q[2];
rz(-0.5813501) q[2];
sx q[2];
rz(0.7996847) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7379511) q[1];
sx q[1];
rz(-0.80437901) q[1];
sx q[1];
rz(0.40108313) q[1];
x q[2];
rz(-1.8676742) q[3];
sx q[3];
rz(-1.6921077) q[3];
sx q[3];
rz(-0.47308009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0200218) q[2];
sx q[2];
rz(-2.3131804) q[2];
sx q[2];
rz(1.4309481) q[2];
rz(2.6089148) q[3];
sx q[3];
rz(-2.1055652) q[3];
sx q[3];
rz(-1.7238341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1581887) q[0];
sx q[0];
rz(-2.991365) q[0];
sx q[0];
rz(1.0839373) q[0];
rz(-3.0896507) q[1];
sx q[1];
rz(-0.40819326) q[1];
sx q[1];
rz(3.1165677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1442102) q[0];
sx q[0];
rz(-2.6591427) q[0];
sx q[0];
rz(0.061003322) q[0];
x q[1];
rz(0.45212169) q[2];
sx q[2];
rz(-2.2703301) q[2];
sx q[2];
rz(1.4750862) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96341771) q[1];
sx q[1];
rz(-1.1421848) q[1];
sx q[1];
rz(0.18842285) q[1];
rz(1.7134361) q[3];
sx q[3];
rz(-1.148292) q[3];
sx q[3];
rz(1.6137992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0673125) q[2];
sx q[2];
rz(-1.68579) q[2];
sx q[2];
rz(0.81706268) q[2];
rz(-0.95064154) q[3];
sx q[3];
rz(-1.0167511) q[3];
sx q[3];
rz(-1.954621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0364712) q[0];
sx q[0];
rz(-0.28577411) q[0];
sx q[0];
rz(-2.5407319) q[0];
rz(0.034424456) q[1];
sx q[1];
rz(-1.333678) q[1];
sx q[1];
rz(-1.6011802) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9118155) q[0];
sx q[0];
rz(-0.86965771) q[0];
sx q[0];
rz(2.3100353) q[0];
rz(2.1911084) q[2];
sx q[2];
rz(-0.99603876) q[2];
sx q[2];
rz(-1.7983914) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.34791495) q[1];
sx q[1];
rz(-1.131212) q[1];
sx q[1];
rz(-1.5949834) q[1];
x q[2];
rz(-1.1320587) q[3];
sx q[3];
rz(-2.0733548) q[3];
sx q[3];
rz(-2.9729321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2711266) q[2];
sx q[2];
rz(-1.8693962) q[2];
sx q[2];
rz(1.4938483) q[2];
rz(-2.3624524) q[3];
sx q[3];
rz(-1.4804163) q[3];
sx q[3];
rz(-0.65472764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7036024) q[0];
sx q[0];
rz(-2.0605189) q[0];
sx q[0];
rz(-2.3533452) q[0];
rz(1.2429271) q[1];
sx q[1];
rz(-1.5595167) q[1];
sx q[1];
rz(-16/(15*pi)) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96974715) q[0];
sx q[0];
rz(-1.7632428) q[0];
sx q[0];
rz(3.0234973) q[0];
rz(-pi) q[1];
rz(0.43717904) q[2];
sx q[2];
rz(-2.533982) q[2];
sx q[2];
rz(0.89113441) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.371468) q[1];
sx q[1];
rz(-1.9731908) q[1];
sx q[1];
rz(1.2312908) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3655534) q[3];
sx q[3];
rz(-2.0668012) q[3];
sx q[3];
rz(2.1797594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.01123151) q[2];
sx q[2];
rz(-1.9731382) q[2];
sx q[2];
rz(1.2569024) q[2];
rz(1.0907178) q[3];
sx q[3];
rz(-1.9382449) q[3];
sx q[3];
rz(-0.91712657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85482368) q[0];
sx q[0];
rz(-1.6463065) q[0];
sx q[0];
rz(-2.0538034) q[0];
rz(2.4913359) q[1];
sx q[1];
rz(-1.4573263) q[1];
sx q[1];
rz(0.021171721) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7268459) q[0];
sx q[0];
rz(-1.124589) q[0];
sx q[0];
rz(-1.3394974) q[0];
x q[1];
rz(-2.6708605) q[2];
sx q[2];
rz(-1.3952878) q[2];
sx q[2];
rz(-0.99711768) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.52439201) q[1];
sx q[1];
rz(-0.78689303) q[1];
sx q[1];
rz(1.8700429) q[1];
rz(3.0983119) q[3];
sx q[3];
rz(-2.1028404) q[3];
sx q[3];
rz(2.8325641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0040969) q[2];
sx q[2];
rz(-1.7986412) q[2];
sx q[2];
rz(0.98319483) q[2];
rz(-0.99299562) q[3];
sx q[3];
rz(-2.0296622) q[3];
sx q[3];
rz(2.907311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8773593) q[0];
sx q[0];
rz(-2.3491884) q[0];
sx q[0];
rz(1.0607251) q[0];
rz(-2.0296774) q[1];
sx q[1];
rz(-0.30525515) q[1];
sx q[1];
rz(0.12259604) q[1];
rz(-2.0531102) q[2];
sx q[2];
rz(-1.7822722) q[2];
sx q[2];
rz(1.7581802) q[2];
rz(3.1226782) q[3];
sx q[3];
rz(-2.6185642) q[3];
sx q[3];
rz(-2.0279283) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
