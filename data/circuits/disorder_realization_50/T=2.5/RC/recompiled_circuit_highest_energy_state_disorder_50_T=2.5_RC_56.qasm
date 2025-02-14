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
rz(1.7293575) q[0];
sx q[0];
rz(3.7339551) q[0];
sx q[0];
rz(10.629551) q[0];
rz(-2.0545948) q[1];
sx q[1];
rz(2.6462061) q[1];
sx q[1];
rz(10.612499) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76554322) q[0];
sx q[0];
rz(-0.57224579) q[0];
sx q[0];
rz(-2.8577198) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56894033) q[2];
sx q[2];
rz(-1.1747588) q[2];
sx q[2];
rz(-0.64096156) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0443327) q[1];
sx q[1];
rz(-2.9625434) q[1];
sx q[1];
rz(-1.8392842) q[1];
x q[2];
rz(1.2212444) q[3];
sx q[3];
rz(-2.7531525) q[3];
sx q[3];
rz(-0.74439186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54186934) q[2];
sx q[2];
rz(-2.5399127) q[2];
sx q[2];
rz(1.0428693) q[2];
rz(-0.045182191) q[3];
sx q[3];
rz(-2.9729645) q[3];
sx q[3];
rz(-1.6056304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.0586108) q[0];
sx q[0];
rz(-3.0284212) q[0];
sx q[0];
rz(-0.97852069) q[0];
rz(1.9784031) q[1];
sx q[1];
rz(-2.1470943) q[1];
sx q[1];
rz(1.8189583) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0019313) q[0];
sx q[0];
rz(-1.0884415) q[0];
sx q[0];
rz(0.84114065) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2873994) q[2];
sx q[2];
rz(-1.3958418) q[2];
sx q[2];
rz(-0.31084479) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6860477) q[1];
sx q[1];
rz(-1.3551222) q[1];
sx q[1];
rz(-2.2310745) q[1];
x q[2];
rz(1.421881) q[3];
sx q[3];
rz(-0.88203871) q[3];
sx q[3];
rz(-1.0721579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3375552) q[2];
sx q[2];
rz(-0.20122169) q[2];
sx q[2];
rz(0.03841722) q[2];
rz(-1.5086959) q[3];
sx q[3];
rz(-1.326606) q[3];
sx q[3];
rz(-1.647515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9134193) q[0];
sx q[0];
rz(-0.67263022) q[0];
sx q[0];
rz(-2.9056554) q[0];
rz(1.963223) q[1];
sx q[1];
rz(-1.447568) q[1];
sx q[1];
rz(0.49740121) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59374124) q[0];
sx q[0];
rz(-1.3914131) q[0];
sx q[0];
rz(-2.9602838) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3942986) q[2];
sx q[2];
rz(-1.9579534) q[2];
sx q[2];
rz(-0.55768572) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5798074) q[1];
sx q[1];
rz(-2.7473463) q[1];
sx q[1];
rz(0.057478776) q[1];
rz(-2.892832) q[3];
sx q[3];
rz(-1.8210871) q[3];
sx q[3];
rz(0.32969013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1824789) q[2];
sx q[2];
rz(-1.4664058) q[2];
sx q[2];
rz(1.2314931) q[2];
rz(-1.69151) q[3];
sx q[3];
rz(-1.3030038) q[3];
sx q[3];
rz(2.2836397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0602144) q[0];
sx q[0];
rz(-2.2980818) q[0];
sx q[0];
rz(-2.9486935) q[0];
rz(1.3554205) q[1];
sx q[1];
rz(-1.5602292) q[1];
sx q[1];
rz(-2.9706109) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9381704) q[0];
sx q[0];
rz(-1.6294384) q[0];
sx q[0];
rz(-1.0542271) q[0];
rz(-pi) q[1];
rz(0.73604446) q[2];
sx q[2];
rz(-1.9878329) q[2];
sx q[2];
rz(-1.2936178) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.65588412) q[1];
sx q[1];
rz(-7*pi/15) q[1];
sx q[1];
rz(1.1397362) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8754575) q[3];
sx q[3];
rz(-2.4551292) q[3];
sx q[3];
rz(1.0047508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96044797) q[2];
sx q[2];
rz(-1.9713216) q[2];
sx q[2];
rz(-2.81874) q[2];
rz(1.3747831) q[3];
sx q[3];
rz(-2.1026473) q[3];
sx q[3];
rz(2.3006181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5508995) q[0];
sx q[0];
rz(-1.0223848) q[0];
sx q[0];
rz(2.2160231) q[0];
rz(-0.12807056) q[1];
sx q[1];
rz(-1.6431199) q[1];
sx q[1];
rz(0.099460348) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2735705) q[0];
sx q[0];
rz(-1.5708459) q[0];
sx q[0];
rz(-0.8542866) q[0];
x q[1];
rz(-0.36014135) q[2];
sx q[2];
rz(-1.8966873) q[2];
sx q[2];
rz(-1.2831836) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0316098) q[1];
sx q[1];
rz(-1.3760202) q[1];
sx q[1];
rz(2.1094448) q[1];
rz(1.8811536) q[3];
sx q[3];
rz(-0.25729968) q[3];
sx q[3];
rz(-0.33724393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1818992) q[2];
sx q[2];
rz(-1.5965896) q[2];
sx q[2];
rz(-0.37364513) q[2];
rz(0.89961189) q[3];
sx q[3];
rz(-0.74266946) q[3];
sx q[3];
rz(-1.1348178) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9099092) q[0];
sx q[0];
rz(-2.9250356) q[0];
sx q[0];
rz(-2.5496971) q[0];
rz(-2.9369211) q[1];
sx q[1];
rz(-2.1494631) q[1];
sx q[1];
rz(3.0904904) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6986324) q[0];
sx q[0];
rz(-1.3784096) q[0];
sx q[0];
rz(-1.7257742) q[0];
rz(-pi) q[1];
rz(0.25801664) q[2];
sx q[2];
rz(-0.92890255) q[2];
sx q[2];
rz(-2.8947865) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6430059) q[1];
sx q[1];
rz(-0.47517794) q[1];
sx q[1];
rz(0.84334405) q[1];
rz(-pi) q[2];
rz(-2.3971812) q[3];
sx q[3];
rz(-0.6767513) q[3];
sx q[3];
rz(1.4686506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11703141) q[2];
sx q[2];
rz(-2.0270963) q[2];
sx q[2];
rz(-2.047211) q[2];
rz(-2.7545605) q[3];
sx q[3];
rz(-0.47457591) q[3];
sx q[3];
rz(-0.3064557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7759906) q[0];
sx q[0];
rz(-0.87823534) q[0];
sx q[0];
rz(-2.8940417) q[0];
rz(1.6869102) q[1];
sx q[1];
rz(-1.2431815) q[1];
sx q[1];
rz(0.036458485) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3943485) q[0];
sx q[0];
rz(-2.3009662) q[0];
sx q[0];
rz(-1.9565814) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22561947) q[2];
sx q[2];
rz(-2.4258412) q[2];
sx q[2];
rz(0.72474397) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.16462886) q[1];
sx q[1];
rz(-2.0348624) q[1];
sx q[1];
rz(1.0010757) q[1];
x q[2];
rz(3.0667449) q[3];
sx q[3];
rz(-1.6597431) q[3];
sx q[3];
rz(1.68881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.054691943) q[2];
sx q[2];
rz(-0.96618235) q[2];
sx q[2];
rz(1.7894233) q[2];
rz(-1.8933206) q[3];
sx q[3];
rz(-1.5636445) q[3];
sx q[3];
rz(-1.1917535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9098814) q[0];
sx q[0];
rz(-0.24869643) q[0];
sx q[0];
rz(1.6491718) q[0];
rz(0.17852783) q[1];
sx q[1];
rz(-1.2622204) q[1];
sx q[1];
rz(0.55007225) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15353841) q[0];
sx q[0];
rz(-3.1044391) q[0];
sx q[0];
rz(-2.1376993) q[0];
rz(1.3230611) q[2];
sx q[2];
rz(-1.6322517) q[2];
sx q[2];
rz(-2.2814192) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7272721) q[1];
sx q[1];
rz(-1.7309233) q[1];
sx q[1];
rz(-0.54267197) q[1];
x q[2];
rz(-0.86217238) q[3];
sx q[3];
rz(-2.5579427) q[3];
sx q[3];
rz(0.72008163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7940346) q[2];
sx q[2];
rz(-1.4057691) q[2];
sx q[2];
rz(-2.3883635) q[2];
rz(-2.8473162) q[3];
sx q[3];
rz(-3.0615276) q[3];
sx q[3];
rz(3.0729955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2733961) q[0];
sx q[0];
rz(-2.8559339) q[0];
sx q[0];
rz(-1.4519325) q[0];
rz(-2.946335) q[1];
sx q[1];
rz(-2.0611019) q[1];
sx q[1];
rz(1.6498227) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30376745) q[0];
sx q[0];
rz(-1.2653102) q[0];
sx q[0];
rz(-2.6432493) q[0];
rz(-1.4597462) q[2];
sx q[2];
rz(-0.80886474) q[2];
sx q[2];
rz(2.8946271) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6386913) q[1];
sx q[1];
rz(-1.4152621) q[1];
sx q[1];
rz(-3.0763094) q[1];
rz(-3.1141485) q[3];
sx q[3];
rz(-1.5641091) q[3];
sx q[3];
rz(-2.1777976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1008272) q[2];
sx q[2];
rz(-0.46697524) q[2];
sx q[2];
rz(-0.85774285) q[2];
rz(-1.6929251) q[3];
sx q[3];
rz(-1.492615) q[3];
sx q[3];
rz(-0.15374507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1345074) q[0];
sx q[0];
rz(-0.15468287) q[0];
sx q[0];
rz(2.3585228) q[0];
rz(-1.1231517) q[1];
sx q[1];
rz(-1.6936561) q[1];
sx q[1];
rz(-3.0812841) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9674112) q[0];
sx q[0];
rz(-0.36752146) q[0];
sx q[0];
rz(1.307748) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10580982) q[2];
sx q[2];
rz(-2.1544837) q[2];
sx q[2];
rz(0.90487827) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.26877182) q[1];
sx q[1];
rz(-0.67292456) q[1];
sx q[1];
rz(-0.10867837) q[1];
rz(-pi) q[2];
rz(-1.6844682) q[3];
sx q[3];
rz(-1.9316564) q[3];
sx q[3];
rz(-1.833503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9874838) q[2];
sx q[2];
rz(-2.2769112) q[2];
sx q[2];
rz(1.3101428) q[2];
rz(-1.3335258) q[3];
sx q[3];
rz(-1.3430877) q[3];
sx q[3];
rz(1.3069299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46351984) q[0];
sx q[0];
rz(-1.1165883) q[0];
sx q[0];
rz(-0.22620871) q[0];
rz(-0.75640596) q[1];
sx q[1];
rz(-0.49497985) q[1];
sx q[1];
rz(2.3351647) q[1];
rz(0.43037597) q[2];
sx q[2];
rz(-1.5793565) q[2];
sx q[2];
rz(-1.013365) q[2];
rz(1.9202833) q[3];
sx q[3];
rz(-1.3987473) q[3];
sx q[3];
rz(-1.7940298) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
